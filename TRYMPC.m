classdef TRYMPC < handle

   %%%%%%%%%%%%%%%%%%%%%%%%%%%% Read-Only:
   % Basic class-instance properties:
   properties(SetAccess = immutable)

      % Admin:
      Name (1,1) string
      ID (1,1) string
      
   end

   % Structs:
   properties(SetAccess = private)
      cas (1,1) struct % This is a struct that will contain all casadi-defined object (casadi symbolic variables, casadi functions, expressions, etc.) This struct can then be deleted, such that the class instance can be saved to a file
      dim (1,1) struct % contains the dimension of each symbolic array (state, algeb, input, param)
      names (1,1) struct % contains the names of all variables, both code names and display names, also contains the indices that the names appear in within the vectors
      integrator (1,1) struct % contains information about the integrator
      horizon (1,1) struct % contains information about the horizon (access to this is limited by its set-function)
      display (1,1) struct % contains information about the problem that is defined, such as number of optimization variables, number of constraints, etc.
      simulation (1,1) struct % contains information about the simulation; live initial condition, results, etc.
      % current (1,1) struct % contains information about the current optimization, relevant for mpc schemes. (current guess, settings, etc.)

      % "SQP_settings" contains the settigns of the SQP algorithm / QP subproblems:
      SQP_settings (1,1) struct = struct("tolerance_lagrangian",5,...
                                         "tolerance_equality",10^-6,...
                                         "tolerance_inequality",10^-6,...
                                         "QP_solver","quadprog_IP",...
                                         "quadprog_AS_options",optimoptions("quadprog","Algorithm","active-set","Display","none",MaxIterations=200,ConstraintTolerance = 10^-8),...
                                         "quadprog_IP_options",optimoptions("quadprog","Algorithm","interior-point-convex","Display","none",MaxIterations=200,ConstraintTolerance = 10^-8),...
                                         "quadprog_TR_options",optimoptions("quadprog","Algorithm","trust-region-reflective","Display","none",MaxIterations=200,ConstraintTolerance = 10^-8),...
                                         "max_N_iterations",100,...
                                         "merit_function",@(decision,lambda_next,mu_next) [],...
                                         "linesearch_method","backtracking",...
                                         "backtracking_rate",0.5,...
                                         "backtracking_min_step_size",10^-8);

      archive (1,1) struct = struct('optimizations',[],'simulations',[]) % contains simulation results (for all simulations that are performed)
   end




   %%%%%%%%%%%%%%%%%%%%%%%%%%%% Hidden
   % Miscellaneous
   properties(SetAccess = private,Hidden)
      var_types (1,:) string % lists the variable types present in the system (state, algebraic states, inputs, and parameters)
      var_types_notpar (1,:) string % same as var_types, but without "param"
      expr_types (1,:) string % types of expressions included in model (dynamics, algebraics)
      
      % flags:
      flag_has_algeb (1,1) logical = false;
      flag_has_input (1,1) logical = true;
      flag_has_param (1,1) logical = false;
   end
   properties(SetAccess = private,Hidden)
      cost_types (1,:) string % types of cost expressions in objective (stage_cost,Q,Z,R,dQ,dZ,dR)
      terminal_cost_types (1,:) string % types of terminal cost expressions in objective (Q_terminal,Z_terminal)

      flag_dynamics_defined (1,1) logical = false;
      flag_objective_defined (1,1) logical = false;
      flag_integrator_defied (1,1) logical = false;
      flag_stage_constraints_defined (1,1) logical = false;
      flag_terminal_constraints_defined (1,1) logical = false;
      flag_horizon_defied (1,1) logical = false;
      flag_problem_defined (1,1) logical = false;

      flag_has_bounds (1,1) logical = false;
      flag_has_terminal_bounds (1,1) logical = false;
      flag_has_inequality (1,1) logical = false;

      internal_sim (1,1) struct % contains information about the simulation, but that is not useful to the user. ("simulation" contains user-relevant info)
      internal_mpc (1,1) struct % contains info about the running optimization scheme (MPC), which are not relevant for user. % not sure yet: ("current" contains user-relevant info)
      SQP_info (1,1) struct % contains the logged information about an SQP solve routine temporarily, before being stored in an approperiate archive location. (mey vary for single solves and mpc runs)

      % Keep track of changing numerical values:
      flag_parameters_changed (1,1) logical = true;
      flag_quadratic_cost_changed (1,1) logical = true;
      flag_bounds_changed (1,1) logical = true;
      flag_terminal_bounds_changed (1,1) logical = true;
      flag_ref_changed (1,1) logical = true;
      flag_Dt_changed (1,1) logical = true;

      % Try to be efficient?:
      speed (1,1) logical = false;


      %%% Restoration:
      restore_cell (:,1) cell = {struct}; % for holding inputs to all "def_..." functions, so they can be regenerated. (does not work if external expressions are used. F.ex. manually providing the stage cost, or terminal bounds, etc...)

      % Horizon definition
      horizon_decider (1,1) string {mustBeMember(horizon_decider,["Dt","T"])} = "Dt";
   end
   properties(Constant,Hidden)
      cost_to_var = struct('Q',"state", ...
                          'dQ',"state", ...
                           'Z',"algeb",...
                          'dZ',"algeb",...
                           'R',"input",...
                          'dR',"input",...
                          'Q_terminal',"state",...
                          'Z_terminal',"algeb");


      allowable_DAE_solver = ["ode15s","ode23t","ode15i"];

   end
   % Immutables:
   properties(SetAccess = immutable,Hidden)
      
   end
  

   properties(Hidden,Dependent)
      % Evaluate objective value:
      current_objective
      
      % Evaluate constraints
      current_equality % ALL equality constraints
      current_inequality % ALL inequality constraints

      % Evaluate Jacobians (w.r.t. decision variables)
      current_J_objective
      current_J_equality
      current_J_inequality

      % Evaluate Hessian (w.r.t. decision variables)
      current_H_objective
      current_H_Lagrangian

      % Evaluate KKT matrix (Hessian of Lagrangian w.r.t. primal-dual vector)
      current_KKT_matrix

      % Keep track of any changes in numerical values:
      flag_numerical_values_changed (1,1) logical
   end




   %%%%%%%%%%%%%%%%%%%%%%%%%%%% Read and Write:
   % Properties that can be adjusted at any time:
   properties
      parameters (:,1) structor = structor("default_mix","separated") % numerical values for parameters, these values are queried whener the parameters are needed, thus chaning the parameters here, instantly changes the parameters of you system.
      quadratic_cost (1,1) struct = struct;
      bounds (1,1) struct = struct; % holds the upper and lower bounds of each individual variable. Ex: "C.bounds.upper.th = 2"
      terminal_bounds (1,1) struct = struct; % holds the upper and lower bounds of each terminal variable. Ex: "C.terminal_bounds.upper.ux = 0.2"
      ref (1,1) struct % shold contain a reference trajectory for states, algebs, and inputs
      initial_state (1,1) structor = structor;

      plotting (1,1) struct % contains informations about desired grafic settings for the various parameters
   end



   %%%% Defining the class, model and problem:
   methods(Access = public)

      %%%%%%%%%%%%%%%%%%%%%%%%%% CONSTRUCTOR:
      function C = TRYMPC(class_instance_name,names)
         arguments
            class_instance_name (1,1) string
            names.state (:,1) string
            names.algeb (:,1) string
            names.input (:,1) string
            names.param (:,1) string
         end

         fprintf(['Created a TRYMPC instance: "', char(class_instance_name),'"... ']);
         create_time = tic;


         % Apply name:
         C.Name = class_instance_name;

         % Create unique ID for class intance:
         C.ID = C.generate_id; % probably unique...

         % Generate variables
         C.def_instance(names)

         % Set default initial state
         for name = C.names.code.state'
            C.initial_state.str.(name) = 0;
         end
         
         create_time = toc(create_time);
         disp(['done.  ',sec2str(create_time)]);
         disp( '            -- Use "TRYMPC.help" for further instructions.')
      end
      function def_instance(C,names)
         % Prepare restoration
         C.restore_cell{end}.def_instance = {names};


         % Characterize types:
         C.var_types = string(fieldnames(names))';
         C.var_types_notpar = C.var_types(C.var_types ~= "param");

         if C.var_types(1) ~= "state"
            C.usererror('must provide a list of state names.')
         end
         if ~any(C.var_types == "input")
            C.flag_has_input = false;
            % C.usererror('must provide a list of input names.')
         end

         C.expr_types = "dynamics";
         if isfield(names,"algeb")
            C.flag_has_algeb = true;
            C.expr_types(end+1) = "algebraics";
         end


         % Define variables:
         C.cas.var = struct;
         for type = C.var_types
            C.def_vars(type,names.(type))
         end

         % Define basic plotting-properteis:
         counter = 1;
         for type = C.var_types
            for name = names.(type)'
               C.plotting.display_names.(type).(name) = name;
               C.plotting.color.(type).(name) = GetColorCode(counter);
               C.plotting.linestyle.(type).(name) = '-';
               C.plotting.linewidth.(type).(name) = 1.5;
               counter = counter + 1;
            end
         end


         % Create parameter structor:
         if isfield(names,"param")
            C.flag_has_param = true;
            for name = C.names.code.param'
               C.parameters.str.(name) = nan;
            end
         end
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%% VARIABLES:
      function def_vars(C,type,names)
         arguments
            C
            type (1,1) string {mustBeMember(type,["state","algeb","input","param"])}
            names (:,1) string
         end

         C.dim.(type) = length(names);

         % Create all variables and add to a vector
         C.cas.var.(type) = casadi.SX.sym(char(names(1)));
         if size(names,1) > 1
            for name = names(2:end)'
               C.cas.var.(type) = [C.cas.var.(type); casadi.SX.sym(char(name))];
            end
         end

         if type ~= "param"
            C.cas.ref.(type) = casadi.SX.sym([char(names(1)),'_ref']);
            C.cas.next.(type) = casadi.SX.sym([char(names(1)),'_next']);
            C.cas.init.(type) = casadi.SX.sym([char(names(1)),'_init']);
            if size(names,1) > 1
               for name = names(2:end)'
                  C.cas.ref.(type) = [C.cas.ref.(type); casadi.SX.sym([char(name),'_ref'])];
                  C.cas.next.(type) = [C.cas.next.(type); casadi.SX.sym([char(name),'_next'])];
                  C.cas.init.(type) = [C.cas.init.(type); casadi.SX.sym([char(name),'_init'])];
               end
            end
         end

         % Add all new variables to separate fields for easy access
         for i = 1:C.dim.(type)
            C.cas.(type).(names(i)) = C.cas.var.(type)(i);
            C.names.ind.(type).(names(i)) = i;
         end
         C.names.code.(type) = names;

      end


      %%%%%%%%%%%%%%%%%%%%%%%%%% MODEL:
      function def_dynamics(C,expr_dynamics,expr_algebraics)
         arguments
            C
            expr_dynamics (:,1) function_handle % an anonymous function (dyn = @(s,a,i,p) *...*) that generates a CasASi symbolic (SX) expression for the system dynamics (it will be called by dyn(C.cas.state,C.cas.algeb,C.cas.input,C.cas.param))
            expr_algebraics (:,1) function_handle = @(s,a,i,p) []; % like dynamics, should generate Casadi symbolic expression for the system algebraic equations
         end
         % Prepare restoration
         C.restore_cell{end}.def_dynamics = {expr_dynamics,expr_algebraics};

         fprintf('Defining dynamics...    ');
         def_time = tic;

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % immediately convert to casadi.SX expressions:
         if C.flag_has_algeb
            a = C.cas.algeb;
         else
            a = struct;
         end
         if C.flag_has_input
            i = C.cas.input;
         else
            i = struct;
         end
         if C.flag_has_param
            p = C.cas.param;
         else
            p = struct;
         end
         expr_dynamics = expr_dynamics(C.cas.state,a,i,p);
         expr_algebraics = expr_algebraics(C.cas.state,a,i,p);
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         if string(class(expr_dynamics)) ~= "casadi.SX"
            error('USER ERROR: dynamics function_handle must take (state,algeb,input,param) structs as arguments, and output a "casadi.SX" expression.')
         elseif (size(expr_dynamics,1) ~= C.dim.state)
            error('USER ERROR: dynamics vector not same size as state vector')
         end

         if C.flag_has_algeb
            if string(class(expr_algebraics)) ~= "casadi.SX"
               error('USER ERROR: algebraic function_handle must take (state,algeb,input,param) structs as arguments, and output a "casadi.SX" expression. (Algebraics: g(x,z,u,p) = 0).')
            end
         end


         % Create CasADi Function for dynamics:
         Arg = C.cas.var;
         in_fields = fieldnames(Arg);

         expr_master = {expr_dynamics, expr_algebraics};
         expr_type_master = C.expr_types;
         while ~isempty(expr_type_master)
            expr_type = expr_type_master(1);
            expr = expr_master{1};
            
            % Add expression to cas:
            C.cas.(expr_type).expr = expr;

            Arg.out = expr;
            C.cas.(expr_type).F = casadi.Function(['F_',char(expr_type)],Arg,in_fields,{'out'});

            % Create Jacobian (wrt states):
            for type = C.var_types
               C.cas.(expr_type).jacobian.(type).expr = jacobian(expr,C.cas.var.(type));
               Arg.out = C.cas.(expr_type).jacobian.(type).expr;
               C.cas.(expr_type).jacobian.(type).F    = casadi.Function(['F_',char(expr_type),'_',char(type)],Arg,in_fields,{'out'});
            end

            expr_type_master(1) = [];
            expr_master(1) = [];
         end

         C.flag_dynamics_defined = true;
         def_time = toc(def_time);
         disp(['done.  ',sec2str(def_time)])
      end


      %%%%%%%%%%%%%%%%%%%%%%%%%% COST:
      function def_objective(C,options)
         arguments
            C

            options.quadratic (1,:) string {mustBeMember(options.quadratic,["Q","Z","R","dQ","dZ","dR","Q_terminal","Z_terminal"])}
            
            % %%% All the terms below are added to the total objective, even
            % %%% when both stage_cost and Q, R are given.
            % options.Q        (1,:) double  % if a quadratic cost is used, simply provide weights on each state here
            % options.Z        (1,:) double  % the same type of cost can be added for algebraic variables
            % options.R        (1,:) double  % ... and input cost here
            % 
            % options.dQ       (1,:) double  % if a cost on the difference in the state variable between stages
            % options.dZ       (1,:) double  % algebraic difference
            % options.dR       (1,:) double  % difference in input

            options.stage_cost (1,1) casadi.SX % if an expression is provided manually
            options.terminal_cost (1,1) casadi.SX % if an expression is provided manually
         end
         % Prepare restoration
         fields = fieldnames(options);
         values = struct2cell(options);
         nameValuePairs = [fields'; values'];
         C.restore_cell{end}.def_objective = nameValuePairs;

         fprintf('Defining objective...   ');
         def_time = tic;


         C.cost_types = [];
         C.terminal_cost_types = [];
         if isfield(options,'quadratic')

            % Create casadi.SX variables to represent the various weights
            for weight_type = options.quadratic
                  var_type = C.cost_to_var.(weight_type);
                  if ismember(var_type,["algeb","input"]) && ~C.(['flag_has_',char(var_type)])
                     warning(['Warning: you are trying to add ',char(weight_type),'-type stage-cost terms, but no ',char(var_type),'-variables have been defined. ignoring...'])
                  else
                     if ~ismember(weight_type,["Q_terminal","Z_terminal"])
                        C.cost_types = [C.cost_types weight_type];
                     else
                        C.terminal_cost_types = [C.terminal_cost_types weight_type];
                     end
                     C.cas.objective.weights.(weight_type) = [];
                     for name = C.names.code.(var_type)'
                        C.cas.objective.weights.(weight_type) = [C.cas.objective.weights.(weight_type);  casadi.SX.sym(['weight_',char(weight_type),'_',char(name)]) ];
                     end
                  end

            end

         end
         if isfield(options,'stage_cost')
            C.cost_types = [C.cost_types "stage_cost"];
         end
         if isfield(options,'terminal_cost')
            C.terminal_cost_types = [C.cost_types "terminal_cost"];
         end
         if isempty(C.cost_types)
            C.usererror('no objective types were selectd.')
         end

         

         %%%% Create a casadi function to represent stage wise objective terms:
         % First create temporary versions that don't include the reference
         % (stage difference does not need reference)
         for weight_type = [C.cost_types C.terminal_cost_types]
            Arg = struct;
            switch weight_type
               case {"Q", "Z", "R", "Q_terminal", "Z_terminal"}
                  var_type = C.cost_to_var.(weight_type);

                  % Input fields:

                  Arg.(var_type) = C.cas.var.(var_type);
                  Arg.([char(var_type),'_ref']) = C.cas.ref.(var_type);
                  Arg.(weight_type) = C.cas.objective.weights.(weight_type);
                  in_fields = fieldnames(Arg);

                  % % Add to total argument struct:
                  % Arg_total.(var_type) = C.cas.var.(var_type);
                  % Arg_total.([char(var_type),'_ref']) = C.cas.ref.(var_type);
                  % Arg_total.(weight_type) = C.cas.objective.weights.(weight_type);

                  % output fields:
                  Arg.out = casadi.SX.zeros;
                  for name = C.names.code.(var_type)'
                     ind = C.names.ind.(var_type).(name);
                     Arg.out = Arg.out + C.cas.objective.weights.(weight_type)(ind)*(C.cas.(var_type).(name) - C.cas.ref.(var_type)(ind))^2;
                  end

               case {"dQ", "dZ", "dR"}
                  var_type = C.cost_to_var.(weight_type);

                  % Input fields:
                  Arg.(var_type) = C.cas.var.(var_type);
                  Arg.(weight_type) = C.cas.objective.weights.(weight_type);
                  Arg.([char(var_type),'_next']) = C.cas.next.(var_type);
                  in_fields = fieldnames(Arg);

                  % % Add to total argument struct:
                  % Arg_total.(var_type) = C.cas.var.(var_type); % might be done twice becuuse of non-difference version of the same var_type, but that's fine.
                  % Arg_total.(weight_type) = C.cas.objective.weights.(weight_type);
                  % Arg_total.([char(var_type),'_next']) = C.cas.next.(var_type);

                  % output fields:
                  Arg.out = casadi.SX.zeros;
                  for name = C.names.code.(var_type)'
                     ind = C.names.ind.(var_type).(name);
                     Arg.out = Arg.out + C.cas.objective.weights.(weight_type)(ind)*(C.cas.(var_type).(name) - C.cas.next.(var_type)(ind))^2;
                  end


               case {"stage_cost","terminal_cost"}

                  for var_type = C.var_types
                     Arg.(var_type) = C.cas.var.(var_type);
                     if var_type ~= "param"
                        Arg.([char(var_type),'_ref'])  = C.cas.ref.(var_type);
                        if weight_type == "stage_cost"
                           Arg.([char(var_type),'_next']) = C.cas.next.(var_type);
                        end
                     end
                  end
                  in_fields = fieldnames(Arg);
                  Arg.out = options.(weight_type);

            end
            C.cas.objective.cost.(weight_type).F = casadi.Function(['F_objective_',char(weight_type)],Arg,in_fields,{'out'});
            C.cas.objective.cost.(weight_type).in_fields = in_fields;
         end


         %%%% TOTAL STAGE COST:

         % Get all possible variable types and weights:
         Arg = C.cas.objective.weights;
         for var_type = C.var_types
            Arg.(var_type) = C.cas.var.(var_type);
            if var_type ~= "param"
               Arg.([char(var_type),'_ref'])  = C.cas.ref.(var_type);
               Arg.([char(var_type),'_next']) = C.cas.next.(var_type);
            end
         end
         in_fields = fieldnames(Arg);

         out = casadi.SX.zeros;
         for cost_type = C.cost_types
            arg = struct;
            for field = string(C.cas.objective.cost.(cost_type).in_fields)'
               arg.(field) = Arg.(field);
            end
            out = out + C.cas.objective.cost.(cost_type).F.call(arg).out;
         end
         Arg.out = out;

         C.cas.objective.cost.stage_total.F = casadi.Function('F_stage_objective',Arg,in_fields,{'out'});
         C.cas.objective.cost.stage_total.in_fields = in_fields;


         %%%% TERMINAL COST:

         % Get all possible variable types and weights:
         Arg = C.cas.objective.weights;
         for var_type = C.var_types
            Arg.(var_type) = C.cas.var.(var_type);
            if var_type ~= "param"
               Arg.([char(var_type),'_ref'])  = C.cas.ref.(var_type);
            end
         end
         in_fields = fieldnames(Arg);

         out = casadi.SX.zeros;
         for cost_type = C.terminal_cost_types
            arg = struct;
            for field = string(C.cas.objective.cost.(cost_type).in_fields)'
               arg.(field) = Arg.(field);
            end
            out = out + C.cas.objective.cost.(cost_type).F.call(arg).out;
         end
         Arg.out = out;

         C.cas.objective.cost.terminal_total.F = casadi.Function('F_terminal_objective',Arg,in_fields,{'out'});
         C.cas.objective.cost.terminal_total.in_fields = in_fields;



         C.flag_objective_defined = true;
         def_time = toc(def_time);
         disp(['done.  ',sec2str(def_time)])
      end
   

      %%%%%%%%%%%%%%%%%%%%%%%%%% INTEGRATOR:
      function def_integrator(C,method,options)
         arguments
            C
            method (1,1) string {mustBeMember(method,["Explicit Euler","Implicit Euler","ERK4","ERK4 (simultaneous)","Implicit Midpoint","Crank-Nicolson (Implicit)","IRK4 (L-stable)","Gauss-Legendre (4. order)","Gauss-Legendre (6. order)","collocation","external","explicit butcher tableau","implicit butcher tableau"])}

            % General settings:
            options.n_increments (1,1) double {mustBePositive,mustBeInteger} = 1;

            % Collocaiton Specific settigns:
            options.collocation_polynomial_order (1,1) double {mustBePositive,mustBeInteger} = 2;
            options.collocation_polynomial_type (1,1) string {mustBeMember(options.collocation_polynomial_type,["legendre","radau"])} = "legendre";

            options.bucher_tableau_b (1,:) double {mustBeReal}
            options.bucher_tableau_A (:,:) double {mustBeReal}
         end
         % Prepare restoration
         fields = fieldnames(options);
         values = struct2cell(options);
         nameValuePairs = [fields'; values'];
         C.restore_cell{end}.def_integrator = [{method},reshape(nameValuePairs,1,[])];

         fprintf('Defining integrator...  ');
         def_time = tic;
         if ~C.flag_dynamics_defined
            C.usererror('the dynamics must be defined before one can define an integrator. Use: ".def_dynamics" to define dynamics.')
         end

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Clear the integrator in case another integrator has been defined
         % already.
         C.cas.integrator = [];
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         C.integrator = struct("order",[],"shooting_method",[],"integrator",[],"n_increments",[],"has_aux",false);
         C.integrator.integrator = method;
         C.integrator.n_increments = options.n_increments;
         
         C.cas.integrator.Dt = casadi.SX.sym('internal_Dt');

         % define integrator:
         switch method
            case "Explicit Euler"
               C.integrator.order = 1;
               C.integrator.shooting_method = "multiple shooting";
               
               % Butcher Tableau:
               C.integrator.BT_b = 1;
               C.integrator.BT_A = 0;

               C.ERK_builder

            case "Implicit Euler"
               C.integrator.has_aux = true;

               C.integrator.order = 1;
               C.integrator.shooting_method = "simultaneous shooting";

               % Butcher Tableau:
               C.integrator.BT_b = 1;
               C.integrator.BT_A = 1;

               C.IRK_builder

            case "ERK4"
               C.integrator.order = 4;
               C.integrator.shooting_method = "multiple shooting";

               % Butcher Tableau:
               C.integrator.BT_b = [1 2 2 1]/6;
               C.integrator.BT_A = zeros(4);
               C.integrator.BT_A(2,1) = 1/2;
               C.integrator.BT_A(3,2) = 1/2;
               C.integrator.BT_A(4,3) = 1;

               C.ERK_builder

            case "ERK4 (simultaneous)"
               C.integrator.order = 4;
               C.integrator.shooting_method = "simultaneous shooting";
               C.integrator.has_aux = true;

               % Butcher Tableau:
               C.integrator.BT_b = [1 2 2 1]/6;
               C.integrator.BT_A = zeros(4);
               C.integrator.BT_A(2,1) = 1/2;
               C.integrator.BT_A(3,2) = 1/2;
               C.integrator.BT_A(4,3) = 1;

               C.IRK_builder

            case "Implicit Midpoint"
               C.integrator.has_aux = true;

               C.integrator.order = 2;
               C.integrator.shooting_method = "simultaneous shooting";

               % Butcher Tableau:
               C.integrator.BT_b = 1;
               C.integrator.BT_A = 1/2;

               C.IRK_builder

            case "Crank-Nicolson (Implicit)"
               C.integrator.has_aux = true;

               C.integrator.order = 2;
               C.integrator.shooting_method = "simultaneous shooting";

               % Butcher Tableau:
               C.integrator.BT_b = [1 1]/2;
               C.integrator.BT_A = [0 0; 1 1]/2;

               C.IRK_builder

            case "Gauss-Legendre (4. order)"
               C.integrator.has_aux = true;

               C.integrator.order = 4;
               C.integrator.shooting_method = "simultaneous shooting";

               % Butcher Tableau:
               C.integrator.BT_b = [1 1]/2;
               C.integrator.BT_A = [   1/4         1/4-sqrt(3)/6 ;
                  1/4+sqrt(3)/6      1/4       ];

               C.IRK_builder

            case "Gauss-Legendre (6. order)"
               C.integrator.has_aux = true;

               C.integrator.order = 6;
               C.integrator.shooting_method = "simultaneous shooting";

               % Butcher Tableau:
               C.integrator.BT_b = [5/18 4/9 5/18];
               C.integrator.BT_A = [5/36              2/9-sqrt(15)/15   5/36-sqrt(15)/30;
                                    5/36+sqrt(15)/24      2/9           5/36-sqrt(15)/24;
                                    5/36+sqrt(15)/30  2/9+sqrt(15)/15   5/36];

               C.IRK_builder

            case "IRK4 (L-stable)"
               C.integrator.order = 3;
               C.integrator.shooting_method = "simultaneous shooting";
               C.integrator.has_aux = true;

               % Butcher Tableau:
               C.integrator.BT_b = [3 -3 1 1]/2;
               C.integrator.BT_A = [ 1   0   0   0  ;
                                    1/3  1   0   0  ;
                                    -1   1   1   0  ;
                                     3  -3   1   1  ] /2;
               % C.integrator.BT_A(1,1) =  1/2;
               % C.integrator.BT_A(2,1) =  1/6;
               % C.integrator.BT_A(3,1) = -1/2;
               % C.integrator.BT_A(4,1) =  3/2;
               % C.integrator.BT_A(2,2) =  1/2;
               % C.integrator.BT_A(3,2) =  1/2;
               % C.integrator.BT_A(4,2) = -3/2;
               % C.integrator.BT_A(3,3) =  1/2;
               % C.integrator.BT_A(4,3) =  1/2;
               % C.integrator.BT_A(4,4) =  1/2;

               C.IRK_builder

            case "explicit butcher tableau"
               if ~isfield(options,'bucher_tableau_b')
                  C.usererror('the b-vector of a butcher tableau must be defined to use a custom butcher tableau. Try f.ex. "ERK4" for a predefined RK method.')
               end
               if ~isfield(options,'bucher_tableau_A')
                  C.usererror('the A-matrix of a butcher tableau must be defined to use a custom butcher tableau. Try f.ex. "ERK4" for a predefined RK method.')
               end

               C.integrator.order = "unknown";
               C.integrator.shooting_method = "multiple shooting";

               % Butcher Tableau:
               C.integrator.BT_b = options.bucher_tableau_b;
               C.integrator.BT_A = options.bucher_tableau_A;

               C.ERK_builder

            case "implicit butcher tableau"
               if ~isfield(options,'bucher_tableau_b')
                  C.usererror('the b-vector of a butcher tableau must be defined to use a custom butcher tableau. Try f.ex. "IRK4" for a predefined RK method.')
               end
               if ~isfield(options,'bucher_tableau_A')
                  C.usererror('the A-matrix of a butcher tableau must be defined to use a custom butcher tableau. Try f.ex. "IRK4" for a predefined RK method.')
               end

               C.integrator.order = "unknown";
               C.integrator.shooting_method = "simultaneous shooting";
               C.integrator.has_aux = true;

               % Butcher Tableau:
               C.integrator.BT_b = options.bucher_tableau_b;
               C.integrator.BT_A = options.bucher_tableau_A;

               C.IRK_builder

            case "collocation"

               % save collocation settings:
               C.integrator.collocation.d = options.collocation_polynomial_order;
               C.integrator.collocation.polynomial_type = options.collocation_polynomial_type;
               
               C.integrator.shooting_method = "direct collocation";
               C.integrator.has_aux = true;
               
               switch C.integrator.collocation.polynomial_type
                  case "legendre"
                     C.integrator.order = C.integrator.collocation.d*2;
                  case "radau"
                     C.integrator.order = C.integrator.collocation.d*2 - 1;
                     error('DEVELOPER ERROR: ops, I am unsure if I have built the collocation scheme specifically for Legendre, of if simply choosing Radau point will procude the correct Radau collocation shceme. - Trym')
                  otherwise
                     error('DEVELOPER ERROR: an invalid colloocation polynomial type was selected')
               end
               C.collocation_builder
         end
         

         C.flag_integrator_defied = true;
         def_time = toc(def_time);
         disp(['done.  ',sec2str(def_time)])
      end
   

      %%%%%%%%%%%%%%%%%%%%%%%%%% STAGE CONSTRAINTS:
      function def_stage_constraints(C,options)
         arguments
            C
         end
         arguments

            % include the names of variables that have bounds, the bound values are then given via the "C.bounds" struct
            options.lower_bounds (1,:) string 
            options.upper_bounds (1,:) string
            
            % General stage constraints:
            options.equality   (1,1) struct % general stage-wise equality constraints given as anonymous functions: @(s,a,i,p) (g() = 0)
            options.inequality (1,1) struct % general stage-wise inequlity constraints given as anonymous functions: @(s,a,i,p) (h() >= 0)
         end

         % First verify that the equality and inequality options are
         % structs of function handles
         for type = ["equality" "inequality"]
            if isfield(options,type)
               for name = string(fieldnames(options.(type)))'
                  if ~isa(options.(type).(name),'function_handle')
                     C.usererror('The equality and inequality optional arguments must be structs with fields containing only funciton handles on the form: @(s,a,i,p) (...).')
                  end
               end
            end
         end

         % Prepare restoration
         fields = fieldnames(options);
         values = struct2cell(options);
         nameValuePairs = [fields'; values'];
         C.restore_cell{end}.def_stage_constraints = nameValuePairs;

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % immediately convert to casadi.SX expressions:
         if C.flag_has_algeb
            a = C.cas.algeb;
         else
            a = struct;
         end
         if C.flag_has_input
            i = C.cas.input;
         else
            i = struct;
         end
         if C.flag_has_param
            p = C.cas.param;
         else
            p = struct;
         end
         if isfield(options,'equality')
            for name = string(fieldnames(options.equality))'
               options.equality.(name)   = options.equality.(name)(C.cas.state,a,i,p);
            end
         end
         if isfield(options,'inequality')
            for name = string(fieldnames(options.inequality))'
               options.inequality.(name)   = options.inequality.(name)(C.cas.state,a,i,p);
            end
         end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         fprintf('Defining stage constraints... ');
         def_time = tic;
         if C.flag_horizon_defied
            C.usererror('The stage-wise constraints must be defined BEFORE the horizon (.def_horizon()).')
         end

         Arg = C.cas.var;
         infields = fieldnames(Arg);

         if isfield(options,'equality')
            equality = options.equality;
            for name = string(fieldnames(equality))'
               if isa(equality.(name),'casadi.SX')
                  if size(equality.(name),2) ~= 1
                     C.usererror('the equality and inequality constraints must be column-vector expressions. That is: size(g()) is (:,1).')
                  end
                  C.cas.constraints.equality.(name).expr = equality.(name);
                  Arg.out = equality.(name);
                  C.cas.constraints.equality.(name).F = casadi.Function(['F_equality_',char(name)],Arg,infields,{'out'});
               else
                  C.usererror(['equalities must be entered as structs with fields that are "casadi.SX" expressions (g() = 0), consisting of variables defined in myTRYMPC.cas.state / myTRYMPC.cas.algeb / myTRYMPC.cas.input / myTRYMPC.cas.param ... you have a "',class(equality.(name)),'"'])
               end
            end
         end
         if isfield(options,'inequality')
            C.flag_has_inequality = true;
            inequality = options.inequality;
            for name = string(fieldnames(inequality))'
               if isa(inequality.(name),'casadi.SX')
                  if size(inequality.(name),2) ~= 1
                     C.usererror('the equality and inequality constraints must be column-vector expressions. That is: size(h()) is (:,1).')
                  end
                  C.cas.constraints.inequality.(name).expr = inequality.(name);
                  Arg.out = inequality.(name);
                  C.cas.constraints.inequality.(name).F = casadi.Function(['F_inequality_',char(name)],Arg,infields,{'out'});
               else
                  C.usererror(['inequalities must be entered as structs with fields that are "casadi.SX" expressions (h() >= 0), consisting of variables defined in myTRYMPC.cas.state / myTRYMPC.cas.algeb / myTRYMPC.cas.input / myTRYMPC.cas.param ... you have a "',class(inequality.(name)),'"'])
               end
            end
         end

         % Form lower bounds:
         expr = [];
         for bound_type = ["lower", "upper"]
            if isfield(options,[char(bound_type),'_bounds'])
               C.flag_has_bounds = true;
               for name = options.([char(bound_type),'_bounds'])
                  invalid_name = true;
                  for type = C.var_types_notpar
                     if ismember(name,C.names.code.(type))
                        invalid_name = false;
                        C.cas.bounds.(bound_type).(name) = casadi.SX.sym([char(bound_type),'_bound_',char(name)]);
                        if bound_type == "lower"
                           expr = [expr; -C.cas.bounds.(bound_type).(name) + C.cas.(type).(name)]; %#ok<AGROW>
                        elseif bound_type == "upper"
                           expr = [expr;  C.cas.bounds.(bound_type).(name) - C.cas.(type).(name)]; %#ok<AGROW>
                        end
                     end
                  end
                  if invalid_name
                     warning(['WARNING: the name provided: "',char(name),'" is not a variable name. No bound is created for this name.'])
                  end
               end
            end
         end
         if C.flag_has_bounds
            C.flag_has_inequality = true;
            C.cas.constraints.bounds.expr = expr;
            Arg.out = C.cas.constraints.bounds.expr;
            C.cas.constraints.bounds.F = casadi.Function('variable_bounds',Arg,infields,{'out'});
         end

         C.flag_stage_constraints_defined = true;
         def_time = toc(def_time);
         disp(['done.  ',sec2str(def_time)])
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%% TERMINAL CONSTRAINTS:
      function def_terminal_constraint(C,options)
         arguments
            C

            % include the names of variables that have bounds, the bound values are then given via the "C.bounds" struct
            options.lower_bounds (1,:) string 
            options.upper_bounds (1,:) string
            
            % General stage constraints:
            options.equality   (1,1) struct % general stage-wise equality constraints (g() = 0)
            options.inequality (1,1) struct % a general stage-wise inequlity constraints (h() >= 0)
         end

         % First verify that the equality and inequality options are
         % structs of function handles
         for type = ["equality" "inequality"]
            if isfield(options,type)
               for name = string(fieldnames(options.(type)))'
                  if ~isa(options.(type).(name),'function_handle')
                     C.usererror('The equality and inequality optional arguments must be structs with fields containing only funciton handles on the form: @(s,a,i,p) (...).')
                  end
               end
            end
         end

         % ' In tutorial 2 - I try adding terminal equality constraints, but they dont seem to be fully satisfied, and the C.cas.horizon.constraints.equality.str does not have a "stage_50", which is where these constraints should be found....
         % Prepare restoration
         fields = fieldnames(options);
         values = struct2cell(options);
         nameValuePairs = [fields'; values'];
         C.restore_cell{end}.def_terminal_constraints = nameValuePairs;

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % immediately convert to casadi.SX expressions:
         if C.flag_has_algeb
            a = C.cas.algeb;
         else
            a = struct;
         end
         if C.flag_has_input
            i = C.cas.input;
         else
            i = struct;
         end
         if C.flag_has_param
            p = C.cas.param;
         else
            p = struct;
         end
         if isfield(options,'equality')
            for name = string(fieldnames(options.equality))'
               options.equality.(name)   = options.equality.(name)(C.cas.state,a,i,p);
            end
         end
         if isfield(options,'inequality')
            for name = string(fieldnames(options.inequality))'
               options.inequality.(name)   = options.inequality.(name)(C.cas.state,a,i,p);
            end
         end
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


         fprintf('Defining terminal constraints... ');
         def_time = tic;
         if C.flag_horizon_defied
            C.usererror('The terminal constraints must be defined BEFORE the horizon (.def_horizon()).')
         end

         Arg = C.cas.var;
         infields = fieldnames(Arg);

         if isfield(options,'equality')
            equality = options.equality;
            for name = string(fieldnames(equality))'
               if isa(equality.(name),'casadi.SX')
                  if size(equality.(name),2) ~= 1
                     C.usererror('the equality and inequality constraints must be column-vector expressions. That is: size(g()) is (:,1).')
                  end
                  C.cas.constraints.terminal.equality.(name).expr = equality.(name);
                  Arg.out = equality.(name);
                  C.cas.constraints.terminal.equality.(name).F = casadi.Function(['F_terminal_equality_',char(name)],Arg,infields,{'out'});
               else
                  C.usererror(['terminal equalities must be entered as structs with fields that are "casadi.SX" expressions (g() = 0), consisting of variables defined in myTRYMPC.cas.state / myTRYMPC.cas.algeb / myTRYMPC.cas.input / myTRYMPC.cas.param ... you have a "',class(equality.(name)),'"'])
               end
            end
         end
         if isfield(options,'inequality')
            C.flag_has_inequality = true;
            inequality = options.inequality;
            for name = string(fieldnames(inequality))'
               if isa(inequality.(name),'casadi.SX')
                  if size(inequality.(name),2) ~= 1
                     C.usererror('the equality and inequality constraints must be column-vector expressions. That is: size(h()) is (:,1).')
                  end
                  C.cas.constraints.terminal.inequality.(name).expr = inequality.(name);
                  Arg.out = inequality.(name);
                  C.cas.constraints.terminal.inequality.(name).F = casadi.Function(['F_terminal_inequality_',char(name)],Arg,infields,{'out'});
               else
                  C.usererror(['terminal inequalities must be entered as structs with fields that are "casadi.SX" expressions (h() >= 0), consisting of variables defined in myTRYMPC.cas.state / myTRYMPC.cas.algeb / myTRYMPC.cas.input / myTRYMPC.cas.param ... you have a "',class(inequality.(name)),'"'])
               end
            end
         end

         % Form lower bounds:
         expr = [];
         for bound_type = ["lower", "upper"]
            if isfield(options,[char(bound_type),'_bounds'])
               C.flag_has_terminal_bounds = true;
               for name = options.([char(bound_type),'_bounds'])
                  invalid_name = true;
                  for type = C.var_types_notpar
                     if ismember(name,C.names.code.(type))
                        invalid_name = false;
                        C.cas.terminal_bounds.(bound_type).(name) = casadi.SX.sym([char(bound_type),'_terminal_bound_',char(name)]);
                        if bound_type == "lower"
                           expr = [expr; -C.cas.terminal_bounds.(bound_type).(name) + C.cas.(type).(name)]; %#ok<AGROW>
                        elseif bound_type == "upper"
                           expr = [expr;  C.cas.terminal_bounds.(bound_type).(name) - C.cas.(type).(name)]; %#ok<AGROW>
                        end
                     end
                  end
                  if invalid_name
                     warning(['WARNING: the name provided: "',char(name),'" is not a variable name. No bound is created for this name.'])
                  end
               end
            end
         end
         if C.flag_has_terminal_bounds
            C.flag_has_inequality = true;
            C.cas.constraints.terminal.bounds.expr = expr;
            Arg.out = C.cas.constraints.terminal.bounds.expr;
            C.cas.constraints.terminal.bounds.F = casadi.Function('terminal_bounds',Arg,infields,{'out'});
         end

         C.flag_terminal_constraints_defined = true;
         def_time = toc(def_time);
         disp(['done.  ',sec2str(def_time)])
      end
   
      %%%%%%%%%%%%%%%%%%%%%%%%%% HORIZON:
      function def_horizon(C,N,options)
         arguments
            C
            N  (1,1) double {mustBePositive,mustBeInteger}
            % options.Dt (1,1) double {mustBePositive}
            options.primaldual (1,1) string {mustBeMember(options.primaldual,["stacked","sorted"])} = "stacked";
         end
         
         % Prepare restoration
         fields = fieldnames(options);
         values = struct2cell(options);
         nameValuePairs = [fields'; values'];
         C.restore_cell{end}.def_horizon = [{N},reshape(nameValuePairs,1,[])];

         fprintf('Defining horizon...     ');
         def_time = tic;
         if ~C.flag_integrator_defied
            C.usererror('an integrator must be defined before one can define the horizon. Use: ".def_integrator" to define an integrator.')
         end

         C.horizon.N = N;
         if ~isfield(C.horizon,'Dt') && ~isfield(C.horizon,'T')
            C.set_Dt(1);
         else
            C.(['set_',char(C.horizon_decider)])(C.horizon.(C.horizon_decider));
         end

         %%%%%% Generate all decision variables on the horizon:
         C.cas.horizon.decision = structor("default_mix","TRYMPC_horizon");
         for var_type = C.var_types
            if var_type ~= "param"
               C.cas.horizon.decision.str.(var_type) = [];
               C.cas.horizon.ref.(var_type) = [];
               for k = 0: (N - (var_type == "input") )
                  C.cas.horizon.decision.str.(var_type) =  [C.cas.horizon.decision.str.(var_type) casadi.SX.sym(['k',num2str(k),'_',char(var_type)],[C.dim.(var_type),1])];
                  C.cas.horizon.ref.(var_type) =  [C.cas.horizon.ref.(var_type) casadi.SX.sym(['k',num2str(k),'_',char(var_type),'_ref'],[C.dim.(var_type),1])];
               end
            end
         end
         if C.integrator.has_aux
            for k = 0:(N-1)
               C.cas.horizon.decision.str.aux.(['stage_',num2str(k)]) = copy_aux(C.cas.integrator.aux.var.str);
            end
         end
         %%%%%% Internal function to help copy the auxiliary structor
         function S_copy = copy_aux(S)
            S_copy = struct;
            for name = string(fieldnames(S))'
               if isstruct(S.(name))
                  S_copy.(name) = copy_aux(S.(name));
               else
                  if numel(S.(name)) ~= 1
                     s1 = S.(name)(1);
                     var_name = s1.name;
                     ind = find(var_name=='_');
                     var_name = var_name(1:ind(end)-1);
                  else
                     var_name = char(name);
                  end
                  var_name = ['k',num2str(k),'_',var_name]; %#ok<AGROW>
                  S_copy.(name) = casadi.SX.sym(var_name,size(S.(name)));
               end
            end
         end




         %%%%%%%%%%% CONSTRAINTS:

         % Prepare a structor to contain all constraints:
         C.cas.horizon.constraints.equality = structor("default_mix","separated");
         if C.flag_has_inequality
            C.cas.horizon.constraints.inequality = structor("default_mix","separated");
         end

         
         %%%%%% Initial Value Embedding:
         C.cas.horizon.constraints.equality.str.stage_0.initial_value = C.cas.init.state - C.cas.horizon.decision.str.state(:,1);


         if C.flag_stage_constraints_defined
            %%%%%% Apply stage-wose constraints on horizon:
            Arg = struct;
            Arg.param = C.cas.var.param;
            for k = 0:C.horizon.N-1
               for var_type = C.var_types_notpar
                  Arg.(var_type) = C.cas.horizon.decision.str.(var_type)(:,k+1);
               end

               % Algebraic equations:
               if C.flag_has_algeb
                  C.cas.horizon.constraints.equality.str.(['stage_',num2str(k)]).algebraic = C.cas.algebraics.F.call(Arg).out;
               end

               % variable bounds
               if isfield(C.cas.constraints,'bounds')
                  C.cas.horizon.constraints.inequality.str.(['stage_',num2str(k)]).bounds = C.cas.constraints.bounds.F.call(Arg).out;
               end

               % general stage-wise constraints:
               if isfield(C.cas.constraints,'equality')
                  for eq = string(fieldnames(C.cas.constraints.equality))'
                     C.cas.horizon.constraints.equality.str.(['stage_',num2str(k)]).general.(eq) = C.cas.constraints.equality.(eq).F.call(Arg).out;
                  end
               end
               if isfield(C.cas.constraints,'inequality')
                  for in = string(fieldnames(C.cas.constraints.inequality))'
                     C.cas.horizon.constraints.inequality.str.(['stage_',num2str(k)]).general.(in) = C.cas.constraints.inequality.(in).F.call(Arg).out;
                  end
               end

            end
         end

         %%%%%% Build dynamic constraints on horizon:
         Arg = struct;
         Arg.Dt = C.cas.integrator.Dt;
         Arg.param = C.cas.var.param;
         counter = 0;
         for k = 1:C.horizon.N
            for var_type = C.var_types_notpar
               Arg.(var_type) = C.cas.horizon.decision.str.(var_type)(:,k);
            end
            if C.integrator.has_aux
               % Arg.aux = C.cas.horizon.decision.subvec(C.cas.horizon.decision.str.aux.(['stage_',num2str(k-1)]),"separated");
               Arg.aux = structor.subvec(C.cas.horizon.decision,C.cas.horizon.decision.str.aux.(['stage_',num2str(k-1)]),"separated");
            end

            % Auxiliary Integration constraints: (collocation constraints, IRK, ...)
            if C.integrator.has_aux
               C.cas.horizon.constraints.equality.str.(['stage_',num2str(k-1)]).dynamic.auxiliary = C.cas.integrator.aux.F.call(Arg).out;
               counter = counter + length(C.cas.horizon.constraints.equality.str.(['stage_',num2str(k-1)]).dynamic.auxiliary);
            end

            % Continuation Constraints: (closing the shooting gaps)
            C.cas.horizon.constraints.equality.str.(['stage_',num2str(k-1)]).dynamic.continuation = C.cas.horizon.decision.str.state(:,k+1) - C.cas.integrator.next.F.call(Arg).out;
         end
         C.display.problem.n_auxiliary_constraints = counter;
         




         %%%%%% Terminal constraints
         if C.flag_terminal_constraints_defined
            % C.cas.horizon.constraints.terminal.equality = structor("default_mix","TRYMPC_horizon");
            % C.cas.horizon.constraints.terminal.inequality = structor("default_mix","TRYMPC_horizon");
            Arg = struct;
            Arg.param = C.cas.var.param;
            k = C.horizon.N;
            for var_type = C.var_types_notpar(C.var_types_notpar ~= "input")
               Arg.(var_type) = C.cas.horizon.decision.str.(var_type)(:,k+1);
            end
            
            
            % variable bounds
            if isfield(C.cas.constraints.terminal,'bounds')
               C.cas.horizon.constraints.inequality.str.(['stage_',num2str(k)]).bounds = C.cas.constraints.terminal.bounds.F.call(Arg).out;
            end

            % general stage-wise constraints:
            if isfield(C.cas.constraints.terminal,'equality')
               for eq = string(fieldnames(C.cas.constraints.terminal.equality))'
                  C.cas.horizon.constraints.equality.str.(['stage_',num2str(k)]).general.(eq) = C.cas.constraints.terminal.equality.(eq).F.call(Arg).out;
               end
            end
            if isfield(C.cas.constraints.terminal,'inequality')
               for in = string(fieldnames(C.cas.constraints.terminal.inequality))'
                  C.cas.horizon.constraints.inequality.str.(['stage_',num2str(k)]).general.(in) = C.cas.constraints.terminal.inequality.(in).F.call(Arg).out;
               end
            end
         end











         C.flag_horizon_defied = true;
         def_time = toc(def_time);
         disp(['done.  ',sec2str(def_time)])

         C.def_problem(options.primaldual)
         C.find_stats




      end
   
   end



   %%%% Methods that help build an integrator:
   methods(Access = private,Hidden)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUNGE-KUTTA:
      function ERK_builder(C)
         
         % Buthcer Tableau:
         b = reshape(C.integrator.BT_b,[],1);
         A = C.integrator.BT_A;
         stages = length(b);

         % Prepare variables to simplify creation of integrator:
         in_vars = C.cas.var;
         next_state = C.cas.var.state;
         
         % Integrate system using RK method described by Butcher Tableau:
         for inc = 1:C.integrator.n_increments
            f = C.cas.dynamics.F.call(in_vars).out;
            for s = 2:stages   
                in_vars.state = next_state;
               for a = 1:s-1
                  in_vars.state = in_vars.state + A(s,a)*f(:,a)*(C.cas.integrator.Dt/C.integrator.n_increments);
               end
               f = [f C.cas.dynamics.F.call(in_vars).out]; %#ok<AGROW>
            end
            next_state = next_state + (f*b).*(C.cas.integrator.Dt/C.integrator.n_increments);
         end

         % Save expression for next state based:
         C.cas.integrator.next.expr = next_state;

         C.create_integrator_F;
      end


      function IRK_builder(C)
         
         % Buthcer Tableau:
         b = reshape(C.integrator.BT_b,[],1);
         A = C.integrator.BT_A;
         samp = length(b); % the number of sample points of the RK scheme (a.k.a. "stages", but that can be confused with stages as in each discrete time point on the prediction horizon)

         % First loop over and create relevant variables:
         C.cas.integrator.aux.var = structor;
         for inc = 1:C.integrator.n_increments
            for samp = 1:samp
               C.cas.integrator.aux.var.str.(['f_inc_',num2str(inc)]).(['samp',num2str(samp)]) = casadi.SX.sym(['RK_f_inc',num2str(inc),'_samp',num2str(samp)],[C.dim.state,1]);
            end
         end

         % Prepare variables to simplify creation of integrator:
         in_vars = C.cas.var;
         next_state = C.cas.var.state;
         C.cas.integrator.aux.expr = structor;
         
         % Integrate system using RK method described by Butcher Tableau:
         for inc = 1:C.integrator.n_increments
            for samp = 1:samp   
                in_vars.state = next_state;
               for a = 1:samp
                  in_vars.state = in_vars.state + A(samp,a)*C.cas.integrator.aux.var.str.(['f_inc_',num2str(inc)]).(['samp',num2str(a)])*(C.cas.integrator.Dt/C.integrator.n_increments);
               end
               C.cas.integrator.aux.expr.str.(['RK_f_expr_inc',num2str(inc),'_samp',num2str(samp)]) = C.cas.integrator.aux.var.str.(['f_inc_',num2str(inc)]).(['samp',num2str(a)]) - C.cas.dynamics.F.call(in_vars).out;
            end

            next_state = next_state + (reshape(structor.subvec(C.cas.integrator.aux.var,C.cas.integrator.aux.var.str.(['f_inc_',num2str(inc)])),C.dim.state,samp)*b) .* (C.cas.integrator.Dt/C.integrator.n_increments);
         end

         % Save expression for next state based:
         C.cas.integrator.next.expr = next_state;

         C.create_integrator_F;         
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ... END RUNGE-KUTTA.
%
%
%
%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COLLOCATION:
      function collocation_builder(C)
         C.Make_Collocation_Polynomial_Coefficients

         % Create symbolic variables
         C.cas.integrator.aux.var = structor;
         for inc = 1:C.integrator.n_increments
            for c = 1:C.integrator.collocation.d
               C.cas.integrator.aux.var.str.(['increment_',num2str(inc)]).(['aux_',num2str(c)]) = casadi.SX.sym(['collocation_var_inc',num2str(inc),'_samp',num2str(c)],[C.dim.state,1]);
            end
         end

         %%%%%%%%%%% Create symbolic collocation polynomial to use on each interval:
         C.cas.integrator.collocation.tau = casadi.SX.sym('internal_tau');
         tau_powers = transpose(C.cas.integrator.collocation.tau.^(C.integrator.collocation.d:-1:0));

         % Prepare Arg:
         Arg = struct;
         Arg.state = C.cas.var.state;
         % Arg.aux = C.cas.integrator.aux.var.subvec(C.cas.integrator.aux.var.str.increment_1);
         Arg.aux = structor.subvec(C.cas.integrator.aux.var,C.cas.integrator.aux.var.str.increment_1);

         % Use collocation variables of first step as collocation variables
         % col_vars = C.cas.integrator.aux.var.subvec(C.cas.integrator.aux.var.str.increment_1);
         col_vars = structor.subvec(C.cas.integrator.aux.var,C.cas.integrator.aux.var.str.increment_1);
         col_vars = reshape(col_vars,C.dim.state,C.integrator.collocation.d);
         col_vars = [C.cas.var.state col_vars];
         
         % Create collocation polynomial p
         C.cas.integrator.collocation.p.expr = col_vars * C.integrator.collocation.Li * tau_powers;
         
         % Define casadi function
         Arg.tau = C.cas.integrator.collocation.tau;
         infields = fieldnames(Arg);
         Arg.out = C.cas.integrator.collocation.p.expr;
         C.cas.integrator.collocation.p.F = casadi.Function('F_collocation_p',Arg,infields,{'out'});


         %%%%%%%%%%%% Also create the derivative of p, i.e: dp
         tau_powers = transpose(C.cas.integrator.collocation.tau.^(C.integrator.collocation.d-1:-1:0));

         % Create collocation polynomial derivative; dp
         C.cas.integrator.collocation.dp.expr = col_vars * C.integrator.collocation.dLi * tau_powers;

         % Create casadi function
         Arg.out = C.cas.integrator.collocation.dp.expr;
         C.cas.integrator.collocation.dp.F = casadi.Function('F_collocation_dp',Arg,infields,{'out'});


         %%%%%%%%%% Create collocation expressions for each step
         C.cas.integrator.aux.expr = structor;
         Arg = struct;
         Arg_dyn = C.cas.var;
         for inc = 1:C.integrator.n_increments
            % Arg.aux = C.cas.integrator.aux.var.subvec(C.cas.integrator.aux.var.str.(['increment_',num2str(inc)]));
            Arg.aux = structor.subvec(C.cas.integrator.aux.var,C.cas.integrator.aux.var.str.(['increment_',num2str(inc)]));
            Arg.tau = 0;
            Arg.state = Arg_dyn.state;
            % C.cas.integrator.aux.expr.str.(['step_',num2str(n)]).stage_1 = Arg_dyn.state - C.cas.integrator.collocation.p.F.call(Arg).out;
            for c = 1:C.integrator.collocation.d
               Arg.tau = C.integrator.collocation.tau(c+1);
               Arg_dyn.state = C.cas.integrator.aux.var.str.(['increment_',num2str(inc)]).(['aux_',num2str(c)]);
               C.cas.integrator.aux.expr.str.(['increment_',num2str(inc)]).(['aux_',num2str(c)]) = C.cas.dynamics.F.call(Arg_dyn).out*(C.cas.integrator.Dt/C.integrator.n_increments) - C.cas.integrator.collocation.dp.F.call(Arg).out;
            end
            Arg.tau = 1;
            Arg_dyn.state = C.cas.integrator.collocation.p.F.call(Arg).out;
         end
         % Finally, the final/('next') state variable is:
         C.cas.integrator.next.expr = Arg_dyn.state;

         %%%%%%%%%%% Create casadi functions of the integration step:
         C.create_integrator_F

      end

      function Make_Collocation_Polynomial_Coefficients(C)
         d = C.integrator.collocation.d; % order of integration and order of polynomial

         C.integrator.collocation.tau = [0 casadi.collocation_points(d, char(C.integrator.collocation.polynomial_type))];
         C.integrator.collocation.Li = [];
         C.integrator.collocation.dLi = [];
         for j=1:d+1
            coeff = 1;
            for r=1:d+1
               if r ~= j
                  coeff = conv(coeff, [1, -C.integrator.collocation.tau(r)]);
                  coeff = coeff / (C.integrator.collocation.tau(j)-C.integrator.collocation.tau(r));
               end
            end
            C.integrator.collocation.Li(j,:) = coeff;
            C.integrator.collocation.dLi(j,:) = polyder(coeff);
         end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ... END COLLOCATION.


      % Generate the casadi functions for auxiliary expression and next
      % state:
      function create_integrator_F(C)
         if ~isfield(C.cas,'integrator') || ~isfield(C.cas.integrator,'next') || ~isfield(C.cas.integrator.next,'expr')
            error('Developer Error: do not call this funciton if the expression for the next state is not defined.')
         end

         if C.integrator.has_aux
            Arg = C.cas.var;
            Arg.Dt = C.cas.integrator.Dt;
            Arg.aux = C.cas.integrator.aux.var.vec;
            infields = fieldnames(Arg);
            Arg.out = C.cas.integrator.aux.expr.vec;
            C.cas.integrator.aux.F = casadi.Function('F_aux',Arg,infields,{'out'});
         end

         Arg = C.cas.var;
         Arg.Dt = C.cas.integrator.Dt;
         if C.integrator.has_aux
            Arg.aux = C.cas.integrator.aux.var.vec;
         end
         infields = fieldnames(Arg);
         Arg.out = C.cas.integrator.next.expr;
         C.cas.integrator.next.F = casadi.Function('F_next_state',Arg,infields,{'out'});
      end
   end




   %%%% Internal methods:
   methods(Access = private)
      function def_problem(C,primaldual_type)
         fprintf('Defining problem...     ');
         def_time = tic;
         % This function defines a problem formulation based on the horizon
         % that was created by "def_horizon"

         C.cas.problem.decision = C.cas.horizon.decision; % this is a handle class, so the values are not copied. (think of it as a pointer)

         %%%%%%%%%%%%%%%%%%%% Generate total cost on horizon:
         C.cas.problem.objective.expr = casadi.SX.zeros;
         Arg_objective_stage = C.cas.objective.weights;
         Arg_objective_stage.param = C.cas.var.param;
         for k = 1:C.horizon.N
            for var_type = C.var_types_notpar
               Arg_objective_stage.(var_type) = C.cas.horizon.decision.str.(var_type)(:,k);
               Arg_objective_stage.([char(var_type),'_ref']) = C.cas.horizon.ref.(var_type)(:,k);
               Arg_objective_stage.([char(var_type),'_next']) = C.cas.horizon.decision.str.(var_type)(:,k + ~(k == C.horizon.N && var_type == "input") );
            end
            C.cas.problem.objective.expr = C.cas.problem.objective.expr + C.cas.objective.cost.stage_total.F.call(Arg_objective_stage).out;
         end
         %%%%% Add terminal cost:
         k = C.horizon.N+1;
         Arg_objective_terminal= C.cas.objective.weights;
         for var_type = C.var_types_notpar
            if var_type ~= "input"
               Arg_objective_terminal.(var_type) = C.cas.horizon.decision.str.(var_type)(:,k);
               Arg_objective_terminal.([char(var_type),'_ref']) = C.cas.horizon.ref.(var_type)(:,k);
            end
         end
         C.cas.problem.objective.expr = C.cas.problem.objective.expr + C.cas.objective.cost.terminal_total.F.call(Arg_objective_terminal).out;

         
         % Create casadi functions for objective:
         % Arg_objective = C.cas.objective.weights;
         % Arg_objective.param = C.cas.var.param;
         % for var_type = C.var_types_notpar
         %    Arg_objective.([char(var_type),'_ref_trajectory'])  = C.cas.horizon.ref.(var_type);
         % end
         % Arg_objective.decision = C.cas.problem.decision.vec;
         % in_fields = fieldnames(Arg_objective);
         [Arg_objective,in_fields] = C.create_Arg_symbolic_objective;

         C.cas.problem.objective.in_fields = in_fields;
         Arg_objective.out = C.cas.problem.objective.expr;
         C.cas.problem.objective.F = casadi.Function('F_objective',Arg_objective,in_fields,{'out'});



         %%%%%%%%%%%%%%%%%%%%% Constraints:
         % Arg_constraints = struct;
         % Arg_constraints.Dt = C.cas.integrator.Dt;
         % Arg_constraints.initial_state = C.cas.init.state;
         % if C.flag_has_param
         %    Arg_constraints.param = C.cas.var.param;
         % end
         % if C.flag_has_bounds
         %    for bound_type = string(fieldnames(C.cas.bounds))'
         %       for name = string(fieldnames(C.cas.bounds.(bound_type)))'
         %          Arg_constraints.([char(bound_type),'_',char(name)]) = C.cas.bounds.(bound_type).(name);
         %       end
         %    end
         % end
         % if C.flag_has_terminal_bounds
         %    for bound_type = string(fieldnames(C.cas.terminal_bounds))'
         %       for name = string(fieldnames(C.cas.terminal_bounds.(bound_type)))'
         %          Arg_constraints.([char(bound_type),'_terminal_',char(name)]) = C.cas.terminal_bounds.(bound_type).(name);
         %       end
         %    end
         % end
         % Arg_constraints.decision = C.cas.problem.decision.vec;
         % in_fields = fieldnames(Arg_constraints);
         [Arg_constraints,in_fields] = C.create_Arg_symbolic_constraints;

         %%% Equality:
         C.cas.problem.equality.in_fields = in_fields;
         C.cas.problem.equality.expr = C.cas.horizon.constraints.equality.vec;

         Arg_constraints.out = C.cas.problem.equality.expr;
         C.cas.problem.equality.F = casadi.Function('F_equality',Arg_constraints,in_fields,{'out'});

         %%% Inequality:
         if C.flag_has_inequality
            C.cas.problem.inequality.in_fields = in_fields;
            C.cas.problem.inequality.expr = C.cas.horizon.constraints.inequality.vec;

            Arg_constraints.out = C.cas.problem.inequality.expr;
            C.cas.problem.inequality.F = casadi.Function('F_inequality',Arg_constraints,in_fields,{'out'});
         end




         %%%%%%%%%%%%%%%%%%%%% Derivatives:

         %%%%% Objective:

         % Jacobian
         C.cas.problem.objective.J.expr = jacobian(C.cas.problem.objective.expr,C.cas.problem.decision.vec);
         Arg_objective.out = C.cas.problem.objective.J.expr;
         C.cas.problem.objective.J.F = casadi.Function('F_J_objective',Arg_objective,C.cas.problem.objective.in_fields,{'out'});
         % Hessian
         C.cas.problem.objective.H.expr = hessian(C.cas.problem.objective.expr,C.cas.problem.decision.vec);
         Arg_objective.out = C.cas.problem.objective.H.expr;
         C.cas.problem.objective.H.F = casadi.Function('F_H_objective',Arg_objective,C.cas.problem.objective.in_fields,{'out'});

         %%% Constraints:

         % Equality - jacobian:
         C.cas.problem.equality.J.expr = jacobian(C.cas.problem.equality.expr,C.cas.problem.decision.vec);
         Arg_constraints.out = C.cas.problem.equality.J.expr;
         C.cas.problem.equality.J.F = casadi.Function('F_J_equality',Arg_constraints,C.cas.problem.equality.in_fields,{'out'});
         
         % Inequality - jacobian:
         if C.flag_has_inequality
            C.cas.problem.inequality.J.expr = jacobian(C.cas.problem.inequality.expr,C.cas.problem.decision.vec);
            Arg_constraints.out = C.cas.problem.inequality.J.expr;
            C.cas.problem.inequality.J.F = casadi.Function('F_J_inequality',Arg_constraints,C.cas.problem.inequality.in_fields,{'out'});
         end
         



         %%%%%%%%%%%%%%%%%%%%% Lagrangian:

         C.generate_multipliers
         lambda = C.cas.problem.lambda.vec;
         if C.flag_has_inequality
            mu = C.cas.problem.mu.vec;
         else
            mu = [];
         end

         % Prepare arguments
         Arg_lagrangian = mergestructs(Arg_objective,Arg_constraints);
         Arg_lagrangian.lambda = lambda;
         if C.flag_has_inequality
            Arg_lagrangian.mu = mu;
         end
         Arg_lagrangian = rmfield(Arg_lagrangian,'out');
         C.cas.problem.Lagrangian.in_fields = fieldnames(Arg_lagrangian);


         % create LAgrangian expression (L = obj + lamnda'*g - mu'*h) where: h >= 0, and mu >= 0.
         if C.flag_has_inequality
            C.cas.problem.Lagrangian.expr = C.cas.problem.objective.expr + transpose(lambda)*C.cas.problem.equality.expr - transpose(mu)*C.cas.problem.inequality.expr;
         else
            C.cas.problem.Lagrangian.expr = C.cas.problem.objective.expr + transpose(lambda)*C.cas.problem.equality.expr;
         end
         Arg_lagrangian.out = C.cas.problem.Lagrangian.expr;
         C.cas.problem.Lagrangian.F    = casadi.Function('F_Lagrangian',Arg_lagrangian,C.cas.problem.Lagrangian.in_fields,{'out'});

         %%%% Jacobian - This is the jacobian of the lagrangian, which should be zero to satisfy the 'stationarity condition'
         C.cas.problem.Lagrangian.J.expr = jacobian(C.cas.problem.Lagrangian.expr,C.cas.problem.decision.vec);
         Arg_lagrangian.out = C.cas.problem.Lagrangian.J.expr;
         C.cas.problem.Lagrangian.J.F    = casadi.Function('F_J_Lagrangian',Arg_lagrangian,C.cas.problem.Lagrangian.in_fields,{'out'});

         %%%% Hessian:
         C.cas.problem.Lagrangian.H.expr = hessian(C.cas.problem.Lagrangian.expr,C.cas.problem.decision.vec);
         Arg_lagrangian.out = C.cas.problem.Lagrangian.H.expr;
         C.cas.problem.Lagrangian.H.F    = casadi.Function('F_H_Lagrangian',Arg_lagrangian,C.cas.problem.Lagrangian.in_fields,{'out'});

         %%%%%%%%%%%%%%%%%%%% KKT system:

         % create primal-dual vector
         switch primaldual_type
            case "stacked"
               C.cas.problem.primaldual = structor("default_mix","separated");
               C.cas.problem.primaldual.str.primal = C.cas.problem.decision.vec;
               C.cas.problem.primaldual.str.lambda.all = C.cas.problem.lambda.vec;
               if C.flag_has_inequality
                  C.cas.problem.primaldual.str.mu.all = C.cas.problem.mu.vec;
               end
               primaldual = C.cas.problem.primaldual.vec;
            case "sorted"
               C.cas.problem.primaldual = structor("default_mix","TRYMPC_horizon");
               C.cas.problem.primaldual.str = C.cas.problem.decision.str;
               C.cas.problem.primaldual.str.lambda = C.cas.problem.lambda.str;
               if C.flag_has_inequality
                  C.cas.problem.primaldual.str.mu = C.cas.problem.mu.str;
               end
               primaldual = C.cas.problem.primaldual.vec;
         end


         %%%% Jacobian - This is the KKT vector (the KKT conditions are: KKT-vector = 0)
         C.cas.problem.KKT.vector.expr = jacobian(C.cas.problem.Lagrangian.expr,primaldual);
         Arg_lagrangian.out = C.cas.problem.KKT.vector.expr;
         C.cas.problem.KKT.vector.F    = casadi.Function('F_KKT_vector',Arg_lagrangian,C.cas.problem.Lagrangian.in_fields,{'out'});

         %%%% Hessian - This is the KKT matrix (used to perform newton steps in the KKT vector, in order to find its roots (zeros))
         C.cas.problem.KKT.matrix.expr = hessian(C.cas.problem.Lagrangian.expr,primaldual);
         Arg_lagrangian.out = C.cas.problem.KKT.matrix.expr;
         C.cas.problem.KKT.matrix.F    = casadi.Function('F_H_Lagrangian',Arg_lagrangian,C.cas.problem.Lagrangian.in_fields,{'out'});

         
         
         C.flag_problem_defined = true;
         def_time = toc(def_time);
         disp(['done.  ',sec2str(def_time)])
      end
   
      function find_stats(C)
         %%%%%%%%%% Prepare display values:
         
         %%% Decision Variables:
         C.display.problem.n_decision = C.cas.problem.decision.len;
         C.display.problem.n_state    = numel(C.cas.problem.decision.str.state);
         if C.flag_has_algeb
            C.display.problem.n_algeb    = numel(C.cas.problem.decision.str.algeb);
         end
         if C.flag_has_input
            C.display.problem.n_input    = numel(C.cas.problem.decision.str.input);
         end
         if C.integrator.has_aux
            % C.display.problem.n_aux      = numel(C.cas.problem.decision.subvec(C.cas.problem.decision.str.aux));
            C.display.problem.n_aux      = numel(structor.subvec(C.cas.problem.decision,C.cas.problem.decision.str.aux));
         end

         %%% Constraints
         
         % ineq vs. eq:
         C.display.problem.n_equality = length(C.cas.problem.equality.expr);
         if C.flag_has_inequality
            C.display.problem.n_inequality = length(C.cas.problem.inequality.expr);
         else
            C.display.problem.n_inequality = 0;
         end

         % dynamic constraints
         % v = C.cas.horizon.constraints.equality.subvec(C.cas.horizon.constraints.equality.str.dynamic);
         % v = structor.subvec(C.cas.horizon.constraints.equality,C.cas.horizon.constraints.equality.str.steps);
         C.display.problem.n_continuation_constraints = C.horizon.N*C.dim.state;
         C.display.problem.n_dynamic_constraints = C.display.problem.n_continuation_constraints + C.display.problem.n_auxiliary_constraints;

         % algebraic constraints
         if C.flag_has_algeb
            % v = C.cas.horizon.constraints.equality.subvec(C.cas.horizon.constraints.equality.str.algebraic);
            v = structor.subvec(C.cas.horizon.constraints.equality,C.cas.horizon.constraints.equality.str.algebraic);
            C.display.problem.n_algebraic_constraints = length(v);
         end

         % other constraints
         if C.flag_stage_constraints_defined && isfield(C.cas.constraints,'equality')
            % v = C.cas.horizon.constraints.equality.subvec(C.cas.horizon.constraints.equality.str.general);
            v = structor.subvec(C.cas.horizon.constraints.equality,C.cas.horizon.constraints.equality.str.general);
            C.display.problem.n_other_constraints = length(v);
         end

         % Stage-wise constraints:
         C.display.problem.stage_constraints.n_equality = struct;
         C.display.problem.stage_constraints.n_inequality = struct;
         if isfield(C.cas,'constraints')
            if isfield(C.cas.constraints,'equality')
               eq_names = string(fieldnames(C.cas.constraints.equality))';
               if ~isempty(eq_names)
                  for name = eq_names
                     C.display.problem.stage_constraints.n_equality.(name) = size(C.cas.constraints.equality.(name).expr,1)*C.horizon.N;
                  end
               end
            end
            if isfield(C.cas.constraints,'inequality')
               in_names = string(fieldnames(C.cas.constraints.inequality))';
               if ~isempty(in_names)
                  for name = in_names
                     C.display.problem.stage_constraints.n_inequality.(name) = size(C.cas.constraints.inequality.(name).expr,1)*C.horizon.N;
                  end
               end
            end
         end



         %%% Jacobians
         for name = ["objective","equality","inequality"]
            if ((name ~= "inequality") + (C.flag_has_inequality))
               C.display.problem.(['J_',char(name),'_nnz'])             = nnz(C.cas.problem.(name).J.expr);
               C.display.problem.(['J_',char(name),'_numel'])           = numel(C.cas.problem.(name).J.expr);
               C.display.problem.(['J_',char(name),'_sparsity'])        = (C.display.problem.(['J_',char(name),'_numel']) - C.display.problem.(['J_',char(name),'_nnz']))/C.display.problem.(['J_',char(name),'_numel']);
               C.display.problem.(['J_',char(name),'_linear_density'])  = (C.display.problem.(['J_',char(name),'_nnz']))/sqrt(C.display.problem.(['J_',char(name),'_numel']));
            end
         end



         %%% Objective types:
         C.display.problem.cost_types = strjoin(C.cost_types, ', ');
      end
   
      function generate_multipliers(C)
         
         % equality multipliers:
         C.cas.problem.lambda = structor("default_mix",C.cas.horizon.constraints.equality.default_mix);
         C.cas.problem.lambda.str = loop_structure(C.cas.horizon.constraints.equality.str,'Lag_mult_lambda');

         % inequality multipliers:
         if C.flag_has_inequality
            C.cas.problem.mu = structor("default_mix",C.cas.horizon.constraints.inequality.default_mix);
            C.cas.problem.mu.str = loop_structure(C.cas.horizon.constraints.inequality.str,'Lag_mult_mu');
         end

         function S_out= loop_structure(S,Name)
            S_out = struct;
            for name = string(fieldnames(S))'
               Name_new = [Name,'_',char(name)];
               if isstruct(S.(name))
                  S_out.(name) = loop_structure(S.(name),Name_new);
               else
                  S_out.(name) = casadi.SX.sym(Name_new,size(S.(name)));
               end
            end
         end
      end
   end











   %%%%%%%%%%% Simulation
   methods(Access = public)
      function simulate(C,duration,init_state,options,mpc_options)
         arguments
            C
            
            duration (1,1) double {mustBePositive}
            init_state (:,1) {mustBeA(init_state, {'double', 'struct'})} = C.initial_state.vec

            % Simulator settings:
            options.simulator (1,1) string {mustBeMember(options.simulator,["ode45","ode23","ode113","ode78","ode89","ode15s","ode23s","ode23t","ode23tb","ode15i"])} = "ode45";
            options.ode_options (1,1) struct = odeset('RelTol',10^(-7),'AbsTol',10^(-7));
            options.start_time (1,1) double {mustBeReal} = 0;

            % controller settings:
            options.sampling_time (1,1) double {mustBePositive} = inf;
            options.controller_type (1,1) string {mustBeMember(options.controller_type,["none","constant","LMPC","NMPC_ipopt","NMPC_sqp","asNMPC","RTI","aRTI","aRTIshp","custom"])} = "none";
            options.QP_solver (1,1) string {mustBeMember(options.QP_solver,["quadprog (interior point)", "quadprog (active set)", "quadprog (trust region)", "custom"])} = "quadprog (interior point)";
            options.control_delay_scaler (1,1) double {mustBeNonnegative} = 0;

            %%% Specific Controller Settings: 
            % (constant controller)
            options.controller_constant (:,1) double {mustBeReal}
            % (custom controller)
            options.controller_custom (1,1) function_handle
            
            % nlpsol (ipopt) options
            options.ipopt_nlpsol_options (1,1) struct = struct('print_time',0);


            % Disturbances
            options.input_disturbance (1,1) function_handle = @(~,~,~) inf;  % input disturbance (additive) (on the form: @(C,t,x) (...) )

            % Miscellaneous
            options.speed (1,1) logical = false; % Setting this to true will make it disregard storing data along the way, and avoid doing unnecessary things in the controller and simulator, such that the simulation will finish faster. The down-side is that little information about the simulation is stored, such that debugging and analyzing becomes harder.
            options.halting_condition (1,1) function_handle = @(s,i,p) 0; % if this outputs a prositive value, than we halt the simulation. (arguments s,i,p are vectors, not structs)

            %%% MPC options:
            mpc_options.initial_guess_primal (:,1) double {mustBeReal}
            mpc_options.initial_guess_dual_eq (:,1) double {mustBeReal}
            mpc_options.initial_guess_dual_in (:,1) double {mustBeReal}
            % mpc_options.state_ref (:,:) double
            % mpc_options.algeb_ref (:,:) double 
            % mpc_options.input_ref (:,:) double
         end

         if ~C.flag_dynamics_defined
            C.usererror('must define dynamics before simulating.')
         end

         % try to be more efficient?:
         C.speed = options.speed;



         % Store settings:
         C.simulation.simulator     = options.simulator;
         C.simulation.initial_state = init_state;
         C.simulation.duration      = duration;
         C.simulation.start_time    = options.start_time;
         C.simulation.controller    = options.controller_type;
         C.simulation.ode_options   = options.ode_options;
         C.simulation.sampling_time = options.sampling_time;
         C.simulation.control_delay_scaler = options.control_delay_scaler;
         if isinf(options.input_disturbance(C,options.start_time,init_state))
            C.simulation.input_disturbance = "none";
         else
            C.simulation.input_disturbance = options.input_disturbance;
         end

         % Assume extarnal controller, and set to true if internal is used:
         C.internal_mpc.internal_controller = false;

         % Configure controller:
         switch options.controller_type
            case "none"
               C.simulation.controller_handle = @(C,t,x) zeros(C.dim.input,1);
            case "constant"
               if C.dim.input == length(options.controller_constant)
                  C.simulation.controller_handle = @(C,t,x) options.controller_constant;
               else
                  C.usererror('constant controller must be same dimension as the input vector.')
               end
            case "custom"
               C.simulation.controller_handle = options.controller_custom; % this handle takes the class instance as an argument, and should use the C.set_contorl_delay() to set the contorl delay
            otherwise
               
               % Prepare mpc_options:
               mpc_options.initial_state = init_state;
               if ~isfield(mpc_options,'initial_guess_primal')
                  mpc_options.initial_guess_primal = zeros(C.display.problem.n_decision,1);
               end
               if ~isfield(mpc_options,'initial_guess_dual_eq')
                  mpc_options.initial_guess_dual_eq = zeros(C.display.problem.n_equality,1);
               end
               if ~isfield(mpc_options,'initial_guess_dual_in')
                  mpc_options.initial_guess_dual_in = zeros(C.display.problem.n_inequality,1);
               end

               % make sure to start fresh, uncorreupted by previous
               % simulations/optimizations:
               C.mpc_start_fresh(mpc_options)

               % Add initial fields:
               C.internal_mpc.internal_controller = true;
               C.internal_mpc.solve_time.total = 0;
               C.internal_mpc.archive_time = 0;


               switch options.controller_type
                  case "NMPC_ipopt"
                     if ~isfield(options.ipopt_nlpsol_options,'ipopt')
                        options.ipopt_nlpsol_options.ipopt.print_level = 0;
                        options.ipopt_nlpsol_options.ipopt.max_iter = 100;

                        % Default values:
                        options.ipopt_nlpsol_options.ipopt.tol = 1e-4; % default: 1e-8;
                        options.ipopt_nlpsol_options.ipopt.acceptable_tol = 1e-2; % default: 1e-6;
                        options.ipopt_nlpsol_options.ipopt.compl_inf_tol = 1e-4;
                        options.ipopt_nlpsol_options.ipopt.constr_viol_tol = 1e-4;
                        options.ipopt_nlpsol_options.ipopt.dual_inf_tol = 1e-4;

                        %%% Some other ipopt-options:
                        % options.ipopt_nlpsol_options.ipopt.hessian_approximation      = 'limited-memory';
                        % options.ipopt_nlpsol_options.ipopt.limited_memory_update_type = 'bfgs';
                        % options.ipopt_nlpsol_options.ipopt.linear_solver              = 'mumps';
                        % options.ipopt_nlpsol_options.ipopt.linear_system_scaling      = 'none';
                     end
                     C.internal_mpc.nlpsol_options = options.ipopt_nlpsol_options;
                     C.pre_initialize_ipopt
                     C.initialize_ipopt
                  case "LMPC"
                  case "NMPC_sqp"
                  case "asNMPC" 
                  case "RTI" 
                  case "aRTI" 
                  case "aRTIshp"
               end
               C.simulation.controller_handle = @(C,t,x) C.(['controller_',char(options.controller_type)])(t,x);
         end



         %%%%%%% INITIALIZE SIMULATION:
         C.internal_sim.simulator_handle = str2func(C.simulation.simulator);
         C.internal_sim.time = C.simulation.start_time;
         C.internal_sim.state = C.simulation.initial_state;
         Arg = struct;
         arg_param = @() setfield(Arg,'param',C.parameters.vec);
         arg_state = @(state) setfield(arg_param(),'state',state);
         if C.flag_has_algeb
            arg_algeb = @(state,algeb) setfield(arg_state(state),'algeb',algeb);
         else
            arg_algeb = @(state,~) arg_state(state);
         end
         if C.flag_has_input
            arg_input = @(state,algeb,input) setfield(arg_algeb(state,algeb),'input',input);
         else
            arg_input = @(state,algeb,~) arg_algeb(state,algeb);
         end
         
         C.internal_sim.arg = @(x,u) arg_input(x(1:C.dim.state),x((C.dim.state+1):end),u);


         if C.flag_has_algeb
            if ~ismember(C.simulation.simulator,C.allowable_DAE_solver)
               C.usererror(['You have a differential algebraic equation (DAE), thus you have to choose a simulator that can handle DAEs. Try either of: ',char(strjoin(C.allowable_DAE_solver,', '))])
            end

            Arg = struct;
            Arg.param = C.parameters.vec;
            Arg.input = zeros(C.dim.input,1);
            arg = @(y) cell2struct([struct2cell(Arg); {y(1:C.dim.state); y((C.dim.state+1):end)}], [fieldnames(Arg); {'state'; 'algeb'}]);
            odefun = @(t,y,yp) full([yp(1:C.dim.state) - C.cas.dynamics.F.call(arg(y)).out;
                                     C.cas.algebraics.F.call(arg(y)).out]);

            initial_algeb = zeros(C.dim.algeb,1);
            initial_algeb_dot = zeros(C.dim.algeb,1);
            initial_y = [C.simulation.initial_state; initial_algeb];
            y_fixed =   [ones(C.dim.state,1); zeros(C.dim.algeb,1)];
            initial_yp = full([C.cas.dynamics.F.call(arg(initial_y)).out;
                          initial_algeb_dot]);
            C.simulation.initial_algeb = decic(odefun,C.simulation.start_time,initial_y,y_fixed,initial_yp,[]);
            C.internal_sim.Arg.algeb = C.simulation.initial_algeb;

            % add mass matrix to tell the solver that this is a DAE:
            C.simulation.ode_options.Mass = diag(y_fixed);
         end


         

         %%%%%%% Before simulating, prepare the sim_results struct
         C.archive.simulations{end+1} = C.simulation;
         C.archive.simulations{end}.ID = C.generate_id;
         C.archive.simulations{end}.stage = struct([]);
         C.archive.simulations{end}.sim.time = C.simulation.start_time;
         C.archive.simulations{end}.sim.state = C.simulation.initial_state;
         C.archive.simulations{end}.sim.input_controller = [];
         C.archive.simulations{end}.optimizations = {};

         %%%% Prepare stuff:

         % Find integration interval when control delay is present:
         until_sample = @(control_delay) min(C.simulation.sampling_time - control_delay,  C.simulation.duration - (C.internal_sim.time - C.simulation.start_time));

         % if the control signal is to be computed continuously (at each point within the ode-integrator)
         if isa(C.simulation.input_disturbance,'function_handle')
            if isinf(C.simulation.sampling_time)
               C.internal_sim.ode_inputs = @(t,x) C.simulation.controller_handle(C,t,x) + C.simulation.input_disturbance(C,t,x);
            else
               C.internal_sim.ode_inputs = @(t,x) C.internal_sim.u_current + C.simulation.input_disturbance(C,t,x);
            end
         else
            if isinf(C.simulation.sampling_time)
               C.internal_sim.ode_inputs = @(t,x) C.simulation.controller_handle(C,t,x);
            else
               C.internal_sim.ode_inputs = @(t,x) C.internal_sim.u_current;
            end
         end

         if C.simulation.control_delay_scaler
            C.internal_sim.delay_period = false; % "delay_period":  "true" - time from a state measurement to control action is applied, "false" - time after control is updated to the next measurement
         end
         
         
         %%% Prepare displays
         fprintf(['Simulating (',char(datetime('now','Format','HH:mm')),') ... ']);
         C.internal_sim.linelength = 0;
         total_dot_length = 3; dot_length = 0;
         sim_time = tic;

         %%%%%%% SIMULATE:
         while C.internal_sim.time < C.simulation.duration + C.simulation.start_time
            
            if options.halting_condition(C.internal_sim.state,C.internal_sim.u_current,C.parameters.vec) > 0
               break;
            end
            display_progress

            if C.simulation.control_delay_scaler
               % If new state is measured, find suitable control action based on current state
               if C.flag_has_input
                  C.internal_sim.delay_period = ~C.internal_sim.delay_period;

                  if C.internal_sim.delay_period % go here in the first pass
                     %{
                      A NEW MEASUREMENT HAS BEEN MADE, AND WE WILL COMPUTE
                      THE CORRESPONDING CONTROL ACTION.
                      BEFORE APPLYING, WE WILL SIMULATE FOR THE AMOUNT OF
                      TIME SPENT COMPUTING THE CONTORL.
                     %}

                     % Find new input:
                     C.internal_sim.u_controller = C.simulation.controller_handle(C,C.internal_sim.time,C.internal_sim.state);

                     % Handle control delay:
                     control_delay = C.internal_sim.control_delay;
                     C.internal_sim.simulation_increment = control_delay;
                  else
                     %{
                      THE CONTROLLER IS NOW DONE COMPUTING, AND WE CAN
                      SIMULATE WITH THE UPDATED CONTROLLER UNTIL THE NEXT
                      SAMPLE.
                     %}

                     % Apply the previously computed control action: ('which just now finished computing, and is ready to be applied'...)
                     C.internal_sim.u_current = C.internal_sim.u_controller; % the "current" u is the current signal that the controller! wants.

                     % Simulate forward until the next sampling (or until the end of simulation if that come first...)
                     C.internal_sim.simulation_increment = until_sample(control_delay);
                  end

                  % At the start of simulation, we ignore control delay:
                  if C.internal_sim.time == C.simulation.start_time
                     C.internal_sim.u_current = C.internal_sim.u_controller; % Start off with the computed control action
                  end
               end

            else
               if ~isinf(C.simulation.sampling_time)
                  C.internal_sim.u_current = C.simulation.controller_handle(C,C.internal_sim.time,C.internal_sim.state);
               end
               C.internal_sim.simulation_increment = until_sample(0);
            end


            if  C.internal_sim.simulation_increment
               % Simulate the current increment
               [C.internal_sim.t_sim,C.internal_sim.x_sim] = C.internal_sim.simulator_handle( ...
                  @(t,x) full(C.cas.dynamics.F.call(C.internal_sim.arg(x,C.internal_sim.ode_inputs(t,x))).out), ...  % system dynamics
                  C.internal_sim.time+[0,C.internal_sim.simulation_increment],...  % time interval of integration
                  C.internal_sim.state,...  % initial state (current state of the greater simulation)
                  C.simulation.ode_options); % ode options


               % log data:
               C.log_simulation
            end


            C.internal_sim.time = C.internal_sim.t_sim(end);
            C.internal_sim.state = C.internal_sim.x_sim(end,:)';
         end

         if options.halting_condition(C.internal_sim.state,C.internal_sim.u_current,C.parameters.vec) > 0
            disp(' ')
            disp('... Halting Condition met. ')
         else
            display_progress
         end

         function display_progress
            dot_length = mod(dot_length,total_dot_length) + 1;
            C.internal_sim.linelength = fprintf([repmat('\b',1,C.internal_sim.linelength),repmat('.',1,dot_length),repmat(' ',1,total_dot_length-dot_length),' %3.2f%%   '], 100*(C.internal_sim.time - C.simulation.start_time)/C.simulation.duration ) - C.internal_sim.linelength;
            if C.internal_mpc.internal_controller
               C.internal_sim.linelength = C.internal_sim.linelength + fprintf('| solve time: (%2.5fs) | archive time: (%2.5fs)    ', C.internal_mpc.solve_time.total , C.internal_mpc.archive_time);
            end
         end
         %%%%%%%%%%%%%%%%% simulation done..
         C.archive.simulations{end}.simulation_time = toc(sim_time);
         disp(' ')
         disp(['Simulation Finished. ',sec2str(C.archive.simulations{end}.simulation_time)])


         % reocompute and log all inputs if relevant:
         if C.flag_has_input
            % fprintf(' - archiving input signals... ')
            % archive_time = tic;
            C.archive.simulations{end}.sim.input_disturbance = [];
            % C.archive.simulations{end}.sim.input_controller = [C.archive.simulations{end}.stage(1).input_controller C.archive.simulations{end}.sim.input_controller];
            x_sim = C.archive.simulations{end}.sim.state;

            % append algebraic variables if there are any
            if C.flag_has_algeb
               x_sim = [x_sim; C.archive.simulations{end}.sim.algeb];
            end

            if isa(C.simulation.input_disturbance,'function_handle')
               % Loop over all time instances of simulation, and recalculate
               % the disturbance (note that this disturbance can vary within the matlab "ode" functions, thus we cannot simply record the values at each stage. I.e. it is not piecewise constant)
               fprintf(' - archiving input signals... ')
               archive_time = tic;
               for i = 1:length(C.archive.simulations{end}.sim.time)
                  t = C.archive.simulations{end}.sim.time(i);
                  x = x_sim(:,i);
                  C.archive.simulations{end}.sim.input_disturbance(:,i) = C.simulation.input_disturbance(C,t,x);
               end
               C.archive.simulations{end}.disturbance_archive_time = toc(archive_time);
               disp(['done. ',sec2str(C.archive.simulations{end}.disturbance_archive_time)])
            elseif isa(C.simulation.input_disturbance,'string') && C.simulation.input_disturbance ~= "none"
               error('input_disturbance is neither "none" nor a funciton_handle... ')
            end
            
            % Either the stage-wise control-signals are already recorded,
            % and we simply need to add the very first control signal at
            % the start,,,
            if ~isinf(C.simulation.sampling_time)
               C.archive.simulations{end}.sim.input_controller = [C.archive.simulations{end}.stage(1).input_controller C.archive.simulations{end}.sim.input_controller];
            else
               % ...or we need to loop over all time instances of
               % simulation, to get the non-piece-wise control signal 
               % (in the case of "sampling_time" = inf, the controller re-evaluated at each time instance within the "ode" function, and we need to reconstruct the signal similarly to what we do for the disturbance)
               fprintf(' - archiving input signals... ')
               archive_time = tic;
               for i = 1:length(C.archive.simulations{end}.sim.time)
                  t = C.archive.simulations{end}.sim.time(i);
                  x = x_sim(:,i);
                  C.archive.simulations{end}.sim.input_controller(:,i) = C.simulation.controller_handle(C,t,x);
               end
               C.archive.simulations{end}.input_archive_time = toc(archive_time);
               disp(['done. ',sec2str(C.archive.simulations{end}.input_archive_time)])
            end

            if isa(C.simulation.input_disturbance,'function_handle')
               % Combine controller signal and disturbance to generate the
               % effective input to the system
               C.archive.simulations{end}.sim.input_effective = C.archive.simulations{end}.sim.input_controller + C.archive.simulations{end}.sim.input_disturbance;
            else
               C.archive.simulations{end}.sim.input_effective = C.archive.simulations{end}.sim.input_controller;
            end

            
         end

         % Turn the speed off:
         C.speed = false;

         
      end

   end


   %%% Internal simulation functions
   methods(Access = private,Hidden)

      function log_simulation(C)

         C.archive.simulations{end}.sim.time = [C.archive.simulations{end}.sim.time C.internal_sim.t_sim(2:end)'];
         C.archive.simulations{end}.sim.state = [C.archive.simulations{end}.sim.state C.internal_sim.x_sim(2:end,1:C.dim.state)'];
         if C.flag_has_algeb
            C.archive.simulations{end}.sim.algeb = [C.archive.simulations{end}.sim.algeb C.internal_sim.x_sim(2:end,(C.dim.state+1):end)'];
         end
         if C.flag_has_input && ~isinf(C.simulation.sampling_time)
            %{
            note that u_current and u_controller are both the signal
            produced by the controller, where u_controller holds the most
            recently produced signal until the control delay period is
            over, after which the sigal is applied to the system via
            u_current.
            In conslusion, u:current is the contorller signal that
            should be stored, as it is the actual applied control signal
            (not including disturbance)
            for the current time period.
            %}
               C.archive.simulations{end}.sim.input_controller = [C.archive.simulations{end}.sim.input_controller C.internal_sim.u_current.*ones(1,length(C.internal_sim.t_sim)-1)];          
         end

         % Log current stage-info
         C.archive.simulations{end}.stage(end+1).time             = C.internal_sim.time;
         C.archive.simulations{end}.stage(end).state              = C.internal_sim.state;
         if C.flag_has_input && ~isinf(C.simulation.sampling_time)
            C.archive.simulations{end}.stage(end).input_controller   = C.internal_sim.u_current;
            if C.simulation.control_delay_scaler
               C.archive.simulations{end}.stage(end).delay_period = C.internal_sim.delay_period;
               C.archive.simulations{end}.stage(end).control_delay = C.internal_sim.control_delay;
            else
               C.archive.simulations{end}.stage(end).control_delay = 0;
            end
         end
         

      end
   end






   %%%% Controllers
   methods(Access = private)
      
      % basic NMPC controller
      function u = controller_NMPC_ipopt(C,t,state)

         C.internal_mpc.initial_state = state;

         %%%% Temporary:
         C.internal_ipopt
         u =  full(C.internal_mpc.solution_ipopt.x(C.cas.problem.decision.ind.input(:,1)));
 
         if ~C.speed
            C.extract_solution_ipopt
            tic;
            C.save_archive_optimization("simulation",C.internal_mpc.solution,"ipopt",C.internal_mpc.solver.stats.return_status,C.internal_mpc.solver.stats.success,C.internal_mpc.solver.stats.iter_count)
            C.internal_mpc.archive_time = toc;
         else
            C.internal_mpc.decision = C.internal_mpc.solution_ipopt.x;
         end
      end
      
      function u = controller_LMPC(C,t,state)
         
      end

      function u = controller_asNMPC(C,t,state)
         
      end

      function u = controller_RTI(C,t,state)
         
      end

      function u = controller_aRTI(C,t,state)
         
      end

      function u = controller_aRTIshp(C,t,state)
         
      end

   end





   %%%% Optimization Algorithms (accessible by users):
   methods
      function solution = solve(C,solver,options)
         arguments
            C 
            solver (1,1) string {mustBeMember(solver,["ipopt","sqp"])} = "ipopt"

            % Options
            options.initial_state (:,1) double = C.initial_state.vec;
            options.initial_guess_primal (:,1) double {mustBeReal} = zeros(C.display.problem.n_decision,1);
            options.initial_guess_dual_eq (:,1) double {mustBeReal} = zeros(C.display.problem.n_equality,1);
            options.initial_guess_dual_in (:,1) double {mustBeReal} = zeros(C.display.problem.n_inequality,1);
            
            options.max_iterations (1,1) {mustBeInteger,mustBePositive}
            
            options.state_ref (:,:) double 
            options.algeb_ref (:,:) double 
            options.input_ref (:,:) double 

            options.display_result (1,1) logical = true;
            
            % Solver specific
            options.ipopt_nlpsol_options (1,1) struct = struct('print_time',0);
         end
         
         C.mpc_start_fresh(options)

         % solve:
         disp(['Solving (',char(solver),') ... ']);
         solution = C.(['solve_',char(solver)])(options);
      end
   end
   methods(Access = private,Hidden)
      function solution = solve_ipopt(C,options)
         arguments
            C
            options (1,1) struct % contains the optional fields form the general "solve()" function
         end
         
         if ~isfield(options.ipopt_nlpsol_options,'ipopt')
            options.ipopt_nlpsol_options.ipopt.print_level = 0;
            options.ipopt_nlpsol_options.ipopt.max_iter = 300;

            % Default values:
            options.ipopt_nlpsol_options.ipopt.tol = 1e-8;
            options.ipopt_nlpsol_options.ipopt.acceptable_tol = 1e-6;
            options.ipopt_nlpsol_options.ipopt.compl_inf_tol = 1e-4;
            options.ipopt_nlpsol_options.ipopt.constr_viol_tol = 1e-4;
            options.ipopt_nlpsol_options.ipopt.dual_inf_tol = 1e-4;

            %%% Some other ipopt-options:
            % options.ipopt_nlpsol_options.ipopt.hessian_approximation      = 'limited-memory';
            % options.ipopt_nlpsol_options.ipopt.limited_memory_update_type = 'bfgs';
            % options.ipopt_nlpsol_options.ipopt.linear_solver              = 'mumps';
            % options.ipopt_nlpsol_options.ipopt.linear_system_scaling      = 'none';

         end

         if isfield(options,'max_iterations')
            options.ipopt_nlpsol_options.ipopt.max_iter = options.max_iterations;
         end
         
         C.internal_mpc.nlpsol_options = options.ipopt_nlpsol_options;

         C.pre_initialize_ipopt
         C.initialize_ipopt
         C.internal_ipopt

         C.extract_solution_ipopt


         % make input trajectory the length of the others, to simplify plot
         % if C.flag_has_input
         %    C.internal_mpc.solution.decision.str.input(:,end+1) = nan(C.dim.input,1);
         % end

         C.save_archive_optimization("optimization",C.internal_mpc.solution,"ipopt",C.internal_mpc.solver.stats.return_status,C.internal_mpc.solver.stats.success,C.internal_mpc.solver.stats.iter_count)
         
         % Convergence Measure:
         C.internal_mpc.decision = C.internal_mpc.solution_ipopt.x;
         [stationarity,equality,inequality] = C.convergence;

         if options.display_result
            disp( '================== IPOPT return message ===========')
            disp(['      -- return status: ',C.internal_mpc.solver.stats.return_status])
            disp(['      --    solve time: ',sec2str(C.internal_mpc.solve_time.total)])
            disp(['      -- N. iterations: ',num2str(C.internal_mpc.solver.stats.iter_count)])
            disp( '   Convergence Measure: ')
            disp(['              -- stationarity: ',num2str(stationarity)]) % C.archive.optimizations{end}.convergence.stationarity(end)
            disp(['              --     equality: ',num2str(equality)]) % C.archive.optimizations{end}.convergence.equality(end)
            disp(['              --   inequality: ',num2str(inequality)]) % C.archive.optimizations{end}.convergence.inequality(end)
            disp( '====================================================')
         end
         
         solution = C.internal_mpc.solution;
      end
   
      function solution = solve_sqp(C,options)
         arguments
            C
            options (1,1) struct % contains the optional fields form the general "solve()" function
         end

         if isfield(options,'max_iterations')
            C.set_SQP_settings("max_N_iterations",options.max_iterations)
         end

         C.initialize_sqp
         solve_time = tic;
         C.internal_sqp
         C.internal_mpc.solve_time.total = toc(solve_time);

         % Extract Solution to use for next initial guess
         solution.decision  = C.cas.problem.decision.retrieve(C.internal_mpc.decision);
         primaldual = C.cas.problem.primaldual.retrieve([C.internal_mpc.decision; C.internal_mpc.lambda; C.internal_mpc.mu]);
         solution.dual_eq = structor.subvec(C.cas.problem.primaldual,primaldual.str.lambda);
         if C.flag_has_inequality
            solution.dual_in = structor.subvec(C.cas.problem.primaldual,primaldual.str.mu);
         end

         % % make input trajectory the length of the others, to simplify plot
         % if C.flag_has_input
         %    solution.decision.str.input(:,end+1) = nan(C.dim.input,1);
         % end



         %%%%%%%%%% Save to archive:
         C.save_archive_optimization("optimization",solution,"sqp", C.SQP_info.return_status,C.SQP_info.success_flag,C.SQP_info.n_iterations)
         C.archive.optimizations{end}.solver_specific.primal_iterations = C.SQP_info.primal_iterations;
         C.archive.optimizations{end}.solver_specific.dual_eq_iterations = C.SQP_info.dual_eq_iterations;
         C.archive.optimizations{end}.solver_specific.dual_in_iterations = C.SQP_info.dual_in_iterations;
         C.archive.optimizations{end}.solver_specific.QP_sol = C.SQP_info.QP_sol;
         C.archive.optimizations{end}.solver_specific.step_size = C.SQP_info.step_size;
         C.archive.optimizations{end}.solver_specific.SQP_settings  = C.SQP_settings;


         % convergence measure:
         [stationarity,equality,inequality] = C.convergence;

         if options.display_result
            disp( '================== SQP return message ===========')
            disp(['      -- return status: ',char(C.archive.optimizations{end}.return_status)])
            disp(['      --    solve time: ',sec2str(C.internal_mpc.solve_time.total)])
            disp(['      -- N. iterations: ',num2str(C.archive.optimizations{end}.n_iterations)])
            disp( '   Convergence Measure: ')
            disp(['              -- stationarity: ',num2str(stationarity)]) % C.archive.optimizations{end}.convergence.stationarity(end)
            disp(['              --     equality: ',num2str(equality)]) % C.archive.optimizations{end}.convergence.equality(end)
            disp(['              --   inequality: ',num2str(inequality)]) % C.archive.optimizations{end}.convergence.inequality(end)
            disp( '====================================================')
         end
      end
   end


   % Internal optimization ruitines
   methods(Access = private)

      %%%%%%%%%%%%%%%% ipopt
      function pre_initialize_ipopt(C)
         C.internal_mpc.S = structor(default_mix="TRYMPC_horizon");
         C.internal_mpc.S.str.equality = C.cas.horizon.constraints.equality.str;
         if C.flag_has_inequality
            C.internal_mpc.S.str.inequality = C.cas.horizon.constraints.inequality.str;
         end
         C.internal_mpc.ipopt_inequalities.expr = C.internal_mpc.S.vec;
         
         % % Verify that the constraints are sorted:
         % spy(jacobian(C.internal_mpc.ipopt_inequalities.expr,C.cas.problem.decision.vec))

         [Arg_constraints,in_fields] = C.create_Arg_symbolic_constraints;
         Arg_constraints.out = C.internal_mpc.ipopt_inequalities.expr;

         C.internal_mpc.ipopt_inequalities.F = casadi.Function('F_ipopt_inequalities',Arg_constraints,in_fields,{'out'});
      end

      function initialize_ipopt(C)

         % decision vector:
         C.internal_mpc.solver_def.x = C.cas.problem.decision.vec;


         %%%%%% Initialize ipopt functions:
         %{
Defining nlpsol takes a lot of time (f.ex. 0.5s), therefore, to save time,
we will change initial conditions to be expressed through the
limits. This way we may update only the limit vector, and not have to
re-define nlpsol every time we get a new initial state (which is every
time...)
Similar tricks van be done for parameters, objective coefficients, bounds, etc.,
but this is less relevant, since these usually don't change that often.
Perhaps the reference trajectory is prone to frequent change in tracking
problems, so maybe I will do this for the reference too at some point.
         %}


         % Prepare Args
         C.create_Arg
         C.internal_mpc.Arg_objective.decision = C.cas.problem.decision.vec;
         C.internal_mpc.Arg_constraints.decision = C.cas.problem.decision.vec;

         % When defining nlpsol, we set initial_state to zero, such that we
         % can instead enforce the initial state through the limits, which
         % is faster to update:
         C.internal_mpc.Arg_constraints.initial_state = zeros(C.dim.state,1);

         %%%%%% Define nlpsol:
         % Objective: (this contains reference, which might be frequently changing. May add ref to the fast update trick later)
         C.internal_mpc.solver_def.f = C.cas.problem.objective.F.call(C.internal_mpc.Arg_objective).out;

         % Inequality Constraint vector  (ipopt style)
         if C.flag_has_inequality
            % stacked version:
            % C.internal_mpc.solver_def.g = [C.cas.problem.equality.F.call(C.internal_mpc.Arg_constraints).out;
            %                                C.cas.problem.inequality.F.call(C.internal_mpc.Arg_constraints).out];

            % Create a vector where the equality and inequality are merged,
            % such that they are sorted by time, rather than stacking them.
            % This creates a diagonal structure of the jacobian w.r.t.
            % primal variables. This may speed up the linear solvers, if
            % using solvers that exploit diagonality.
            % S = structor(default_mix="TRYMPC_horizon");
            % S.str.equality = C.cas.horizon.constraints.equality.str;
            % S.str.inequality = C.cas.horizon.constraints.inequality.str;
            C.internal_mpc.solver_def.g = C.internal_mpc.ipopt_inequalities.F.call(C.internal_mpc.Arg_constraints).out;

            % limits on ipopt's inequality constraint vector (lbg < g < ubg)
            C.internal_mpc.lower_limits =  zeros(C.display.problem.n_equality + C.display.problem.n_inequality,1);
            C.internal_mpc.upper_limits =  zeros(C.display.problem.n_equality + C.display.problem.n_inequality,1);
            if C.flag_has_bounds
               for k = 0:C.horizon.N - ~isfield(C.internal_mpc.S.ind.inequality,['stage_',num2str(C.horizon.N)])
                  for name = string(fieldnames(C.internal_mpc.S.ind.inequality.(['stage_',num2str(k)])))'
                     if isa(C.internal_mpc.S.ind.inequality.(['stage_',num2str(k)]).(name),'double')
                        C.internal_mpc.upper_limits(C.internal_mpc.S.ind.inequality.(['stage_',num2str(k)]).(name)) = inf;
                     else
                        for name_2 = string(fieldnames(C.internal_mpc.S.ind.inequality.(['stage_',num2str(k)]).(name)))'
                           C.internal_mpc.upper_limits(C.internal_mpc.S.ind.inequality.(['stage_',num2str(k)]).(name).(name_2)) = inf;
                        end
                     end
                  end
               end
            end
            % if C.flag_has_terminal_bounds
            %    C.internal_mpc.upper_limits(C.internal_mpc.S.ind.inequality.(['stage_',num2str(C.horizon.N)]).bounds) = inf;
            % end
         else
            C.internal_mpc.solver_def.g = [C.cas.problem.equality.F.call(C.internal_mpc.Arg_constraints).out];

            % limits on ipopt's inequality constraint vector (lbg < g < ubg)
            C.internal_mpc.lower_limits = zeros(C.display.problem.n_equality,1);
            C.internal_mpc.upper_limits = zeros(C.display.problem.n_equality,1);
         end

         % Define nlpsol:
         tic
         C.internal_mpc.solver = casadi.nlpsol('solver','ipopt', C.internal_mpc.solver_def, C.internal_mpc.nlpsol_options);
         disp(['define nlpsol: ',sec2str(toc)])

         C.flag_numerical_values_changed = false;
      end

      function internal_ipopt(C)


         if C.flag_numerical_values_changed
            C.initialize_ipopt
         end
         
         % Set the initial condition:
         C.internal_mpc.lower_limits(1:C.dim.state) = -C.internal_mpc.initial_state;
         C.internal_mpc.upper_limits(1:C.dim.state) = -C.internal_mpc.initial_state;

         % Solve:
         tic;
         C.internal_mpc.solution_ipopt = C.internal_mpc.solver('x0', C.internal_mpc.decision, 'lbg', C.internal_mpc.lower_limits, 'ubg', C.internal_mpc.upper_limits);
         C.internal_mpc.solve_time.total = toc;

      end

      function extract_solution_ipopt(C)
         % Extract Solution to use for next initial guess
         C.internal_mpc.solution.decision  = C.cas.problem.decision.retrieve(full(C.internal_mpc.solution_ipopt.x));
         S_num = C.internal_mpc.S.retrieve(full(C.internal_mpc.solution_ipopt.lam_g));
         C.internal_mpc.solution.dual_eq = structor.subvec(S_num,S_num.str.equality);
         if C.flag_has_inequality
            C.internal_mpc.solution.dual_in = structor.subvec(S_num,S_num.str.inequality);
         end

         % when stacked:
         % C.internal_mpc.solution.dual_eq   = full(C.internal_mpc.solution_ipopt.lam_g(1:C.display.problem.n_equality));
         % C.internal_mpc.solution.dual_in   = full(C.internal_mpc.solution_ipopt.lam_g((C.display.problem.n_equality+1):end));

      end
   
   
      %%%%%%%%%%%%%%%% SQP:
      function initialize_sqp(C)


         % Prepare Args:
         C.create_Arg

         % Prepare archive:
         C.SQP_info.return_status = [];
         C.SQP_info.success_flag  = 0;
         C.SQP_info.n_iterations  = [];
         C.SQP_info.decision      = [];
         C.SQP_info.ref = C.internal_mpc.ref;
         C.SQP_info.primal_iterations = C.internal_mpc.decision;
         C.SQP_info.dual_eq_iterations = C.internal_mpc.lambda;
         C.SQP_info.dual_in_iterations = C.internal_mpc.mu;
         C.SQP_info.QP_sol = {};
         C.SQP_info.step_size = [];
         C.SQP_info.convergence.stationarity = [];
         C.SQP_info.convergence.equality = [];
         C.SQP_info.convergence.inequality = [];
         

      end
      
      function internal_sqp(C)

         progress_tolerance = 10;

         if C.SQP_converged
            C.SQP_info.success_flag = 1;
            C.SQP_info.return_status = "Solve Successful. (No iterations needed)";
         else

            % To monitor progress:
            best_staionarity = inf;
            best_equality    = inf;
            no_prog_stationarity = 0;
            no_prog_equality = 0;
            if C.flag_has_inequality
               best_inequality  = inf;
               no_prog_inequality   = 0;
            else
               no_prog_inequality   = [];
            end
            

            C.SQP_info.return_status = "Maximum number of iterations reached.";
            for iteration = 1:C.SQP_settings.max_N_iterations
               % Find Newton Step (Dx, lambda_next, mu_next):
               C.SQP_Newton_Step

               % Linesearch:
               C.SQP_take_step

               % Log iteration:
               C.SQP_log

               % Check convergence
               if C.SQP_converged
                  C.SQP_info.success_flag = 1;
                  C.SQP_info.return_status = "Solve Successful.";
                  break;
               else % Monitor progress:
                  if best_staionarity <= C.SQP_info.convergence.stationarity(end)
                     no_prog_stationarity = no_prog_stationarity + 1;
                     % dispt('stationarity: (best, current)',[best_staionarity, C.SQP_info.convergence.stationarity(end) , best_staionarity - C.SQP_info.convergence.stationarity(end)])
                     % round((best_staionarity- C.SQP_info.convergence.stationarity(end)) * 1e18 )
                  else
                     best_staionarity = C.SQP_info.convergence.stationarity(end);
                     no_prog_stationarity = 0;
                  end
                  if best_equality <= C.SQP_info.convergence.equality(end)
                     no_prog_equality = no_prog_equality + 1;
                  else
                     best_equality = C.SQP_info.convergence.equality(end);
                     no_prog_equality = 0;
                  end
                  if C.flag_has_inequality
                     if best_inequality <= C.SQP_info.convergence.inequality(end)
                        no_prog_inequality = no_prog_inequality + 1;
                     else
                        best_inequality = C.SQP_info.convergence.inequality(end);
                        no_prog_inequality = 0;
                     end
                  end
                  if min([no_prog_stationarity,no_prog_equality,no_prog_inequality]) >= progress_tolerance
                     C.SQP_info.return_status = ['No progress for ',num2str(progress_tolerance),' iterations...'];
                     break;
                  end
                  % 
                  % dispt('iter: ',iteration)
                  % dispt('stationarity: (best, current)',[best_staionarity, C.SQP_info.convergence.stationarity(end) , best_staionarity - C.SQP_info.convergence.stationarity(end)])
                  % dispt('    equality: (best, current, diff)',[best_equality, C.SQP_info.convergence.equality(end), best_equality - C.SQP_info.convergence.equality(end)])
                  % % [best_equality, C.SQP_info.convergence.equality(end), best_equality - C.SQP_info.convergence.equality(end)]
                  % dispt('  inequality: (best, current)',[best_inequality, C.SQP_info.convergence.inequality(end)])
                  % dispt('no progress: (stationarity) ',no_prog_stationarity)
                  % dispt('no progress:     (equality) ',no_prog_equality)
                  % dispt('no progress:   (inequality) ',no_prog_inequality)

               end
            end
         end

         C.SQP_info.n_iterations = iteration;

      end
   end











   %%%%% SQP methods
   methods(Access = private,Hidden)
      
      %%%% Find Newton Step from QP subproblem:
      function SQP_Newton_Step(C)
         

         %%%%%%%%%% Prepare Objective matrices
         C.internal_mpc.Arg_Lagrangian.decision = C.internal_mpc.decision;
         C.internal_mpc.Arg_objective.decision = C.internal_mpc.decision;
         C.internal_mpc.Arg_constraints.decision = C.internal_mpc.decision;

         % Hessian of Lagrangian w.r.t. primal
         HL = sparse(C.cas.problem.Lagrangian.H.F.call(C.internal_mpc.Arg_Lagrangian).out);
         min_eig = min(eig(HL));
         if min_eig < 1e-8
            HL = HL + speye(C.cas.problem.decision.len)*(abs(min_eig)+1e-8); % ensure positive definiteness
            % dispt('min eig: ',min_eig)
         end

         % Jacobian of objective w.r.t. primal
         Jf = sparse(C.cas.problem.objective.J.F.call(C.internal_mpc.Arg_objective).out);
         



         %%%%%%%%% Equality constraints:
         % jacboian of equality constraints w.r.t. primal
         Aeq =  sparse(C.cas.problem.equality.J.F.call(C.internal_mpc.Arg_constraints).out);
         beq = -sparse(C.cas.problem.equality.F.call(C.internal_mpc.Arg_constraints).out);


         %%%%%%%% Inequality constraints:
         if C.flag_has_inequality
            % Notice the sigs, due to the form h(x) > 0:
            Ain = -sparse(C.cas.problem.inequality.J.F.call(C.internal_mpc.Arg_constraints).out);
            bin =  sparse(C.cas.problem.inequality.F.call(C.internal_mpc.Arg_constraints).out);
         else
            Ain = []; bin = [];
         end



         %%%%%%%%% Solve QP subproblem:
         switch C.SQP_settings.QP_solver
            case {"quadprog_AS","quadprog_IP","quadprog_TR"}
               switch C.SQP_settings.QP_solver
                  case "quadprog_AS"
                     QP_options = C.SQP_settings.quadprog_AS_options;
                  case "quadprog_IP"
                     QP_options = C.SQP_settings.quadprog_IP_options;
                  case "quadprog_TR"
                     QP_options = C.SQP_settings.quadprog_TR_options;
               end
               [C.internal_mpc.QP_sol.Dz,C.internal_mpc.QP_sol.FVAL,C.internal_mpc.QP_sol.EXITFLAG,C.internal_mpc.QP_sol.OUTPUT,C.internal_mpc.QP_sol.LAMBDA] ...
                  = quadprog(HL,Jf',Ain,bin,Aeq,beq,[],[],C.internal_mpc.decision,QP_options);
         end

      end




      %%%% Perform a Linesearch and take an approperiate Newton Step:
      function SQP_take_step(C)
         a = 1; % initial step length

         % just for ease of writing the code:
         z  = C.internal_mpc.decision; 
         Dz = C.internal_mpc.QP_sol.Dz;
         lambda      = C.internal_mpc.lambda;
         lambda_next = C.internal_mpc.QP_sol.LAMBDA.eqlin;
         if C.flag_has_inequality
            mu          = C.internal_mpc.mu;
            mu_next     = C.internal_mpc.QP_sol.LAMBDA.ineqlin;
         else
            mu      = -inf;
            mu_next = -inf;
         end
         
         


         % Find step length: (linesearch)
         switch C.SQP_settings.linesearch_method
            case "none"
            case "backtracking" 
               while 1
                  % new variables if the current step size is used:
                  z_new      = z + a*Dz;
                  lambda_new = (1-a)*lambda + a*lambda_next;
                  mu_new     = (1-a)*mu + a*mu_next;


                  % Find multiplier "ny" for merit function: (based on lagrangian multipliers at the candidate next guess, so has to be recomputed every time)
                  ny = max(abs([lambda_new;mu_new]),[],"all")+1e1; % must be bigger than the biggest lagrangian multiplier at the candidate next guess

                  % Find current merit: (must be recomputed every time because of the updated "ny" muliplier)
                  current_merit = full(C.merit(z,ny));
                  % dispt('Current merit:',current_merit)

                  % Find merit after step a*C.Dx:
                  next_merit = full(C.merit(z_new,ny));
                  % dispt('Next merit:',next_merit)

                  % dispt('Merit decrease:',next_merit-current_merit,'a:',a)
                  if (next_merit < current_merit) || (a*C.SQP_settings.backtracking_rate < C.SQP_settings.backtracking_min_step_size)
                     break;
                  end
                  a = a*C.SQP_settings.backtracking_rate;
               end

            case "best_of_N"
               step = nan(2,C.SQP_settings.best_of_N);
               step(1,:) = (1:C.SQP_settings.best_of_N)/C.SQP_settings.best_of_N;


               for i = 1:C.SQP_settings.best_of_N
                  
                  % step length
                  a = step(1,i);

                  % new variables if the current step size is used:
                  z_new      = z + a*Dz;
                  lambda_new = (1-a)*lambda + a*lambda_next;
                  mu_new     = (1-a)*mu + a*mu_next;

                  % Find multiplier "ny" for merit function: (based on lagrangian multipliers at the candidate next guess, so has to be recomputed every time)
                  ny = max(abs([lambda_new;mu_new]),[],"all")+10^-8; % must be bigger than the biggest lagrangian multiplier at the candidate next guess

                  % Find current merit: (must be recomputed every time because of the updated "ny" muliplier)
                  current_merit = full(C.merit(z,ny));

                  % Find merit after step a*C.Dx:
                  next_merit = full(C.merit(z_new,ny));

                  % record performance of step length:
                  step(2,i) = next_merit - current_merit;
               end

               [~,ind] = min(step(2,:));
               a = step(1,ind);
         end



         % Take Step:
         C.internal_mpc.decision = C.internal_mpc.decision + a*Dz;
         C.internal_mpc.lambda   = (1-a)*lambda + a*lambda_next;
         if C.flag_has_inequality
            C.internal_mpc.mu       = (1-a)*mu + a*mu_next;
         end

         % log step size:
         C.SQP_info.step_size(end+1) = a;

         % delete_me_1 = [z Dz C.internal_mpc.decision];
         % disp('step: [z Dz dec.]')
         % disp(delete_me_1(1:10,:))
         % dispt('a:',a,'Newton step norm:', norm(Dz,2))
      end




      %%%% Default merit function for SQP
      function merit = merit(C,z,ny)

         % prepare objective arguments:
         arg_obj = C.internal_mpc.Arg_objective;
         arg_obj.decision = z;

         % prepare constraint arguments:
         arg_const = C.internal_mpc.Arg_constraints;
         arg_const.decision = z;

         merit =   C.cas.problem.objective.F.call(arg_obj).out...
                 + ny*norm(C.cas.problem.equality.F.call(arg_const).out,1);

         if C.flag_has_inequality
            merit = merit - ny*sum(min(0,C.cas.problem.inequality.F.call(arg_const).out));
         end
      end




      %%%% Log current iteration
      function SQP_log(C)
         C.SQP_info.QP_sol{end+1} = C.internal_mpc.QP_sol;
         C.SQP_info.primal_iterations(:,end+1) = C.internal_mpc.decision;
         C.SQP_info.dual_eq_iterations(:,end+1) = C.internal_mpc.lambda;
         C.SQP_info.dual_in_iterations(:,end+1) = C.internal_mpc.mu;
      end




      %%%% Check if SQP has converged
      function all_converged = SQP_converged(C)
         % Get convergence measures:
         [stationarity,equality,inequality] = C.convergence;
         C.SQP_info.convergence.stationarity(end+1) = stationarity;
         C.SQP_info.convergence.equality(end+1) = equality;
         if C.flag_has_inequality
            C.SQP_info.convergence.inequality(end+1) = inequality;
         end

         flag_Lagrangian_converged = (stationarity <= C.SQP_settings.tolerance_lagrangian);
         flag_equality_converged   = (equality <= C.SQP_settings.tolerance_equality);
         if C.flag_has_inequality
            flag_inequality_converged = (inequality >= -C.SQP_settings.tolerance_inequality);
         else
            flag_inequality_converged = 1;
         end
         
         all_converged = flag_Lagrangian_converged && flag_inequality_converged && flag_equality_converged;
      end
   end












   % General Optimization Methods:
   methods(Access = private,Hidden)


      % Initialize the internal_mpc struct:
      function mpc_start_fresh(C,options)
%{
Relevant options:
            options.initial_state (:,1) double = zeros(C.dim.state,1);
            options.initial_guess_primal (:,1) double {mustBeReal} = zeros(C.display.problem.n_decision,1);
            options.initial_guess_dual_eq (:,1) double {mustBeReal} = zeros(C.display.problem.n_equality,1);
            options.initial_guess_dual_in (:,1) double {mustBeReal} = zeros(C.display.problem.n_inequality,1);

            options.state_ref (:,:) double 
            options.algeb_ref (:,:) double 
            options.input_ref (:,:) double
%}

         % Start fresh:
         C.internal_mpc = struct;

         % add reference:
         for type = C.var_types_notpar
            if isfield(options,[char(type),'_ref']) 
               if numel(options.([char(type),'_ref'])) ~= numel(C.cas.horizon.decision.str.(type))
                  C.usererror(['number of elements of "',char(type),'_ref" is not equal to the number of elements of the ',char(type),' prediction trajectory.'])
               end
               C.internal_mpc.ref.(type) = options.(type);
            else
               C.internal_mpc.ref.(type) = zeros(size(C.cas.horizon.decision.str.(type)));
            end
         end

         % initial guess:
         C.internal_mpc.decision = options.initial_guess_primal;
         C.internal_mpc.lambda   = options.initial_guess_dual_eq;
         C.internal_mpc.mu       = options.initial_guess_dual_in;

         % initial state:
         C.internal_mpc.initial_state = options.initial_state;
      end



      % Convergence measure:
      function [stationarity,equality,inequality] = convergence(C)

         % Generate Numeric Arguments:
         C.create_Arg
         Arg_Lagrangian = C.internal_mpc.Arg_Lagrangian;
         Arg_Lagrangian.decision = C.internal_mpc.decision;
         Arg_constraints = C.internal_mpc.Arg_constraints;
         Arg_constraints.decision = C.internal_mpc.decision;
         
         [stationarity,equality,inequality] = C.internal_convergence(Arg_constraints,Arg_Lagrangian);

      end

      function [stationarity,equality,inequality] = internal_convergence(C,Arg_constraints,Arg_Lagrangian)
         % Evaluate:
         stationarity  = full(max(abs(C.cas.problem.Lagrangian.J.F.call(Arg_Lagrangian).out))); % maximum sationarity violation
         equality      = full(max(abs(C.cas.problem.equality.F.call(Arg_constraints).out))); % maximum equality violation
         if C.flag_has_inequality
            % note that the inequalitites are on the form (h(z) >= 0),
            % meaning that the most negative number is the most severe
            % violation:
            
            % worst_inequality_value = full(min(C.cas.problem.inequality.F.call(Arg_constraints).out));    % maximum inequality violation

            % If the wors value still satisfies the inequality constraints
            % (worst_inequality_value > 0), then we say that the
            % violation is zero:

            % violation = min([0 worst_inequality_value]);

            % We want the convergence measure to be a positive value, that
            % should be as small as possible:
            
            % inequality = -violation;


            %%%%% But we do it all in one line to avoid storing
            %%%%% intermetiate variables:
            inequality = -min([0 full(min(C.cas.problem.inequality.F.call(Arg_constraints).out))]);
         else
            inequality    = 0;
         end
      end

      % Prepare Args for defining casadi constraint Functions:
      function [Arg_objective,in_fields] = create_Arg_symbolic_objective(C)
         Arg_objective = C.cas.objective.weights;
         Arg_objective.param = C.cas.var.param;
         for var_type = C.var_types_notpar
            Arg_objective.([char(var_type),'_ref_trajectory'])  = C.cas.horizon.ref.(var_type);
         end
         Arg_objective.decision = C.cas.problem.decision.vec;
         in_fields = fieldnames(Arg_objective);
      end




      % Prepare Args for defining casadi constraint Functions:
      function [Arg_constraints,in_fields] = create_Arg_symbolic_constraints(C)
         Arg_constraints = struct;
         Arg_constraints.Dt = C.cas.integrator.Dt;
         Arg_constraints.initial_state = C.cas.init.state;
         if C.flag_has_param
            Arg_constraints.param = C.cas.var.param;
         end
         if C.flag_has_bounds
            for bound_type = string(fieldnames(C.cas.bounds))'
               for name = string(fieldnames(C.cas.bounds.(bound_type)))'
                  Arg_constraints.([char(bound_type),'_',char(name)]) = C.cas.bounds.(bound_type).(name);
               end
            end
         end
         if C.flag_has_terminal_bounds
            for bound_type = string(fieldnames(C.cas.terminal_bounds))'
               for name = string(fieldnames(C.cas.terminal_bounds.(bound_type)))'
                  Arg_constraints.([char(bound_type),'_terminal_',char(name)]) = C.cas.terminal_bounds.(bound_type).(name);
               end
            end
         end
         Arg_constraints.decision = C.cas.problem.decision.vec;
         in_fields = fieldnames(Arg_constraints);
      end




      % Prepare Args for calling problem functions:
      function create_Arg(C)
         % Prepare generic arguments:
         Arg = struct;
         Arg.param = C.parameters.vec;


         % prepare Arg for objective evaluatiosn
         Arg_objective = C.quadratic_cost;
         for type = C.var_types_notpar
            Arg_objective.([char(type),'_ref_trajectory']) = C.internal_mpc.ref.(type);
         end
         C.internal_mpc.Arg_objective = mergestructs(Arg_objective,Arg);

         % Prepare arguemnts for constraint evaluation:
         Arg_constraints.Dt = C.horizon.Dt;
         Arg_constraints.initial_state = C.internal_mpc.initial_state;
         if C.flag_has_bounds
            for bound_type = string(fieldnames(C.cas.bounds))'
               for name = string(fieldnames(C.cas.bounds.(bound_type)))'
                  Arg_constraints.([char(bound_type),'_',char(name)]) = C.bounds.(bound_type).(name);
               end
            end
         end
         if C.flag_has_terminal_bounds
            for bound_type = string(fieldnames(C.cas.terminal_bounds))'
               for name = string(fieldnames(C.cas.terminal_bounds.(bound_type)))'
                  Arg_constraints.([char(bound_type),'_terminal_',char(name)]) = C.terminal_bounds.(bound_type).(name);
               end
            end
         end

         C.internal_mpc.Arg_constraints = mergestructs(Arg_constraints,Arg);
         C.internal_mpc.Arg_Lagrangian = mergestructs(Arg_constraints,C.internal_mpc.Arg_objective);
      end
   


      % save optimization info to archive:
      function save_archive_optimization(C,type,solution,solver,return_status,success_flag,n_iterations)
         arguments
            C 
            type (1,1) string {mustBeMember(type,["optimization","simulation"])}
            solution
            solver
            return_status 
            success_flag 
            n_iterations 
         end
         optimization = struct;
         optimization.index   = nan; % this is the index in the cell array it is contained in. Is added below.
         optimization.ID = C.generate_id;
         optimization.solver         = solver; % "ipopt"
         optimization.return_status  = return_status;%C.internal_mpc.solver.stats.return_status;
         optimization.success_flag   = success_flag;%C.internal_mpc.solver.stats.success;
         optimization.n_iterations   = n_iterations;%C.internal_mpc.solver.stats.iter_count;
         optimization.solve_time     = C.internal_mpc.solve_time;
         optimization.integrator     = C.integrator;
         optimization.quadratic_cost = C.quadratic_cost;

         for stage_type = ["def_stage_constraints" "def_terminal_constraints"]
            if isfield( C.restore_cell{end}, stage_type)
               for const_type = ["equality" "inequality"]
                  ind = find(string(C.restore_cell{end}.(stage_type)(1,:)) == const_type,1);
                  if ~isempty(ind)
                     optimization.constraints.(const_type) = C.restore_cell{end}.(stage_type){2,ind};
                  end
               end
            end
         end

         optimization.parameters     = C.parameters;
         optimization.Dt             = C.horizon.Dt;
         optimization.N              = C.horizon.N;
         optimization.T              = C.horizon.T;
         optimization.time           = (0:C.horizon.N).*C.horizon.Dt;
         optimization.ref            = C.internal_mpc.ref;
         optimization.initial_state  = C.internal_mpc.initial_state;
         if C.flag_has_bounds
            optimization.bounds      = C.bounds;
         end
         if C.flag_has_bounds
            optimization.terminal_bounds      = C.terminal_bounds;
         end
         optimization.dual_eq        = solution.dual_eq;
         if C.flag_has_inequality
            optimization.dual_in     = solution.dual_in;
         end
         optimization.decision       = solution.decision;
         optimization.convergence = C.convergence_measure(optimization);
         optimization.solver_specific = struct;

         if type == "optimization"
            C.archive.optimizations{end+1} = optimization;
            C.archive.optimizations{end}.index = length(C.archive.optimizations);
         else
            C.archive.simulations{end}.optimizations{end+1} = optimization;
            C.archive.simulations{end}.optimizations{end}.index = length(C.archive.simulations{end}.optimizations);
         end
      end
   end





   % Evaluate problem functions
   methods
      function arg = Arg_objective(C,optimization)
         arguments
            C 
            optimization (1,1) struct
         end
         arg = optimization.quadratic_cost;
         arg.param = optimization.parameters.vec;
         for type = C.var_types_notpar
            arg.([char(type),'_ref_trajectory']) = optimization.ref.(type);
         end
      end
      function arg = Arg_constraints(C,optimization)
         arguments
            C 
            optimization (1,1) struct
         end
         % Prepare arguemnts for constraint evaluation:
         arg.param = optimization.parameters.vec;
         arg.Dt = optimization.Dt;
         arg.initial_state = optimization.initial_state;
         if C.flag_has_bounds
            for bound_type = string(fieldnames(C.cas.bounds))'
               for name = string(fieldnames(C.cas.bounds.(bound_type)))'
                  arg.([char(bound_type),'_',char(name)]) = optimization.bounds.(bound_type).(name);
               end
            end
         end
         if C.flag_has_terminal_bounds
            for bound_type = string(fieldnames(C.cas.terminal_bounds))'
               for name = string(fieldnames(C.cas.terminal_bounds.(bound_type)))'
                  arg.([char(bound_type),'_terminal_',char(name)]) = optimization.terminal_bounds.(bound_type).(name);
               end
            end
         end
      end
      function [arg,arg_obj,arg_const] = Arg_Lagrangian(C,optimization)
         arguments
            C
            optimization (1,1) struct
         end
         arg_obj   = C.Arg_objective(optimization);
         arg_const = C.Arg_constraints(optimization);
         arg = mergestructs(arg_obj,arg_const);
      end

      % Convergence measure:
      function out = convergence_measure(C,optimization)
         arguments
            C 
            optimization (1,1) struct
         end

         % Generate Numeric Arguments:
         [Arg_Lagrangian,~,Arg_constraints] = C.Arg_Lagrangian(optimization);
         Arg_constraints.decision = optimization.decision.vec;
         Arg_Lagrangian.decision = optimization.decision.vec;

         out = struct;
         [out.stationarity,out.equality,out.inequality] = C.internal_convergence(Arg_constraints,Arg_Lagrangian);
         % out.stationarity = stationarity;
         % out.equality = equality;
         % out.inequality = inequality;
      end
   end










   %%%% Set/Get / Clear (access functions)
   methods

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Set numerical values:
      % These set a flag to indicate that a change has been made, which
      % means that something may have to be recomputed. These are mainly
      % meant to be used inside the optimization algorithms, to save time
      % when simulating closed loop mpc control.
      % The flags must be turned off whenever everything is recomputed,
      % such that the algorithm knows that it can skip those computations
      % until next time they change.
      % OBS: make sure not to rely too heavily on these flags in general
      % funcitonality, because of the difficulity in knowing whether
      % something has reset the flags after a change has been made.
      function set_Dt(C,Dt)
         arguments
            C
            Dt (1,1) double {mustBePositive,mustBeReal}
         end

         C.horizon.Dt = Dt;
         if isfield(C.horizon,'N')
            C.horizon.T = C.horizon.Dt * C.horizon.N;
         end

         C.flag_Dt_changed = true;
         C.horizon_decider = "Dt";
      end

      function set_T(C,T)
         arguments
            C
            T (1,1) double {mustBePositive,mustBeReal}
         end

         C.horizon.T = T;
         if isfield(C.horizon,'N')
            C.horizon.Dt = C.horizon.T / C.horizon.N;
         end

         C.flag_Dt_changed = true;
         C.horizon_decider = "T";
      end
   
      function set.parameters(C,in)
         C.parameters = in;
         C.flag_parameters_changed = true;
      end

      function set.quadratic_cost(C,in)
         C.quadratic_cost = in;
         C.flag_quadratic_cost_changed = true;
      end

      function set.bounds(C,in)
         C.bounds = in;
         C.flag_bounds_changed = true;
      end

      function set.terminal_bounds(C,in)
         C.terminal_bounds = in;
         C.flag_terminal_bounds_changed = true;
      end
      
      function set.ref(C,in)
         C.ref = in;
         C.flag_ref_changed = true;
      end

      function out = get.flag_numerical_values_changed(C)
         out = C.flag_parameters_changed || C.flag_quadratic_cost_changed || C.flag_bounds_changed || C.flag_terminal_bounds_changed || C.flag_ref_changed || C.flag_Dt_changed;
      end

      function set.flag_numerical_values_changed(C,in)
         if ~in
            C.flag_Dt_changed = false;
            C.flag_parameters_changed = false;
            C.flag_quadratic_cost_changed = false;
            C.flag_bounds_changed = false;
            C.flag_terminal_bounds_changed = false;
            C.flag_ref_changed = false;
         end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END


      function set_ref(C,type,trajectory)
         arguments
            C
         end
         arguments(Repeating)
            type (1,1) string {mustBeMember(type,["state","algeb","input"])}
            trajectory (:,:) double
         end
         ' function to define a reference funciton @(t) (...) --> [ ; ; ]
         for i = 1:length(type)

         end
      end



      function set_SQP_settings(C,Field,Value)
         arguments
            C
         end
         arguments(Repeating)
            Field (1,1) string {mustBeMember(Field,["tolerance_lagrangian",...
                                                    "tolerance_equality",...
                                                    "tolerance_inequality",...
                                                    "QP_solver",...
                                                    "quadprog_AS_options",...
                                                    "quadprog_IP_options",...
                                                    "quadprog_TR_options",...
                                                    "max_N_iterations",...
                                                    "merit_function",...
                                                    "linesearch_method",...
                                                    "backtracking_rate",...
                                                    "backtracking_min_step_size",...
                                                    "best_of_N"])}
            Value (1,1)
         end


         for i = 1:length(Field)
            field = Field{i};
            value = Value{i};
            switch field
               case {"tolerance_lagrangian", "tolerance_equality", "tolerance_inequality"}
                  if ~isa(value,'double') || value <= 0 || isinf(value) || isnan(value)
                     C.usererror(['The tolerance: "',char(field),'" must be a positive (non-infinite, non-nan) double. We got: "',class(value),'"'])
                  elseif value > 10
                     warning(['WARNING: you are setting a SQP convergence tolerance: "',char(field),'" to a value greater than 10. (this is may be a very slack tolarance)'])
                  end
               case "QP_solver"
                  QP_solvers = ["quadprog_AS","quadprog_IP","quadprog_TR"];
                  if ~isstring(value) || ~ismember(value,QP_solvers)
                     C.usererror(['To select QP solver, the "QP_sover" field must be set to one of the sollowing strings: ',char(strjoin(QP_solvers,", "))])
                  end
               case {"quadprog_AS_options", "quadprog_IP_options", "quadprog_TR_options"}
                  if ~isa(value,'optim.options.Quadprog')
                     C.usererror(['The quadprog options must be "optim.options.Quadprog" instances. We got "',char(field),'" as a: "',class(value),'"'])
                  end
               case "max_N_iterations"
                  if ~isa(value,'double') || round(value) ~= value || value <= 0
                     C.usererror(['The maximum number of iterations: "',char(field),'" must be a positive integer (double). We got: "',class(value),'"'])
                  end
               case "merit_function"
                  if ~isa(value,'function_handle')
                     C.usererror(['The merit function: "',char(field),'" must be a "function _handle" (i.e. "@(decision,lambda,mu) *expr*). We got: "',class(value),'"'])
                  end
               case "linesearch_method"
                  linesearch_strategies = ["none","backtracking","best_of_N"];
                  if ~isstring(value) || ~ismember(value,linesearch_strategies)
                     C.usererror(['The linesearch method: "',char(field),'" must be a string which is member of ["',char(strjoin(linesearch_strategies,", ")),'"]. We got: "',class(value),'". Set to the empty funciton (@(...) []) to use the default merit function.'])
                  end
               case {"backtracking_rate", "backtracking_min_step_size"}
                  if ~isa(value,'double') || value <= 0 || value >= 1 || isnan(value) || isinf(value)
                     C.usererror(['The backtracking option: "',char(field),'" must be a double in range (0,1) (open interval). We got: "',class(value),'"'])
                  end
               case "best_of_N"
                  if ~isa(value,'double') || value <= 1 || isnan(value) || isinf(value) || ~isint(value)
                     C.usererror(['The best_of_N option must be an integer grater than 1. We got: "',class(value),'"'])
                  end
               otherwise
                  warning(['WARNING while setting SQP_settings: The name "',char(field),'" is not a valied SQP_settings property name.'])
            end

            C.SQP_settings.(field) = value;
         end


      end


      function set_control_delay(C,contorl_delay)
         C.internal_sim.control_delay = contorl_delay*C.simulation.control_delay_scaler;
         if C.internal_sim.control_delay > C.simulation.sampling_time
            % disp(' ')
            warning(['The control delay (',sec2str(C.internal_sim.control_delay),') is greater than the sampling time (',sec2str(C.simulation.sampling_time),'). In a practical setting, this may cause failure. Here; we set the control delay equal to the sampling time.'])
            C.internal_sim.control_delay = min([C.internal_sim.control_delay, C.simulation.sampling_time]);
            C.internal_sim.linelength = 0;
         end
      end


      function clear_archive(C,options)
         arguments
            C
            options.field (1,1) string {mustBeMember(options.field,["simulations","optimization"])}
         end
         if isfield(options,'field')
            C.archive.(field) = {};
         else
            C.archive.simulations = {};
            C.archive.optimizations = {};
         end
      end


      function clear_problem(C)
         C.cas.objective = [];
         C.cas.integrator = [];
         C.cas.constraints = [];
         C.cas.horizon = [];
         C.cas.problem = [];
         
         C.internal_mpc = struct;
         C.internal_sim = struct;

         C.clear_flags
      end
   
      function clear_cas(C,omit_warning)
         arguments
            C
            omit_warning (1,1) logical = false;
         end
         if ~omit_warning
            warning('Externally provided expressions will be lost and not recoverable. F.ex. manually provided stage cost, terminal bounds, etc. will not be possible to restore.')
         end
         C.cas = struct;
         C.clear_flags

         C.restore_cell{end+1} = struct;
      end

      function clear_flags(C)
         C.flag_objective_defined = false;
         C.flag_integrator_defied = false;
         C.flag_stage_constraints_defined = false;
         C.flag_terminal_constraints_defined = false;
         C.flag_horizon_defied = false;
         C.flag_problem_defined = false;
         C.flag_has_bounds = false;
         C.flag_has_terminal_bounds = false;
         C.flag_has_inequality = false;

         C.flag_parameters_changed = true;
         C.flag_quadratic_cost_changed = true;
         C.flag_bounds_changed = true;
         C.flag_terminal_bounds_changed = true;
         C.flag_ref_changed = true;
         C.flag_Dt_changed = true;
      end

      function save(C,filename,save_options)
         arguments
            C 
            filename (1,:) char
            save_options = {};
         end
         C.clear_cas(true)
         C.internal_mpc = struct;
         C.internal_sim = struct;
         fprintf('saving... ')
         save(filename,"C",save_options{:})
         disp('')
         disp(['TRYMPC instance: "',char(C.Name),'" (ID:',char(C.ID),') saved to ',filename,'.mat file. Use C.restore_cas() to restore the CasADi expressions.'])
      end

      function restore(C,index)
         arguments
            C
            index (1,1) double {mustBeInteger, mustBePositive} = 1;
         end
         disp('Restoring CasADi variables, expressions and Functions...')
         for def = string(fieldnames(C.restore_cell{index}))'
            C.(def)(C.restore_cell{index}.(def){:})
         end
         disp('Restoration Complete.')
      end
   end






   % Useful stuff (Usererror, help, etc.)
   methods(Static)
      function usererror(text)
         arguments
            text (1,:) char = ''
         end
         error(['USER ERROR: ',char(text),' -Try TRYMPC.help for help.'])
      end

      function help()
         warning('The help function is not yet implmented. Please poke the developer and tell them to get started with this function, so that others can understand how to use this class.')
      end

      function get_started()

         % Define system:
         def_problem = ['%%%%%%%%%%%%%%%%% Initiate class with variable names:' newline ...
            'C = TRYMPC(''Tester Instance'',...' newline ...
            '    state = ["x" "th" "dx" "dth"], ...' newline ...
            '    input = "ux", ...' newline ...
            '    param = ["L" "g" "mx" "mth"]);' newline newline ...
            '%%%%%%%%%%%%%%%%%  Define dynamics:' newline ...
            's = C.cas.state;' newline ...
            'i = C.cas.input;' newline ...
            'p = C.cas.param;' newline newline ...
            'u_pendulum = 0;' newline newline ...
            '% cart acceleration' newline ...
            'ddx = -(p.L*i.ux + u_pendulum*cos(s.th) + p.L^2*s.dth^2*p.mth*sin(s.th) - p.L*p.g*p.mth*cos(s.th)*sin(s.th))/(p.L*(- p.mth*cos(s.th)^2 + p.mx + p.mth));' newline newline ...
            '% pendulum angular acceleration' newline ...
            'ddth = -(p.mx*u_pendulum+ p.mth*u_pendulum + p.L*p.mth*i.ux*cos(s.th) - p.L*p.g*p.mth^2*sin(s.th) + p.L^2*s.dth^2*p.mth^2*cos(s.th)*sin(s.th) - p.L*p.g*p.mx*p.mth*sin(s.th))/(p.L^2*p.mth*(- p.mth*cos(s.th)^2 + p.mx + p.mth));' newline newline ...
            '% dynamics (dq = f(q,u))' newline ...
            'dynamics = [s.dx;' newline ...
            '            s.dth;' newline ...
            '         	 ddx;' newline ...
            '            ddth];' newline newline ...
            '% Apply to Trympc:' newline ...
            'C.def_dynamics(dynamics)' newline newline ...
            '%%%%%%%%%%%%%%% Define LaTex Display names:' newline ...
            'C.plotting.display_names.state.x   = "$x$";' newline ...
            'C.plotting.display_names.state.th  = "$\theta$";' newline ...
            'C.plotting.display_names.state.dx  = "$\dot{x}$";' newline ...
            'C.plotting.display_names.state.dth = "$\dot{\theta}$";' newline ...
            'C.plotting.display_names.input.ux  = "$u_{x}$";' newline newline ...
            '%%%%%%%%%%%%%%%%%  Objective:' newline ...
            '% Define objective:' newline ...
            'C.def_objective(quadratic=["Q","R","dR"]);' newline ...
            'C.quadratic_cost.Q  = [1 10 0.1 0.1];' newline ...
            'C.quadratic_cost.R  = 1;' newline ...
            'C.quadratic_cost.dR = 1;' newline newline ...
            '%%%%%%%%%%%%%%%%% Define Integrator:' newline ...
            '% C.def_integrator("Explicit Euler","n_increments",1)' newline ...
            '% C.def_integrator("Implicit Euler","n_increments",5)' newline ...
            '% C.def_integrator("ERK4","n_increments",2)' newline ...
            'C.def_integrator("collocation","collocation_polynomial_order",4,  "n_increments",2,   collocation_polynomial_type="legendre")' newline newline ...
            '%%%%%%%%%%%%%%%%% Constraints:' newline ...
            '% Stage' newline ...
            'inequality.sum_of_pos = s.x + s.th + 10;' newline ...
            'C.def_stage_constraints("lower_bounds",["x","ux"],"upper_bounds","ux",inequality=inequality)' newline newline ...
            '% Terminal' newline ...
            'terminal_equality.cart_at_zero = s.x;' newline ...
            'terminal_equality.pendulum_at_zero = s.th;' newline ...
            'terminal_equality.pendulum_at_still = s.dth;' newline ...
            'C.def_terminal_constraint("equality",terminal_equality,"upper_bounds","dx","lower_bounds","dx")' newline newline ...
            '%%%%%%%%%%%%%%%%% horizon:' newline ...
            'horizon_length = 150;' newline ...
            'C.def_horizon(horizon_length,"primaldual","sorted")' newline newline ...
            'C.display_problem(4)' newline ...
            'C.display_sparsity("KKT")'];


         optimize = ['%% Set parameters, bounds, and initial state' newline ...
            'C.parameters.str.L = 1;' newline ...
            'C.parameters.str.g = 9.81;' newline ...
            'C.parameters.str.mx = 5;' newline ...
            'C.parameters.str.mth = 3;' newline newline ...
            'init_state = structor;' newline ...
            'init_state.str.x   = 1;' newline ...
            'init_state.str.th  = 0;' newline ...
            'init_state.str.dx  = 0;' newline ...
            'init_state.str.dth = 0;' newline newline ...
            'C.bounds.lower.x = -0.1;' newline ...
            'C.bounds.lower.ux = -0.2;' newline ...
            'C.bounds.upper.ux = 0.3;' newline newline ...
            'C.terminal_bounds.lower.dx = -0.01;' newline ...
            'C.terminal_bounds.upper.dx = 0.01;' newline newline ...
            '%% Open-Loop optimization (ipopt)' newline ...
            'C.set_T(15);' newline ...
            'sol = C.solve("ipopt","initial_state",init_state.vec);' newline ...
            'C.display_optimization;' newline newline ...
            '%% Open-Loop optimization (SQP)' newline ...
            'C.set_T(15);' newline ...
            'C.set_SQP_settings("max_N_iterations",20,"tolerance_lagrangian",5)' newline ...
            'sol = C.solve("sqp","initial_state",init_state.vec);' newline ...
            'C.display_optimization;'];



         disp('To get started, try optimizing an offset correction maneuver of a cart-pendulum system.')
         disp('Here is some code to get you started:')
         disp(' ')
         disp('% ==========================================================')
         disp('% ============== OFFSET CORRECTION MANEUVER ================')
         disp(' ')
         disp(def_problem)
         disp(optimize)
         disp(' ')
         disp('% ======================== END =============================')
         disp('% ==========================================================')
         disp(' ')
      end

      function vec = var2vec(var)
         % note - var can be different datatypes, such as struct, double
         % array, etc.
         % The idea is that these are different and convenient ways of
         % representing data such as parameters or initial states.
         % This function then converts whatever datatype you use into a
         % vector.
         switch string(class(init_state))
            case "double"
               vec = var;
            case "struct"
               vec = [];
               for name = string(fieldnames(var))'
                  var = [vec var.(name)];
               end
         end
      end
      
      function id = generate_id
         characters = ['A':'Z', 'a':'z', '0':'9'];
         id = string(characters(randi(numel(characters), [1, 6])));
      end
   end







   %%% Display methods:
   methods
      function display_problem(C,print_level)
         arguments
            C
            print_level (1,1) double {mustBeInteger,mustBeInRange(print_level,1,4)} = 2;
         end

         if ~C.flag_problem_defined
            C.usererror('the horizon/problem must be defined before displaying. Use "def_horizon" to define horizon/problem.')
         end

         disp(' ')
         disp( '||----------------------------------------')
         disp( '||-------------- PROBLEM -----------------')
         
         %%% Level 1:
         disp( '|| ')
         disp( '||-------- Level 1: (problem size)')
         disp(['||   -  N.     decision variables: ',num2str(C.display.problem.n_decision)])
         disp(['||         +        (state variables): ',num2str(C.display.problem.n_state)])
         
         if C.flag_has_algeb
         disp(['||         +    (algebraic variables): ',num2str(C.display.problem.n_algeb)])
         end

         if C.flag_has_input
         disp(['||         +        (input variables): ',num2str(C.display.problem.n_input)])
         end
         
         if C.integrator.has_aux
         disp(['||         +    (auxiliary variables): ',num2str(C.display.problem.n_aux)])
         end
         
         disp(['||   -  N.   equality constraints: ',num2str(C.display.problem.n_equality)])
         if C.flag_has_inequality
            disp(['||   -  N. inequality constraints: ',num2str(C.display.problem.n_inequality)])
            disp(['||         +   (total N. constraints): ',num2str(C.display.problem.n_inequality + C.display.problem.n_equality)])
         end
         
         disp( '|| ')
         %%% Level 2:
         if print_level >= 2
            disp( '||-------- Level 2: (integration details)')
            disp(['||   -  shooting approach: ', char(C.integrator.shooting_method)])
            disp(['||   -         integrator: ', char(C.integrator.integrator)])
            disp(['||   -  integration order: ', num2str(C.integrator.order)])
            disp(['||   -      N. increments: ', num2str(C.integrator.n_increments)])
            disp(['||   -     horizon length: ',num2str(C.horizon.N),' (N. stages, excl. 0th)'])
         end


         disp( '|| ')
         %%% Level 3:
         if print_level >= 3
            disp( '||-------- Level 3: (constraint details)')
            disp(['||   -  N.   dynamic constraints: ',num2str(C.display.problem.n_dynamic_constraints)])
            disp(['||         +          (continuation): ',num2str(C.display.problem.n_continuation_constraints)])
            disp(['||         +             (auxiliary): ',num2str(C.display.problem.n_auxiliary_constraints)])
            if C.flag_has_algeb
            disp(['||   -  N. algebraic constraints: ',num2str(C.display.problem.n_algebraic_constraints)])
            end

            if C.flag_stage_constraints_defined && isfield(C.display.problem,'n_other_constraints')
            disp(['||   -  N.     other constraints: ',num2str(C.display.problem.n_other_constraints)])
            end

            eq_names = string(fieldnames(C.display.problem.stage_constraints.n_equality))';
            if ~isempty(eq_names)
                     disp( '||   -  Stage-wise equality constraints: (g() = 0)')
                  for name = eq_names
                     disp(['||         +   N.  ',num2str(C.display.problem.stage_constraints.n_equality.(name)),'    "',char(name),'"'])
                  end
            end

            in_names = string(fieldnames(C.display.problem.stage_constraints.n_inequality))';
            if ~isempty(in_names)
                     disp( '||   -  Stage-wise inequality constraints: (h() >= 0)')
                  for name = in_names
                     disp(['||         +   N.  ',num2str(C.display.problem.stage_constraints.n_inequality.(name)),'    "',char(name),'"'])
                  end
            end
         end

         disp( '|| ')
         %%% Level 4:
         if print_level >= 4
            disp( '||-------- Level 4: (Jacobian details)')

            for name = ["objective","equality","inequality"]
               if ((name ~= "inequality") + (C.flag_has_inequality))
                  disp(['||   -  ',char(name),' Jacobian: '])
                  disp(['||         +      N. non-zeros: ',num2str(C.display.problem.(['J_',char(name),'_nnz']))])
                  disp(['||         +          sparsity: ',num2str(C.display.problem.(['J_',char(name),'_sparsity']))])
                  disp(['||         +    linear density: ',num2str(C.display.problem.(['J_',char(name),'_linear_density']))])
               end
            end
            
         end

         disp( '||----------------------------------------')
         disp(' ')
      end

      function display_sparsity(C,expression)
         arguments
            C
            expression (1,:) string {mustBeMember(expression,["objective","equality","inequality","KKT","Lagrangian"])} = ["equality", "KKT"];
         end

         Tit.objective = "Objective";
         Tit.equality = "Equality";
         Tit.inequality = "Inequality";
         Tit.KKT = "KKT Matrix";
         Tit.Lagrangian = "Lagrangian";

         desc.objective = "(hessian w.r.t. primal)";
         desc.equality = "(jacobian w.r.t. primal)";
         desc.inequality = "(jacobian w.r.t. primal)";
         desc.KKT = "(hessian of Lagrangian w.r.t. primal-dual)";
         desc.Lagrangian = "(hessian w.r.t. primal)";

         Col.objective = 'b';
         Col.equality = 'r';
         Col.inequality = 'g';
         Col.KKT = 'm';
         Col.Lagrangian = 'c';

         figure('Name','Sparsity Pattern ')
         Layout = tiledlayout('flow','TileSpacing','compact','Padding','compact');
         title(Layout,'Sparsity');

         for expr_type = expression

            switch expr_type
               case {"equality", "inequality"}
                  expr = C.cas.problem.(expr_type).J.expr;
               case {"objective", "Lagrangian"}
                  expr = C.cas.problem.(expr_type).H.expr;
               case "KKT"
                  expr = C.cas.problem.(expr_type).matrix.expr;
            end
            tile = nexttile;
            spy(expr,Col.(expr_type))
            title(tile,[char(Tit.(expr_type)),' ',char(desc.(expr_type))])
         end
      end
   
      function [tiles,Layout] = display_simulation(C,options)
         arguments
            C

            % unique options:
            options.simulation_number (1,:) double {mustBeInteger,mustBePositive} = length(C.archive.simulations);

            % common options:
            options.tiledlayout_varargin (1,:) cell;
            options.gridstyle (1,1) string {mustBeMember(options.gridstyle,["flow","vertical","horizontal"])} = "flow";
            options.state (1,:) string
            options.algeb (1,:) string
            options.input (1,:) string
            options.time_order (1,1) string {mustBeMember(options.time_order,["seconds","minutes","hours","days","weeks","months","years"])} = "seconds";
            options.mark_samples (1,1) logical = true;
            options.multiplot (1,1) string {mustBeMember(options.multiplot,["separate","ontop"])} = "ontop";
            options.transparency (1,1) double {mustBeInRange(options.transparency,0,1,"exclude-lower")} = 0.7;
            options.legend (1,:) string
            options.linestyle (1,:) string
            options.linewidth (1,:) string
            options.colors (1,:) % either a cell array of bgr triplets or string array for "GetCOlorCode"
            options.color_match(1,1) string {mustBeMember(options.color_match,["plot","solution"])} = "plot";
         end

         % create figure object to plot in
         figure_text = [char(C.Name), ' : simulation '];
         figure(Name=figure_text)

         % Prepare title text:
         ID_text = append("ID:",C.ID," - sim-ID(index): ");
         contr_text = "Controller: ";
         sampl_text = "Sampling Time: ";
         cd_text = "control delay scaler: ";

         % Prepare data
         data = dictionary;

         for i = options.simulation_number
            ID_text = append(ID_text,"/",C.archive.simulations{i}.ID,"(",string(i),")");
            contr_text = append(contr_text,"/",strrep(C.simulation.controller,'_','\_'));
            sampl_text = append(sampl_text,"/",num2str(C.archive.simulations{i}.sampling_time)); 
            cd_text = append(cd_text,"/",num2str(C.archive.simulations{i}.control_delay_scaler));

            data(i) = C.archive.simulations{i}.sim;
         end

         title_text = [ID_text ; contr_text; sampl_text];


         % Number:
         options.display_number = options.simulation_number;

         [tiles,Layout] = display_trajectories(C,...
            "simulation",...
            data,...
            title_text,...
            options);



      end


      function [tiles,Layout] = display_optimization(C,options)
         arguments
            C

            % unique options:
            options.optimization_number (1,:) double {mustBeInteger,mustBePositive} = length(C.archive.optimizations);

            % common options:
            options.tiledlayout_varargin (1,:) cell
            options.gridstyle (1,1) string {mustBeMember(options.gridstyle,["flow","vertical","horizontal"])} = "flow";
            options.state (1,:) string
            options.algeb (1,:) string
            options.input (1,:) string
            options.time_order (1,1) string {mustBeMember(options.time_order,["seconds","minutes","hours","days","weeks","months","years"])} = "seconds";
            options.mark_samples (1,1) logical = true;
            options.multiplot (1,1) string {mustBeMember(options.multiplot,["separate","ontop"])} = "ontop";
            options.transparency (1,1) double {mustBeInRange(options.transparency,0,1,"exclude-lower")} = 0.7;
            options.legend (1,:) string
            options.linestyle (1,:) string
            options.linewidth (1,:) string
            options.colors (1,:) % either a cell array of bgr triplets or string array for "GetCOlorCode"
            options.color_match(1,1) string {mustBeMember(options.color_match,["plot","solution"])} = "plot";
         end
 


         % create figure object to plot in
         figure_text = [char(C.Name), ' : optimization'];
         figure(Name=figure_text)

         % Prepare title text:
         ID_text = append("ID:",C.ID," - opt-ID(index): ");
         solver_text = "Solver (n.iter/sol.time): ";
         integr_text = "Integrator (order/incr): ";
         return_text = "Return Status: ";

         % Prepare data
         data = dictionary;

         for i = options.optimization_number
            ID_text = append(ID_text," / ",C.archive.optimizations{i}.ID,"(",string(i),")");
            solver_text = append(solver_text," / ",C.archive.optimizations{i}.solver,"(",string(C.archive.optimizations{i}.n_iterations),"/",sec2str(C.archive.optimizations{i}.solve_time.total),")");
            integr_text = append(integr_text," / ",C.archive.optimizations{i}.integrator.integrator,"(",string(C.archive.optimizations{i}.integrator.order),"/",string(C.archive.optimizations{i}.integrator.n_increments),")"); 
            return_status = strrep(C.archive.optimizations{i}.return_status,'_','\_');
            return_text = append(return_text," / ",return_status);

            % data(i) = C.archive.optimizations{i}.sim;
            data(i) = C.archive.optimizations{i}.decision.str;
            data(i).time = C.archive.optimizations{i}.time;
         end

         title_text = [ID_text ; solver_text; integr_text; return_text];


         % Number:
         options.display_number = options.optimization_number;

         [tiles,Layout] = display_trajectories(C,...
            "optimization",...
            data,...
            title_text,...
            options);


      end

      function [tiles,Layout] = display_constraints(C,options)
         arguments
            C

            options.optimization_number (1,:) double {mustBeInteger,mustBePositive} = length(C.archive.optimizations);

            options.tiledlayout_varargin (1,:) cell
            options.gridstyle (1,1) string {mustBeMember(options.gridstyle,["flow","vertical","horizontal"])} = "flow";
            options.equality (1,:) string % specify what equality constraints to plot
            options.inequality (1,:) string % specify what inequality constraints to plot
            options.time_order (1,1) string {mustBeMember(options.time_order,["seconds","minutes","hours","days","weeks","months","years"])} = "seconds";
            options.mark_samples (1,1) logical = true;
            options.multiplot (1,1) string {mustBeMember(options.multiplot,["separate","ontop"])} = "ontop";
            options.transparency (1,1) double {mustBeInRange(options.transparency,0,1,"exclude-lower")} = 0.7;
            options.legend (1,:) string
            options.linestyle (1,:) string
            options.linewidth (1,:) string
            options.colors (1,:) % either a cell array of bgr triplets or string array for "GetCOlorCode"
            options.color_match(1,1) string {mustBeMember(options.color_match,["plot","solution"])} = "plot";
         end

         % create figure object to plot in
         figure_text = [char(C.Name), ' : constraints'];
         figure(Name=figure_text)

         % Prepare title text:
         ID_text = append("ID:",C.ID," - opt-ID(index): ");
         solver_text = "Solver (n.iter/sol.time): ";
         integr_text = "Integrator (order/incr): ";
         return_text = "Return Status: ";

         % Prepare data
         data = dictionary;

         for i = options.optimization_number
            ID_text = append(ID_text," / ",C.archive.optimizations{i}.ID,"(",string(i),")");
            solver_text = append(solver_text," / ",C.archive.optimizations{i}.solver,"(",string(C.archive.optimizations{i}.n_iterations),"/",sec2str(C.archive.optimizations{i}.solve_time.total),")");
            integr_text = append(integr_text," / ",C.archive.optimizations{i}.integrator.integrator,"(",string(C.archive.optimizations{i}.integrator.order),"/",string(C.archive.optimizations{i}.integrator.n_increments),")"); 
            return_status = strrep(C.archive.optimizations{i}.return_status,'_','\_');
            return_text = append(return_text," / ",return_status);

            data(i) = C.archive.optimizations{i}.decision.str;
            data(i).time = C.archive.optimizations{i}.time;
         end

         title_text = [ID_text ; solver_text; integr_text; return_text];


         % Number:
         options.display_number = options.optimization_number;

         [tiles,Layout] = display_trajectories(C,...
            "constraints",...
            data,...
            title_text,...
            options);

      end
   end


   %%% Internal display methods
   methods(Access = private)
      function [tiles,Layout] = display_trajectories(C,display_type,data,title_text,options)
         arguments
            C
            % unique:
            display_type (1,1) string {mustBeMember(display_type,["simulation","optimization","constraints"])}
            data (1,1) dictionary
            title_text string

            options

            % Common optoins:
            % options.tiledlayout_varargin (1,:) cell;
            % options.gridstyle (1,1) string {mustBeMember(options.gridstyle,["flow","vertical","horizontal"])} = "flow";
            % options.state (1,:) string
            % options.algeb (1,:) string
            % options.input (1,:) string
            % (options.equality (1,:) string)
            % (options.inequality (1,:) string)
            % options.time_order (1,1) string {mustBeMember(options.time_order,["seconds","minutes","hours","days","weeks","months","years"])} = "seconds";
            % options.mark_samples (1,1) logical = true;
            %options.transparency (1,1) double {mustBeInRange(options.transparency,0,1,"exclude-lower")} = 0.7;
            % options.legend (1,:) string
            % options.linestyle (1,:) string
            % options.linewidth (1,:) string
            % options.colors (1,:) % either a cell array of bgr triplets or string array for "GetCOlorCode"
            % options.color_match(1,1) string {mustBeMember(options.color_match,["plot","solution"])} = "plot";
         end


         % Tiled layout to structure the figure
         if isfield(options,'tiledlayout_varargin')
            tiled_varargin = options.tiledlayout_varargin;
         else
            tiled_varargin = {options.gridstyle,"TileSpacing","compact","Padding","compact"};
         end
         Layout = tiledlayout(tiled_varargin{:});
         title(Layout,title_text,'FontSize',11,'FontWeight','bold')


         Types = [];
         if display_type == "constraints"
            for const_type = ["equality", "inequality"]
               if isfield(options,const_type)
                  Types = [Types const_type]; %#ok<AGROW>
               end
            end
            if isempty(Types)
               Types = ["equality", "inequality"];
            end
         else
            for type = C.var_types_notpar
               if isfield(options,type)
                  Types = [Types type]; %#ok<AGROW>
               end
            end
            if isempty(Types)
               Types = C.var_types_notpar;
            end
         end


         % Time order
         switch options.time_order
            case "seconds"
               time_scaler = @(t) t;
            case "minutes"
               time_scaler = @(t) t/60;
            case "hours"
               time_scaler = @(t) t/(60*60);
            case "days"
               time_scaler = @(t) t/(60*60*24);
            case "weeks"
               time_scaler = @(t) t/(60*60*24*7);
            case "months"
               time_scaler = @(t) t/(60*60*24*30.44);
            case "years"
               time_scaler = @(t) t/(60*60*24*365.25);
         end


         % Prepare display names (legend):
         if isfield(options,'legend')
            if length(options.legend) ~= length(options.display_number)
               C.usererror('the "legend" argument must be a string array with the same length as the number of solutions to display (siulation_/optimization_number)')
            end
            disp_names = options.legend;
         else
            disp_names = string(options.display_number);
         end


         % hold all axes handles in a struct, and pass as output so that
         % one can edit them outside this function.
         tiles = struct([]);

         % Plot:
         L = length(options.display_number);
         for i = 1:L
            for type = Types
               if ~isfield(options,type) || (isscalar(options.(type)) && options.(type) == "all")
                  if display_type == "constraints"
                     names = [string(fieldnames(C.archive.optimizations{options.display_number(i)}.constraints.equality))',...
                              string(fieldnames(C.archive.optimizations{options.display_number(i)}.constraints.inequality))']; %#ok<*AGROW>
                  else
                     names = C.names.code.(type)';
                  end
               else
                  names = options.(type); %#ok<*PROPLC>
               end

               for name = names
                  if i == 1 || options.multiplot == "separate"% on first pass, create new axes for every plot:
                     tiles(i).(type).(name) = nexttile(Layout);
                     tile = tiles(i).(type).(name);
                     hold(tile,"on");
                     grid(tile,"on");

                     if display_type == "constraints"
                        ylabel(tile,strrep(name,"_","\_"));
                     else
                        ylabel(tile,C.plotting.display_names.(type).(name),Interpreter="latex");
                     end
                     
                     if L > 1
                        legend(tile)
                     end

                     if type == "state"
                        tile.Color = [0.99, 1, 0.96];
                     elseif type == "algeb"
                        tile.Color = [1, 1, 0.96];
                     elseif type == "input"
                        tile.Color = [1, 0.98, 0.99];
                     elseif type == "equality"
                        tile.Color = [0.97, 0.97, 0.99];
                     elseif type == "inequality"
                        tile.Color = [1, 0.95, 0.95];
                     end

                  elseif options.multiplot == "ontop" % for the remaining passes, reuse old axes:
                     tile = tiles(1).(type).(name);
                  else
                     error('DEVELOPER ERROR: not "ontop" nor "separate"... ')
                  end
                  
                  % Prepare data:
                  if display_type == "constraints"
                     Data = data(options.display_number(i));
                     D.algeb = struct;
                     D.input = struct;
                     D.param = C.archive.optimizations{options.display_number(i)}.parameters.str;

                     variable = nan(1,length(Data.time));
                     for j = 1:(length(Data.time)-1)
                        for variable_type = C.var_types_notpar
                           for jj = 1:C.dim.(variable_type)
                              D.(variable_type).(C.names.code.(variable_type)(jj)) = Data.(variable_type)(jj,j);
                           end
                        end
                        variable(j) = C.archive.optimizations{options.display_number(i)}.constraints.(type).(name)(D.state,D.algeb,D.input,D.param);
                     end
                  else
                     if type == "input" && display_type == "simulation"
                        type_2 = "input_effective";
                     else
                        type_2 = type;
                     end
                     variable = data(options.display_number(i)).(type_2)(C.names.ind.(type).(name),:);
                     if type == "input" && display_type ~= "simulation"
                        variable(:,end+1) = nan; % make inputs the same length as the time vector
                     end
                  end

                  % Time:
                  time = time_scaler(data(options.display_number(i)).time);

                  % Manage color
                  if isfield(options,'colors')
                     if isa(options.colors,"string")
                        color = GetColorCode(options.colors(i));
                     elseif isa(options.colors,"double") && length(options.colors{i}) == 3
                        color = options.colors{i};
                     else
                        C.usererror('The colors arguments must either be a string array of arguments for the "GetColorCode" function, or be a cell array of RGB triplets.')
                     end
                  elseif options.color_match == "solution"
                     color = GetColorCode(i);
                  elseif options.color_match == "plot"
                     if display_type == "constraints"
                        if type == "equality"
                           pre_color = GetColorCode('b',1.3);
                        elseif type == "inequality"
                           pre_color = GetColorCode('r');
                        else
                           error('DEVELOPER ERROR: neither "equality" nor "inequality"...')
                        end
                        
                     else
                        pre_color = C.plotting.color.(type).(name);
                     end
                     if L == 1 || options.multiplot == "separate"
                        color = pre_color;
                     elseif options.multiplot == "ontop"
                        color = interp1([0 (L+1)/2 L+1],[0 0 0; pre_color; 1 1 1],i);
                     end
                  end

                  % Menage linestyle (if externally provided linestyle, it is on a solution basis, if not, the linestyle is dictated by the variable preference)
                  if isfield(options,'linestyle')
                     linestyle = options.linestyle(i);
                  else
                     if display_type == "constraints"
                        linestyle = '-';
                     else
                        linestyle = C.plotting.linestyle.(type).(name);
                     end
                  end

                  % Menage linewidth (if externally provided linewidth, it is on a solution basis, if not, the linewidth is dictated by the variable preference)
                  if isfield(options,'linewidth')
                     linewidth = options.linewidth(i);
                  else
                     if display_type == "constraints"
                        linewidth = 1;
                     else
                        linewidth = C.plotting.linewidth.(type).(name);
                     end
                  end
                  

                  if display_type == "constraints"
                     plot(tile,[0 time(end)],[0 0],Color=GetColorCode('e'),LineStyle='--',HandleVisibility='off')
                  end


                  plot(tile,...
                       time,...
                       variable,...
                       Color=[color options.transparency],...
                       LineStyle=linestyle,...
                       LineWidth=linewidth,...
                       DisplayName=disp_names(i));




                  % Mark state sampling times
                  if display_type == "simulation" && ~isinf(C.archive.simulations{options.display_number(i)}.sampling_time) && options.mark_samples && ismember(type,["state","algeb"])
                     if C.simulation.control_delay_scaler
                        inc = 2;
                     else
                        inc = 1;
                     end
                     samp_time = [C.archive.simulations{options.display_number(i)}.stage(1:inc:end).time];
                     samp_values = [C.archive.simulations{options.display_number(i)}.stage(1:inc:end).(type)];
                     samp_values = samp_values(C.names.ind.(type).(name),:);
                     plot(tile,samp_time,samp_values,LineStyle='none',Marker='square',MarkerSize=4,MarkerEdgeColor=color,HandleVisibility='off')
                  end

               end
            end
         end

         xlabel(['Time [',char(options.time_order),']'])
      end

   end
end