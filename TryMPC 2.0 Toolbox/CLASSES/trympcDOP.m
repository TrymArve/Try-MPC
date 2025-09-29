classdef trympcDOP
%{
This class takes a discretized model and produces a Dynamic Optimization
Problem (DOP). This DOP is formulated as a nonlinear program (NLP), and is the optimization problem used in the
following MPC simulations.


INFO:

- inequalities are on form: h(x) >= 0
%}

   properties(SetAccess=immutable)
      Name
      discretizer
      N_horizon % Number of steps on the horizon
   end
   properties(Hidden,SetAccess=immutable)
      ID (1,1) string

      % Other properties
      rDT
      args struct = struct; % a struct that contains the fields required (by default) to evaluate the various casadi-functions (objective, constraints, ...)
   end
   properties(SetAccess=private)
      initial_condition casadi.SX % a variable that represent the initial condition (a parameter)
      T_horizon casadi.SX % the total length of the horizon (in time)
      decision structor = structor("default_mix","TRYMPC_horizon");
      constraints struct = struct;
      objective struct = struct;
      reference struct = struct;
   end



   methods
      function C = trympcDOP(dop_name,discretizer,options)
         arguments
            dop_name (1,1) string
            discretizer (1,1) trympcDISCRETIZER
            
            %%%%%%%% Main Options:

            %%% Select time discretization grid. Must choose ONE of these:
            options.rDT (1,:) double {mustBePositive}
            options.N_horizon (1,1) double {mustBePositive,mustBeInteger} % set the number of integration steps
            %{
Priority:
rDT - ("relative DT") Use this to set each integration step individually
(N_horizon is implied). Each timestep (Dt) is given relative to the full
horizon. The full horizon is given externally for each solve of the NLP.
Each time step i has length = rDT(i)/sum(rDT)*T_horizon (where T_horizon is the horizon given in time)
Setting rDT implies the number of integration steps.

Secondary options:
N_horizon - The number of integration steps. Setting this instead of rDT
means that all time steps are equally long.
            %}


            %%% Add a stage cost: (choose one of these)
            options.quad_Q (:,1) double {mustBeReal} % quadratic cost for states
            options.quad_Z (:,1) double {mustBeReal} % quadratic cost for algebraic variables
            options.quad_R (:,1) double {mustBeReal} % quadratic cost for inputs
            options.quad_Y (:,1) double {mustBeReal} % quadratic cost for outputs
%{
These are the quadratic cost matrices for a single stage, thus:
  - Q \in R^(n_state \times n_state)
  - Z \in R^(n_algeb \times n_algeb)
  - R \in R^(n_input \times n_input)
  - Y \in R^(n_output \times n_output)

Note that the matrices are given as vectors of the diagonal elements, as no
cross-terms are included. (cross terms can be manually added if needed.)
%}


            %%% add variable bounds
            options.bounds_state (:,2) double {mustBeReal} % [lower upper] bound on state vector
            options.bounds_algeb (:,2) double {mustBeReal} % [lower upper] bound on algeb vector
            options.bounds_input (:,2) double {mustBeReal} % [lower upper] bound on input vector
            options.bounds_output(:,2) double {mustBeReal} % [lower upper] bound on output vector

            %%%%%%%%% Other:
            options.include_algebraic_variables_at_stages (1,1) logical = false; % note that this is not necessary for integration/prediction purposes, just if one wants to add constraints or sumthin' to 'em...
         end

         fprintf('Creating Dynamic Optimization Problem (DOP)...  ');
         def_time = tic;

         % Set Properties
         C.Name = dop_name;
         C.ID = TRYMPC2.generate_id;
         C.discretizer = discretizer;
         if isfield(options,'rDT')
            C.rDT = options.rDT./sum(options.rDT);
         else
            C.rDT = ones(1,options.N_horizon)./options.N_horizon;
         end
         C.N_horizon = length(C.rDT);
         C.T_horizon = casadi.SX.sym('T_horizon');
         DT = C.rDT.*C.T_horizon;





         %% Create Decision variables:

         % Stage Variables:
         C.decision.str.state = casadi.SX.sym('decision_state',[C.discretizer.model.dim.state,C.N_horizon]); % contains state_0 -> state_(N-1)
         if options.include_algebraic_variables_at_stages
            C.decision.str.algeb = casadi.SX.sym('decision_algeb',[C.discretizer.model.dim.algeb,C.N_horizon]); % contains algeb_1 -> algeb_(N)
         end
         C.decision.str.input = casadi.SX.sym('decision_input',[C.discretizer.model.dim.input,C.N_horizon]); % contains input_0 -> input_(N-1)
         
         % Reference variables
         C.reference.state = casadi.SX.sym('reference_state',[C.discretizer.model.dim.state,C.N_horizon]); % contains state_0 -> state_(N-1)
         if options.include_algebraic_variables_at_stages
            C.reference.algeb = casadi.SX.sym('reference_algeb',[C.discretizer.model.dim.algeb,C.N_horizon]); % contains algeb_1 -> algeb_(N)
         end
         C.reference.input = casadi.SX.sym('reference_input',[C.discretizer.model.dim.input,C.N_horizon]); % contains input_0 -> input_(N-1)
         if ismember(C.discretizer.model.expression_types,"output")
            C.reference.output = casadi.SX.sym('reference_output',[C.discretizer.model.dim.output,C.N_horizon]); % contains state_0 -> state_(N-1)
         end

         %{
Note that the last state state_N is not part of the decision variables.
(the terminal state is state_N = Next(state_(N-1),input_(N-1)), but it doesn't get it own variable)
Of course, the inputs are neither present at the terminal stage.

We emphazise that the all algebraic variables that are necessary for
propogating the system are present in the "aux" variables already, and are
used for evaluating the dynamics inside the discretizer. Now new algebraic
variables are necessary at the stages of the DOP. Though if you want to
add constraints of costs associated with them, it might be nice to predict
them too. There this option is added.
         %}

         % % Create auxiliary variables for each stage.
         % i = 1;
         % aux_matrix = casadi.SX.sym(['decision_aux_step_',num2str(i)],[C.discretizer.aux_var.len,1]); % (for ease of using them after this)
         % for i = 1:C.N_horizon
         %    aux_vars = casadi.SX.sym(['decision_aux_step_',num2str(i)],[C.discretizer.aux_var.len,1]);
         %    aux_matrix = [aux_matrix aux_vars]; %#ok<AGROW>
         %    aux_vars = C.discretizer.aux_var.retrieve(aux_vars);
         %    C.decision.str.aux.(['step_',num2str(i)]) = aux_vars.str;
         % end
         % aux_matrix = aux_matrix(:,2:end);

         % Create auxiliary variables for each stage.
         aux_matrix = casadi.SX.zeros([C.discretizer.aux_var.len,C.N_horizon]);
         for i = 1:C.N_horizon
            aux_vars = casadi.SX.sym(['decision_aux_step_',num2str(i)],[C.discretizer.aux_var.len,1]);
            aux_matrix(:,i) = aux_vars;
            aux_vars = C.discretizer.aux_var.retrieve(aux_vars);
            C.decision.str.aux.(['step_',num2str(i)]) = aux_vars.str;
         end



         %% Dynamic constraints:

         %{
 ensure 
   ---- state_(k+1) = discretizer.Next(state_k,input_k,param,aux_var,Dt)
and
   ---- 0 = discretizer.Algebraics(state_k,input_k,param,aux_var,Dt)

at each stage k \in {0,...,N-1}
         %}

         % Prepare structor and arguments
         args = struct;
         C.constraints.equality = structor("default_mix","vector_mixed");

         % initial condition constraint:
         C.initial_condition = casadi.SX.sym('initial_condition',[C.discretizer.model.dim.state,1]);
         C.constraints.equality.str.initial_condition = C.decision.str.state(:,1) - C.initial_condition;

         for k = 0:(C.N_horizon-1) % two because the last state does not have a varaible, thus no shooting gap at the end
            args.state = C.decision.str.state(:,k+1);
            args.input = C.decision.str.input(:,k+1);
            args.aux_var = aux_matrix(:,k+1);
            args.param = C.discretizer.model.param.vec;
            args.Dt = DT(k+1);
            if ~isempty(C.discretizer.Algebraics)
               C.constraints.equality.str.(['step_',num2str(k+1)]).algebraic = C.discretizer.Algebraics.call(args).out;
            end
            if k <= C.N_horizon-2
               C.constraints.equality.str.(['step_',num2str(k+1)]).dynamic   = C.decision.str.state(:,k+2) + C.discretizer.Next.call(args).out;
            end
            
         end

         % keep track of fields required:
         C.args.equality.decision = [];
         C.args.equality.initial_condition = [];
         C.args.equality.T_horizon = [];
         C.args.equality.parameters = [];
         


         %% Objective
         quad_cost = struct;

         % state
         if isfield(options,'quad_Q')
            if size(options.quad_Q,1) ~= C.discretizer.model.dim.state
               TRYMPC2.usererror('The quadratic cost vector (quad_Q) must have the same length as the state vector.')
            end
            quad_cost.state = casadi.SX.zeros(C.discretizer.model.dim.state,C.N_horizon);
            for k = 1:C.N_horizon
               stage_deviation = C.decision.str.state(:,k) - C.reference.state(:,k);
               quad_cost.state(:,k) = options.quad_Q.*(stage_deviation.^2);
            end
         end

         % algeb
         if isfield(options,'quad_Z')
            if ~options.include_algebraic_variables_at_stages
               TRYMPC2.usererror('Cannot define quadratic cost on algebraic variables unless "include_algebraic_variables_at_stages" is set to true.')
            end
            if size(options.quad_Z,1) ~= C.discretizer.model.dim.algeb
               TRYMPC2.usererror('The quadratic cost vector (quad_Z) must have the same length as the algeb vector.')
            end
            quad_cost.algeb = casadi.SX.zeros(C.discretizer.model.dim.algeb,C.N_horizon);
            for k = 1:C.N_horizon
               stage_deviation = C.decision.str.algeb(:,k) - C.reference.algeb(:,k);
               quad_cost.algeb(:,k) = options.quad_Z.*(stage_deviation.^2);
            end
         end

         % input
         if isfield(options,'quad_R')
            if size(options.quad_R,1) ~= C.discretizer.model.dim.input
               TRYMPC2.usererror('The quadratic cost vector (quad_R) must have the same length as the input vector.')
            end
            quad_cost.input = casadi.SX.zeros(C.discretizer.model.dim.input,C.N_horizon);
            for k = 1:C.N_horizon
               stage_deviation = C.decision.str.input(:,k) - C.reference.input(:,k);
               quad_cost.input(:,k) = options.quad_R.*(stage_deviation.^2);
            end
         end

         % ouput
         if isfield(options,'quad_Y')
            if size(options.quad_Y,1) ~= C.discretizer.model.dim.output
               TRYMPC2.usererror('The quadratic cost vector (quad_Y) must have the same length as the output vector.')
            end
            
            % prepare args:
            args = C.discretizer.model.args;

            quad_cost.output = casadi.SX.zeros(C.discretizer.model.dim.output,C.N_horizon);
            for k = 1:C.N_horizon
               % get output from state vector:
               args.state = C.reference.state(:,k);
               output_k = C.discretizer.model.output.call(args).out;

               stage_deviation = output_k - C.reference.output(:,k);
               quad_cost.output(:,k) = options.quad_Y.*(stage_deviation.^2);
            end
         end

         C.objective.quadratic = quad_cost;
         % %% Objective
         % quad_cost = casadi.SX.zeros;
         % 
         % % state
         % if isfield(options,'quad_Q')
         %    if size(options.quad_Q,1) ~= C.discretizer.model.dim.state
         %       TRYMPC2.usererror('The quadratic cost vector (quad_Q) must have the same length as the state vector.')
         %    end
         %    for k = 1:C.N_horizon
         %       stage_deviation = C.decision.str.state(:,k) - C.reference.state(:,k);
         %       quad_cost = quad_cost + options.quad_Q.*(stage_deviation.^2);
         %    end
         % end
         % 
         % % algeb
         % if isfield(options,'quad_Z')
         %    if ~options.include_algebraic_variables_at_stages
         %       TRYMPC2.usererror('Cannot define quadratic cost on algebraic variables unless "include_algebraic_variables_at_stages" is set to true.')
         %    end
         %    if size(options.quad_Z,1) ~= C.discretizer.model.dim.algeb
         %       TRYMPC2.usererror('The quadratic cost vector (quad_Z) must have the same length as the algeb vector.')
         %    end
         %    for k = 1:C.N_horizon
         %       stage_deviation = C.decision.str.algeb(:,k) - C.reference.algeb(:,k);
         %       quad_cost = quad_cost + options.quad_Z.*(stage_deviation.^2);
         %    end
         % end
         % 
         % % input
         % if isfield(options,'quad_R')
         %    if size(options.quad_R,1) ~= C.discretizer.model.dim.input
         %       TRYMPC2.usererror('The quadratic cost vector (quad_R) must have the same length as the input vector.')
         %    end
         %    for k = 1:C.N_horizon
         %       stage_deviation = C.decision.str.input(:,k) - C.reference.input(:,k);
         %       quad_cost = quad_cost + options.quad_R.*(stage_deviation.^2);
         %    end
         % end
         % 
         % % ouput
         % if isfield(options,'quad_Y')
         %    if size(options.quad_Y,1) ~= C.discretizer.model.dim.output
         %       TRYMPC2.usererror('The quadratic cost vector (quad_Y) must have the same length as the output vector.')
         %    end
         % 
         %    % prepare args:
         %    args = C.discretizer.model.args;
         % 
         %    for k = 1:C.N_horizon
         %       % get output from state vector:
         %       args.state = C.reference.state(:,k);
         %       output_k = C.discretizer.model.output.call(args).out;
         % 
         %       stage_deviation = output_k - C.reference.output(:,k);
         %       quad_cost = quad_cost + options.quad_Y.*(stage_deviation.^2);
         %    end
         % end
         % 
         % C.objective.quadratic = sum(transpose(quad_cost));


         % keep track of fields required:
         C.args.objective.decision = [];
         C.args.objective.reference = [];
         C.args.objective.parameters = [];


         %% Add variable bounds
         C.constraints.inequality = structor("default_mix","TRYMPC_horizon");
         for type = ["state","algeb","input","output"] % OBS! to make the "stage_0" field appear first, the "state" type should be the last type, since this start at stage_1
            bound_type = ['bounds_',char(type)];
            if isfield(options,bound_type)
               if ~isfield(C.constraints.inequality.str,'bounds')
                  C.constraints.inequality.str.bounds.stage_0 = struct; % Make sure that the first stage is stage_0 (this will not happen is "state" has bounds, since the first state bound is k=1)
               end
               if type == "algeb" && ~options.include_algebraic_variables_at_stages
                  TRYMPC2.usererror('To add bounds on the algebraic variable, the option "include_algebraic_variables_at_stages" must be set to true.')
               end
               if isfield(options,bound_type)
                  if size(options.(string(['bounds_',char(type)])),1) ~= C.discretizer.model.dim.(type)
                     TRYMPC2.usererror(['Bounds-vector on ',char(type),'s must be the same size as the ',char(type),'-vector.'])
                  end

                  switch type
                     case {"state","input","algeb"}
                        var_vector = C.decision.str.(type);
                     case "output"
                        % prepare args:
                        args = C.discretizer.model.args;
                        var_vector = casadi.SX.zeros([C.discretizer.model.dim.output,C.N_horizon]);
                        for k = 1:C.N_horizon
                           % get output from state vector:
                           args.state = C.decision.state(:,k);
                           var_vector(:,k)  = C.discretizer.model.output.call(args).out;
                        end

                  end

                  
                  for k = (type == "state"):C.N_horizon-1 % if "state"
                     % add each element bound separately in order to leave out
                     % infinite bounds. (fewer constraints)
                     for i = 1:C.discretizer.model.dim.(type)
                        %%% Lower bounds:
                        bound = options.(bound_type)(i,1);
                        if ~isinf(bound) && ~isnan(bound)% Add bound only if it is finitie and not NaN.
                           C.constraints.inequality.str.bounds.(['stage_',num2str(k)]).(type).lower.(C.discretizer.model.names.(type)(i)) = var_vector(i,k+1) - bound;
                        end
                        %%% Upper bounds:
                        bound = options.(bound_type)(i,2);
                        if ~isinf(bound) && ~isnan(bound)% Add bound only if it is finitie and not NaN.
                           C.constraints.inequality.str.bounds.(['stage_',num2str(k)]).(type).upper.(C.discretizer.model.names.(type)(i)) = bound - var_vector(i,k+1);
                        end
                     end
                  end
               end
            end
         end

         % keep track of fields required:
         C.args.inequality.decision = [];
         C.args.inequality.parameter = [];

         disp(['done.  ',sec2str(toc(def_time)),' Name: "',char(C.Name),'"'])
      end
   end



   % USER METHODS:
   methods
      function C = def_objective(C,type,stage_cost)
         arguments
            C
            type
            stage_cost
         end

         warning('OOPS: not implemented yet...')
      end




      function [numeric_decision,sol,solver] = solve(C,numeric_model,options)
         arguments
            C

            % Basic:
            numeric_model trympcNUMERIC_MODEL
            options.T_horizon (1,1) double {mustBePositive} = 10;
            options.start_time (1,1) double {mustBeReal} = 0;

            % Settigns:
            options.initial_guess % a vector or a structor
            options.max_iterations (1,1) {mustBeInteger,mustBePositive} % takes priority over the nlpsol options
            options.nlpsol_options (1,1) struct = struct('print_time',0);
            
         end


         if ~isfield(options.nlpsol_options,'ipopt')
            options.nlpsol_options.ipopt.print_level = 5;
            options.nlpsol_options.ipopt.max_iter = 1000;

            % Default values:
            % options.nlpsol_options.ipopt.tol = 1e-8;
            options.nlpsol_options.ipopt.acceptable_tol = 1e-6;
            options.nlpsol_options.ipopt.compl_inf_tol = 1e-4;
            options.nlpsol_options.ipopt.constr_viol_tol = 1e-4;
            options.nlpsol_options.ipopt.dual_inf_tol = 1e-4;

            %%% Some other ipopt-options:
            % options.ipopt_nlpsol_options.ipopt.hessian_approximation      = 'limited-memory';
            % options.ipopt_nlpsol_options.ipopt.limited_memory_update_type = 'bfgs';
            % options.ipopt_nlpsol_options.ipopt.linear_solver              = 'mumps';
            % options.ipopt_nlpsol_options.ipopt.linear_system_scaling      = 'none';

         end

         if isfield(options,'max_iterations')
            options.nlpsol_options.ipopt.max_iter = options.max_iterations;
         end

         if ~isfield(options,'initial_guess')
            initial_guess = zeros(C.decision.len,1);
         elseif isa(options.initial_guess,'double')
            initial_guess = options.initial_guess;
            if length(initial_guess) ~= C.decision.len
               TRYMPC2.usererror('The initial_guess has the wrong size.')
            end
         elseif isa(options.initial_guess,'structor')
            initial_guess = options.initial_guess.vec;
         end


         % ===========================
         % Apply the decision vector:
         solver_def.x = C.decision.vec;



         % ===========================
         %%%%%%% Apply the constraints: 
         % prepare casadi-Function (in order to correctly apply parameter values, etc.)
         Args.decision = C.decision.vec;
         Args.initial_condition = C.initial_condition;
         Args.T_horizon = C.T_horizon;
         Args.parameters = C.discretizer.model.param.vec;
         infields = fieldnames(Args);
         Args.out = [C.constraints.equality.vec; C.constraints.inequality.vec];
         g = casadi.Function('ipopt_constraints',Args,infields,'out');

         % generate expression:
         Args = struct;
         Args.decision = C.decision.vec;
         Args.initial_condition = numeric_model.initial_state.vec;
         Args.T_horizon = options.T_horizon;
         Args.parameters = numeric_model.param.vec;
         solver_def.g = g.call(Args).out;

         % Define bounds (0 < eq(x) < 0) and (0 < in(x) < inf)
         lb_eq = zeros(C.constraints.equality.len,1);
         ub_eq = lb_eq;
         lb_in = zeros(C.constraints.inequality.len,1);
         ub_in = inf(C.constraints.inequality.len,1);




         % =======================
         %&&&&& Apply objective:
         % create casadi-Function
         Args = struct;
         Args.decision = C.decision.vec;
         for type = ["state","algeb","input","output"] 
            if isfield(C.reference,type)
               Args.(['reference_',char(type)]) = C.reference.(type);
            end
         end
         Args.parameters = C.discretizer.model.param.vec;
         infields = fieldnames(Args);
         Args.out = trympcDOP.sum_structure(C.objective);
         f = casadi.Function('ipopt_objective',Args,infields,'out');

         % gererate expression
         Args = struct;
         Args.decision = C.decision.vec;
         for type = ["state","algeb","input","output"]
            if isfield(C.reference,type)
               Args.(['reference_',char(type)]) = numeric_model.(['ref_',char(type)])(options.start_time + cumsum(C.rDT).*options.T_horizon);
            end
         end
         Args.parameters = numeric_model.param.vec;
         solver_def.f = f.call(Args).out;




         % ======================
         %%%%%%%%% Define solver:
         solver = casadi.nlpsol('solver', 'ipopt', solver_def, options.nlpsol_options);




         % ===================
         %%%%%%%% CALL IPOPT:
         tic
         sol = solver('x0',initial_guess,... % Initial guess (of decision variables / primal solution)
                      'lbg', [lb_eq; lb_in],...  % Lower bound on inequality vector "g" in casadi-language
                      'ubg', [ub_eq; ub_in]);    % Upper bound on inequality vector "g" in casadi-language
         solve_time = toc;
         numeric_decision = C.decision.retrieve(full(sol.x));

         disp( '================== IPOPT return message ===========')
         disp(['      -- return status: ',solver.stats.return_status])
         disp(['      --    solve time: ',sec2str(solve_time)])
         disp(['      -- N. iterations: ',num2str(solver.stats.iter_count)])
         % disp( '   Convergence Measure: ')
         % disp(['              -- stationarity: ',num2str(stationarity)]) % C.archive.optimizations{end}.convergence.stationarity(end)
         % disp(['              --     equality: ',num2str(equality)]) % C.archive.optimizations{end}.convergence.equality(end)
         % disp(['              --   inequality: ',num2str(inequality)]) % C.archive.optimizations{end}.convergence.inequality(end)
         disp( '====================================================')
      end



      function Fig = show(C,decision,types)
         arguments
            C
            decision structor % containing the numeric values at the to display
            types (1,:) string {mustBeMember(types,["state","algeb","input","output"])} = ["state","input"];
         end


         for type = types
            Fig.(type) = figure("Name",C.ID+"-"+C.Name+" ("+type+")");
            Tiles = tiledlayout('flow');
            title(Tiles,type)
            for i = 1:C.discretizer.model.dim.(type)
               ax = nexttile(i);
               ylabel(ax,C.discretizer.model.names.(type)(i))
               hold on
               plot(ax,decision.str.(type)(i,:))
            end
         end

      end
   end





   % Internal methods
   methods(Static)
      function out = sum_structure(S)
         out = [];
         for name = string(fieldnames(S))'
            if isa(S.(name),'struct')
               value = trympcDOP.sum_structure(S.(name));
            else
               value = sum(S.(name)(:));
            end
            % Instead of just settign out = 0 at the start, we set it to
            % empty and define the value at the first value.
            % Doing it this way allows compatibility with f.ex. casadi
            % symbols (which don't accept mixing with the 'double' 0).
            if isempty(out)
               out = value;
            else
               out = out + value;
            end
         end
      end

   end


end