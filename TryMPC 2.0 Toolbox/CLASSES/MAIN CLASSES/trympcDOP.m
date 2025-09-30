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
      % args struct = struct; % a struct that contains the fields required (by default) to evaluate the various casadi-functions (objective, constraints, ...)
   end

   % User attainable symbolic variables:
   properties(SetAccess=private)
      
      % Input variables:
      initial_state casadi.SX % a variable that represent the initial condition (a parameter)
      reference struct = struct;
      T_horizon casadi.SX % the total length of the horizon (in time)
      quad_cost struct = struct; % struct that may contain fields "state", "input", "output", each of which is a vector of weights (symbolic)
      decision structor = structor("default_mix","TRYMPC_horizon");

   end
   properties
      % Expressions:
      constraints struct = struct;
      objective struct = struct;
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


            % %%% Possibility of quadratic stage cost
            options.quad_cost (1,:) string {mustBeMember(options.quad_cost,["state","input","output"])}

            %%% add variable bounds
            options.bounds_state (:,2) double {mustBeReal} % [lower upper] bound on state vector
            options.bounds_input (:,2) double {mustBeReal} % [lower upper] bound on input vector
            options.bounds_output(:,2) double {mustBeReal} % [lower upper] bound on output vector

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
         C.decision.str.state = casadi.SX.sym('decision_state',[C.discretizer.model.dim.state,C.N_horizon+1]); % contains state_0 -> state_(N)
         C.decision.str.input = casadi.SX.sym('decision_input',[C.discretizer.model.dim.input,C.N_horizon]); % contains input_0 -> input_(N-1)
         
         % Reference Variables:
         C.reference.state = casadi.SX.sym('reference_state',[C.discretizer.model.dim.state,C.N_horizon]); % contains state_1 -> state_(N)
         C.reference.input = casadi.SX.sym('reference_input',[C.discretizer.model.dim.input,C.N_horizon]); % contains input_0 -> input_(N-1)
         if ismember(C.discretizer.model.expression_types,"output")
            C.reference.output = casadi.SX.sym('reference_output',[C.discretizer.model.dim.output,C.N_horizon+1]); % contains state_1 -> state_(N)
         end

         %{
We emphazise that the all algebraic variables that are necessary for
propogating the system are present in the "aux" variables already, and are
used for evaluating the dynamics inside the discretizer. Now new algebraic
variables are necessary at the stages of the DOP. Though if you want to
add constraints of costs associated with them, it might be nice to predict
them too. There this option is added.
         %}


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
 ensure:        state_(k+1) = discretizer.Next(state_k,input_k,param,aux_var,Dt)
 and:                     0 = discretizer.Algebraics(state_k,input_k,param,aux_var,Dt)
at each stage k \in {0,...,N-1}
         %}

         % Prepare structor and arguments
         args = struct;
         C.constraints.dynamic = structor("default_mix","vector_mixed");

         % initial condition constraint:
         C.initial_state = casadi.SX.sym('initial_state',[C.discretizer.model.dim.state,1]);
         C.constraints.dynamic.str.initial_condition = C.decision.str.state(:,1) - C.initial_state;

         for k = 0:(C.N_horizon-1)
            args.state = C.decision.str.state(:,k+1);
            args.input = C.decision.str.input(:,k+1);
            args.aux_var = aux_matrix(:,k+1);
            args.param = C.discretizer.model.param.vec;
            args.Dt = DT(k+1);
            if ~isempty(C.discretizer.Algebraics)
               C.constraints.dynamic.str.(['step_',num2str(k+1)]).integrator = C.discretizer.Algebraics.call(args).out;
            end
            C.constraints.dynamic.str.(['step_',num2str(k+1)]).continuity   = C.decision.str.state(:,k+2) - C.discretizer.Next.call(args).out;
         end
         


         %% Objective

         % prepare args (for output computation)
         args = C.discretizer.model.args;

         % create quadratic cost terms on horizon
         if isfield(options,'quad_cost')
            for type = options.quad_cost
               C.quad_cost.(type) = casadi.SX.sym(['quad_',char(type)],[C.discretizer.model.dim.(type),1]);
               C.objective.quadratic.(type) = casadi.SX.zeros(C.discretizer.model.dim.(type),C.N_horizon);


               for k = (1:C.N_horizon)
                  if type == "output"
                     % get output from state vector:
                     args.state = C.decision.state(:,k+1);
                     output_k = C.discretizer.model.output.call(args).out;
                     stage_deviation = output_k - C.reference.output(:,k); % note that the reference vector is offset, since the state_0 variable does not have a reference
                  elseif type == "state"
                     stage_deviation = C.decision.str.(type)(:,k+1) - C.reference.(type)(:,k); % note that the reference vector is offset, since the state_0 variable does not have a reference
                  else
                     stage_deviation = C.decision.str.(type)(:,k) - C.reference.(type)(:,k);
                  end
                  C.objective.quadratic.(type)(:,k) = C.quad_cost.(type).*(stage_deviation.^2);
               end

            end
         end




         %% Add variable bounds
         C.constraints.inequality = structor("default_mix","TRYMPC_horizon");
         for type = ["state","input","output"]
            bound_type = ['bounds_',char(type)];
            if isfield(options,bound_type)
               if ~isfield(C.constraints.inequality.str,'bounds')
                  C.constraints.inequality.str.bounds.stage_0 = struct; % Make sure that the first stage is stage_0 (this will not happen is "state" has bounds, since the first state bound is k=1)
               end
               if size(options.(string(['bounds_',char(type)])),1) ~= C.discretizer.model.dim.(type)
                  TRYMPC2.usererror(['Bounds-vector on ',char(type),'s must be the same size as the ',char(type),'-vector.'])
               end

               % prepare variable vector:
               switch type
                  case {"state","input"}
                     var_vector = C.decision.str.(type);
                  case "output"
                     % prepare args:
                     args = C.discretizer.model.args;
                     var_vector = casadi.SX.zeros([C.discretizer.model.dim.output,C.N_horizon]);
                     for k = 1:C.N_horizon
                        % get output from state vector:
                        args.state = C.decision.str.state(:,k);
                        var_vector(:,k)  = C.discretizer.model.output.call(args).out;
                     end
               end


               % add bounds:
               for k = (0:C.N_horizon-1)+ismember(type,["state","output"])
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

         disp(['done.  ',sec2str(toc(def_time)),' Name: "',char(C.Name),'"'])
      end
   end








   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%% USER METHODS:
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods
      function [sol,solver] = solve(C,numeric_model,options)
         arguments
            C

            % Basic:
            numeric_model trympcNUMERIC_MODEL
            options.T_horizon (1,1) double {mustBePositive} = 10;
            options.start_time (1,1) double {mustBeReal} = 0;
            options.F (1,1) struct = C.make_Functions % a struct that contains the casadi-Functions (hence "F") with the problem formulation
            options.quad (1,1) struct % struct containing the quadratic weights if relevant
            %{
Fields of F:
-    "dynamic":   contains the dynamic constraints ( dyn(x) = 0 )
-  "objective":   contains the objective function (to minimize)
- "inequality":   contains the inequality constraints ( h(x) >= 0 )
(equality constraints should be included in the inequalities)

Fields of quad:
-  "state": weight-vector (column) for quadratic states
-  "input": weight-vector (column) for quadratic inputs
- "output": weight-vector (column) for quadratic outputs
            %}

            % Settigns:
            options.initial_guess % a vector or a structor
            options.max_iterations (1,1) {mustBeInteger,mustBePositive} % takes priority over the nlpsol options
            options.nlpsol_options (1,1) struct
            options.print_level (1,1) double {mustBeInteger,mustBeInRange(options.print_level,0,12)} % takes priority over the nlpsol options

         end

         fprintf('Preparing ipopt-solver...  ');
         def_time = tic;


         % =====================================
         %%%%%%%%%% ENSURE PROPER IPOPT OPTIONS:

         % % Default values:
         % % options.nlpsol_options.ipopt.tol = 1e-8;
         % options.nlpsol_options.ipopt.acceptable_tol = 1e-6;
         % options.nlpsol_options.ipopt.compl_inf_tol = 1e-4;
         % options.nlpsol_options.ipopt.constr_viol_tol = 1e-4;
         % options.nlpsol_options.ipopt.dual_inf_tol = 1e-4;

         %%% Some other ipopt-options:
         % options.ipopt_nlpsol_options.ipopt.hessian_approximation      = 'limited-memory';
         % options.ipopt_nlpsol_options.ipopt.limited_memory_update_type = 'bfgs';
         % options.ipopt_nlpsol_options.ipopt.linear_solver              = 'mumps';
         % options.ipopt_nlpsol_options.ipopt.linear_system_scaling      = 'none';

         %%%%%%% IPOPT Suggestion:
         nlpsol_options = struct;
         nlpsol_options.print_time = 0;

         % IPOPT main settings
         ipopt.tol = 1e-6;                % Overall convergence tolerance
         ipopt.dual_inf_tol = 1e-6;       % Dual infeasibility tolerance
         ipopt.constr_viol_tol = 1e-6;    % Constraint violation tolerance
         ipopt.compl_inf_tol = 1e-6;      % Complementarity tolerance

         % Linear solver
         ipopt.linear_solver = 'mumps';   % Robust general-purpose solver
         % Alternatives: 'ma57', 'ma86', 'ma97' (HSL, faster if licensed)

         % Regularization and scaling
         ipopt.hessian_approximation = 'exact';   % Full Hessian by CasADi
         % Alternatives: 'limited-memory' (for large-scale problems)
         % ipopt.limited_memory_update_type = 'bfgs'; % I added this option
         ipopt.fixed_variable_treatment = 'make_constraint';
         ipopt.linear_system_scaling = 'mc19';    % Improves conditioning

         % Iteration control
         ipopt.max_iter = 1000;           % Maximum iterations
         ipopt.print_level = 5;           % IPOPT output verbosity (0–12)
         ipopt.print_timing_statistics = 'no';

         % Initialization and step length
         ipopt.acceptable_tol = 1e-4;     % Early stopping tolerance
         ipopt.acceptable_iter = 5;       % How many acceptable steps before quit
         ipopt.bound_relax_factor = 0.0;  % Do not relax bounds artificially
         ipopt.start_with_resto = 'no';   % Don’t start with restoration phase

         % Initialization & regularization
         ipopt.mu_strategy           = 'adaptive'; % Barrier update strategy
         ipopt.mu_init               = 1e-1;      % Initial barrier parameter
         ipopt.warm_start_init_point = 'yes';     % Warm start if re-solving

         if isfield(options,'nlpsol_options')
            for name = string(fieldnames(options.nlpsol_options))'
               nlpsol_options.(name) = options.nlpsol_options.(name);
            end
            nlpsol_options.ipopt = ipopt;
            if isfield(options.nlpsol_options,'ipopt')
               for name = string(fieldnames(options.nlpsol_options.ipopt))'
                  nlpsol_options.ipopt.(name) = options.nlpsol_options.ipopt.(name);
               end
            end
         end
         for name = ["max_iterations", "print_level"]
            if isfield(options,name)
               nlpsol_options.ipopt.(name) = options.(name);
            end
         end


         % =====================
         %%%%%%%% INITIAL GUESS:
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



         % ===================================================
         % Compute sampling times (time of the various stages)
         sampling_times = options.start_time + cumsum(C.rDT).*options.T_horizon;



         % ===========================
         % Apply the decision vector:
         solver_def.x = C.decision.vec;



         % ===========================
         %%%%%%% Create problem (call functions where all but decision variables are numeric)
         
         % inequality constrainst:
         ineq_args = struct;
         ineq_args.decision = C.decision.vec;
         ineq_args.parameters = numeric_model.parameters.vec;
         if C.constraints.inequality.len > 0
         inequalities = options.F.inequality.call(ineq_args).out;
         ineq_ones = ones(length(inequalities),1);
         else
            inequalities = [];
            ineq_ones = [];
         end

         % Dynamic constraints:
         dyn_args = ineq_args;
         dyn_args.T_horizon = options.T_horizon;
         dyn_args.initial_state = numeric_model.initial_state.vec;
         dynamic_constraints = options.F.dynamic.call(dyn_args).out;
         dyn_zero = zeros(length(dynamic_constraints),1);

         % Objective:
         obj_args = struct;
         obj_args.decision = C.decision.vec;
         obj_args.parameters = numeric_model.parameters.vec;
         for type = string(fieldnames(C.quad_cost))'
            if ~isfield(numeric_model.ref,type)
               TRYMPC.usererror(['The numeric model has not defined a reference for "',char(type),'"'])
            end
            if ~isfield(options.quad,type)
               TRYMPC.usererror(['The quadratic cost vector "quad.',char(type),'" is not provided as argument for "',char(type),'"'])
            end
            obj_args.(['reference_',char(type)]) = numeric_model.ref.(type)(sampling_times);
            obj_args.(['quad_',char(type)]) = reshape(options.quad.(type),[],1);
         end
         obj = options.F.objective.call(obj_args).out;


         solver_def.g = [dynamic_constraints; inequalities];
         solver_def.f = obj;

         lb = [dyn_zero;ineq_ones*0];
         ub = [dyn_zero;ineq_ones*inf];




         % ======================
         %%%%%%%%% Define solver:
         solver = casadi.nlpsol('solver', 'ipopt', solver_def, nlpsol_options);


         

         % ===================
         %%%%%%%% CALL IPOPT:
         disp(['done.  ',sec2str(toc(def_time)),' Name: "',char(C.Name),'"'])
         fprintf('solving...  ');
         solve_time = tic;
         sol_ipopt = solver('x0',initial_guess,... % Initial guess (of decision variables / primal solution)
                      'lbg', lb,...  % Lower bound on inequality vector "g" in casadi-language
                      'ubg', ub);    % Upper bound on inequality vector "g" in casadi-language
         solve_time = toc(solve_time);
         disp(['done.  ',sec2str(solve_time)])
         numeric_decision = full(sol_ipopt.x);




         % ===========================================
         %%%%%%%%% Retrieve solutions and constraints:
         sol.ipopt = sol_ipopt;
         sol.decision = C.decision.retrieve(numeric_decision);

         dyn_args.decision = numeric_decision;
         ineq_args.decision = numeric_decision;
         obj_args.decision = numeric_decision;

         dynamic = full(options.F.dynamic.call(dyn_args).out);
         sol.constraints.dynamic = C.constraints.dynamic.retrieve(dynamic);
         
         ineq = full(options.F.inequality.call(ineq_args).out);
         sol.constraints.inequality = C.constraints.inequality.retrieve(ineq);


         % ========================
         %%%%%%%%%% Display Result:
         disp( '================== IPOPT return message ===========')
         disp(['      -- return status: ',solver.stats.return_status])
         disp(['      --    solve time: ',sec2str(solve_time)])
         disp(['      -- N. iterations: ',num2str(solver.stats.iter_count)])
         disp( '   Constraint Satisfaction: ')
         disp(['              --      dynamic: ',num2str(norm(dynamic,1)),'  ( ||dyn(x)||_1 )']) % C.archive.optimizations{end}.convergence.equality(end)
         disp(['              --   inequality: ',num2str(min(ineq)),'  ( min(ineq(x)) )']) % C.archive.optimizations{end}.convergence.inequality(end)
         disp( '====================================================')
      end




      % +++++++++++++++++++++++++++
      %%%%%%%%%% DISPLAY SOLUTION:
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



      % ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      %%%%%%%%%% GET ARGUMENTS STRUCT FOR CALLING CASADI-FUNCTIONS:

      % OBJECTIVE
      function args = args_objective(C)
         args.decision = C.decision.vec;
         args.parameters = C.discretizer.model.param.vec;
         for type = string(fieldnames(C.quad_cost))'
            args.(['reference_',char(type)]) = C.reference.(type);
            args.(['quad_',char(type)]) = C.quad_cost.(type);
         end
      end

      % INEQUALITY CONSTRAINTS
      function args = args_inequality(C)
         args.decision = C.decision.vec;
         args.parameters = C.discretizer.model.param.vec;
      end

      % DYNAMIC CONSTRAINTS
      function args = args_dynamic(C)
         args.decision = C.decision.vec;
         args.parameters = C.discretizer.model.param.vec;
         args.T_horizon = C.T_horizon;
         args.initial_state = C.initial_state;
      end




      % ++++++++++++++++++++++++++++++++++++++++++++++++
      %%%%%%%%%% CREATE CASADI-FUNCTIONS OF THE PROBLEM:
      function F = make_Functions(C)
         % Dynamic constraints:
         args = C.args_dynamic;
         infields = fieldnames(args);
         args.out = C.constraints.dynamic.vec;
         F.dynamic = casadi.Function('F_dynamic_constraints',args,infields,'out');

         % Objective:
         args = C.args_objective;
         infields = fieldnames(args);
         args.out = trympcDOP.sum_structure(C.objective.quadratic);
         F.objective = casadi.Function('F_objective',args,infields,'out');

         % Inequality constraints:
         args = C.args_inequality;
         infields = fieldnames(args);
         args.out = C.constraints.inequality.vec;
         F.inequality = casadi.Function('F_inequality_constraints',args,infields,'out');
      end



   end 






   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%% Static Methods:
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   methods(Static)

      % ++++++++++++++++++++++++++++++++++++
      %%%%%%%%%% SUM ALL FIELDS OF A STRUCT:
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