classdef NMPC < trympcCONTROLLER
   
   properties(SetAccess = private)
      primal (:,1) % the primal solution (decision)
      dual (:,1) % the dual solution (lagrange multipliers)

      %%% Define problem
      dop (1,1) % the dynamic optimization problem that is solved each iteration of the NMPC simulation
      casadi_solver (1,1)
   end

   properties
      %%% Problem Settings:
      numeric_model trympcNUMERIC_MODEL % contains references and parameters
      T_horizon (1,1) double
      nlpsol_options (1,1) struct
      quad (1,1) struct % struct containing the quadratic weights if relevant

      toggle_archive (1,1) logical = false;
   end

   properties(SetAccess=private)
      %%% Archive to store info of each solve:
      archive (1,1)
   end

   properties(Hidden,SetAccess=private)
      nlp_parameters % struct containing parameters for the nlpsol object (from casadi)
      % args % struct of arguments for the casadi-Functions
      % solver_def % struct used to define the casadi-solver
      lb % a stacked vector of zeros
      ub % a stacked vector of zeros for equalities and inf for inequalities
   end


   methods
      function C = NMPC(Name,dop,set_prop)
         arguments
            Name (1,1) string
            dop (1,1) trympcDOP
            set_prop.numeric_model (1,1) trympcNUMERIC_MODEL
            set_prop.T_horizon (1,1) double {mustBePositive}
            set_prop.quad (1,1) struct
            set_prop.toggle_archive (1,1) logical
         end

         % ==============
         % INSTANCIATE:
         C@trympcCONTROLLER(Name)
         C.dop = dop;

         for name = string(fieldnames(set_prop))'
            C.(name) = set_prop.(name);
         end


         % ===========================
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
         ipopt.print_level = 0;           % IPOPT output verbosity (0–12)
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

         nlpsol_options.ipopt = ipopt;
         C.nlpsol_options = nlpsol_options;
         % ===============================================





      % =================================================================
      % =================================================================
      % INSTANCIATE THE PROBLEM:
      % ------------------------



         % ========================
         % Create casadi-Functions:
         F = C.dop.make_Functions; % a struct that contains the casadi-Functions (hence "F") with the problem formulation




         % ===========================
         % Apply the decision vector:
         solver_def.x = C.dop.decision.vec;




         % ===========================
         %%%%%%% Create problem (call functions where all but decision variables are numeric)

         %%%% Inequality constrainst:
         arg.ineq = struct;
         arg.ineq.decision = C.dop.decision.vec;
         arg.ineq.parameters = C.dop.discretizer.model.param.vec;
         if C.dop.constraints.inequality.len > 0
            inequalities = F.inequality.call(arg.ineq).out;
         else
            inequalities = [];
         end

         %%%% Dynamic constraints:
         arg.dyn = arg.ineq;
         arg.dyn.T_horizon = C.dop.T_horizon;
         arg.dyn.initial_state = C.dop.initial_state.vec;
         dynamic_constraints = F.dynamic.call(arg.dyn).out;

         % add to nlp-parameter list
         C.nlp_parameters.parameters = arg.dyn.parameters;
         C.nlp_parameters.initial_state = arg.dyn.initial_state;
         C.nlp_parameters.T_horizon = arg.dyn.T_horizon;

         %%%% Objective:
         arg.obj = struct;
         arg.obj.decision = C.dop.decision.vec;
         arg.obj.parameters = C.dop.discretizer.model.param.vec;
         for type = string(fieldnames(C.dop.quad_cost))'
            arg.obj.(['reference_',char(type)]) = C.dop.reference.(type);
            arg.obj.(['quad_',char(type)]) = C.dop.quad_cost.(type);
            C.nlp_parameters.(['reference_',char(type)]) = arg.obj.(['reference_',char(type)]);
            C.nlp_parameters.(['quad_',char(type)]) = arg.obj.(['quad_',char(type)]);
         end
         obj = F.objective.call(arg.obj).out;


         solver_def.g = [dynamic_constraints; inequalities];
         solver_def.f = obj;
         solver_def.p = C.stack_parameters;



         dyn_zeros = zeros(prod(F.dynamic.size_out(0)),1);
         ineq_ones = ones(prod(F.inequality.size_out(0)),1);
         C.lb = [dyn_zeros;ineq_ones*0];
         C.ub = [dyn_zeros;ineq_ones*inf];



         % ======================
         %%%%%%%%% Define solver:
         C.casadi_solver = casadi.nlpsol('solver', 'ipopt', solver_def, C.nlpsol_options);




         % ==========================
         %%%%%%% Set an initial guess:
         C.primal = zeros(C.dop.decision.len,1);

      end

   end
   methods(Access=protected)
      function u = internal_control(C,time,state)
         
         % ===================================================
         % Compute sampling times (time of the various stages)
         sampling_times =  time + cumsum(C.dop.rDT).*C.T_horizon;

         % add to nlp-parameter list
         C.nlp_parameters.initial_state = state;
         C.nlp_parameters.parameters = C.numeric_model.parameters.vec;
         C.nlp_parameters.T_horizon = C.T_horizon;
         for type = string(fieldnames(C.dop.quad_cost))'
            C.nlp_parameters.(['reference_',char(type)]) = C.numeric_model.(['reference_',char(type)])(sampling_times);
            C.nlp_parameters.(['quad_',char(type)]) = C.quad.(type);
         end


         % ===================
         %%%%%%%% CALL IPOPT:
         sol_ipopt = C.casadi_solver('x0', C.primal,... % Initial guess (of decision variables / primal solution)
                      'lbg', C.lb,...  % Lower bound on inequality vector "g" in casadi-language
                      'ubg', C.ub,...  % Upper bound on inequality vector "g" in casadi-language
                      'p', C.stack_parameters); % provide model-parameters and references as parameters
         C.primal = full(sol_ipopt.x);


         % ===========================================
         %%%%%%%%% Retrieve solutions and constraints:
         arch_ind = length(C.archive)+1;
         C.archive(arch_ind).decision     = C.primal;
         C.archive(arch_ind).iter_count   = C.casadi_solver.stats.iter_count;
         C.archive(arch_ind).success      = C.casadi_solver.stats.success;
         C.archive(arch_ind).dual         = sol_ipopt.lam_g;
         C.archive(arch_ind).time         = time;
         C.archive(arch_ind).state        = state;




         % ==============================
         %%%%%% EXTRACT CONTROL VARIABLE:
         u_traj = C.primal(C.dop.decision.dim.input);
         u = @(t) interp1(sampling_times',u_traj',t)';

      end



      %%%%%% stack parameters:
      function p_stack = stack_parameters(C)
         fields = fieldnames(C.nlp_parameters);
         param_cells = cellfun(@(f) reshape(C.nlp_parameters.(f),numel(C.nlp_parameters.(f)),1), fields, 'UniformOutput', false);
         p_stack = vertcat(param_cells{:});
      end
   end
end