classdef trympcDISCRETIZER
%{
This class defines a single integration step.
That is, it generates the discrete dynamcis x_next = f_d(x_now), based on
the integration scheme of your choice.
%}
   properties(SetAccess=immutable)
      Name (1,1) string 
      order (1,1) double {mustBeNonnegative} = 0; % set to 0 for custom methods with unknown order...
      plicity (1,1) string {mustBeMember(plicity,["explicit","implicit"])} = "explicit";
      integrator (1,1) string % name of integration scheme (if given a name)
      n_increments (1,1) double {mustBeInteger,mustBePositive} = 1;

      % CasADi functions
      Next % gives the next state value (discrete dynamics)
      Algebraics % these expressions must be zero for "Next" be the correct dynamics
   end
   
   properties(SetAccess=immutable,Hidden)
      ID (1,1) string
      model (1,1) 
      BT_b double
      BT_A double
      aux_var (1,1) structor
      discretizer_args (1,1) struct
      collocation (1,1) struct
   end

   methods
      function C = trympcDISCRETIZER(integrator_name,model,method,options)
         arguments
            integrator_name
            model (1,1) trympcMODEL
            method (1,1) string {mustBeMember(method,["Explicit Euler","Implicit Euler","ERK4","ERK4 (simultaneous)","Implicit Midpoint","Crank-Nicolson (Implicit)","IRK4 (L-stable)","Gauss-Legendre (4. order)","Gauss-Legendre (6. order)","custom collocation","custom explicit butcher tableau","custom implicit butcher tableau"])}

            % General settings:
            options.n_increments (1,1) double {mustBePositive,mustBeInteger} = 1;

            % Collocaiton Specific settigns:
            options.collocation_polynomial_order (1,1) double {mustBePositive,mustBeInteger} = 2;
            options.collocation_polynomial_type (1,1) string {mustBeMember(options.collocation_polynomial_type,["legendre","radau"])} = "legendre";

            options.bucher_tableau_b (1,:) double {mustBeReal}
            options.bucher_tableau_A (:,:) double {mustBeReal}
         end

         fprintf('Creating Discretizer...  ');
         def_time = tic;
         
         % % Call superclass constructor first
         % C@TRYMPC2_overclass(integrator_name);
         C.Name = integrator_name;
         C.ID = TRYMPC2.generate_id;

         % reset, since this lingers from the previous instantiation of a
         % trympcDISCRETIZER object...
         C.aux_var = structor("default_mix","TRYMPC_horizon");

         % Set model:
         C.model = model;

         % Store information:
         C.integrator = method;
         C.n_increments = options.n_increments;

         % Define time-step variables:
         Dt = casadi.SX.sym('Dt');
         dt = (Dt/C.n_increments);

         % Misc:
         aux_expr = structor("default_mix","TRYMPC_horizon");
         next_state = C.model.state.vec;

         %% Select Integrator

         % define integrator:
         switch method
            case "Explicit Euler"
               C.order = 1;
               C.plicity = "explicit";
               
               % Butcher Tableau:
               C.BT_b = 1;
               C.BT_A = 0;

               ERK_builder

            case "Implicit Euler"

               C.order = 1;
               C.plicity = "implicit";

               % Butcher Tableau:
               C.BT_b = 1;
               C.BT_A = 1;

               IRK_builder

            case "ERK4"
               C.order = 4;
               C.plicity = "explicit";

               % Butcher Tableau:
               C.BT_b = [1 2 2 1]/6;
               C.BT_A = zeros(4);
               C.BT_A(2,1) = 1/2;
               C.BT_A(3,2) = 1/2;
               C.BT_A(4,3) = 1;

               ERK_builder

            case "ERK4 (simultaneous)"
               C.order = 4;
               C.plicity = "implicit";

               % Butcher Tableau:
               C.BT_b = [1 2 2 1]/6;
               C.BT_A = zeros(4);
               C.BT_A(2,1) = 1/2;
               C.BT_A(3,2) = 1/2;
               C.BT_A(4,3) = 1;

               IRK_builder

            case "Implicit Midpoint"

               C.order = 2;
               C.plicity = "implicit";

               % Butcher Tableau:
               C.BT_b = 1;
               C.BT_A = 1/2;

               IRK_builder

            case "Crank-Nicolson (Implicit)"

               C.order = 2;
               C.plicity = "implicit";

               % Butcher Tableau:
               C.BT_b = [1 1]/2;
               C.BT_A = [0 0; 1 1]/2;

               IRK_builder

            case "Gauss-Legendre (4. order)"

               C.order = 4;
               C.plicity = "implicit";

               % Butcher Tableau:
               C.BT_b = [1 1]/2;
               C.BT_A = [   1/4         1/4-sqrt(3)/6 ;
                  1/4+sqrt(3)/6      1/4       ];

               IRK_builder

            case "Gauss-Legendre (6. order)"

               C.order = 6;
               C.plicity = "implicit";

               % Butcher Tableau:
               C.BT_b = [5/18 4/9 5/18];
               C.BT_A = [5/36              2/9-sqrt(15)/15   5/36-sqrt(15)/30;
                                    5/36+sqrt(15)/24      2/9           5/36-sqrt(15)/24;
                                    5/36+sqrt(15)/30  2/9+sqrt(15)/15   5/36];

               IRK_builder

            case "IRK4 (L-stable)"

               C.order = 3;
               C.plicity = "implicit";

               % Butcher Tableau:
               C.BT_b = [3 -3 1 1]/2;
               C.BT_A = [ 1   0   0   0  ;
                                    1/3  1   0   0  ;
                                    -1   1   1   0  ;
                                     3  -3   1   1  ] /2;

               IRK_builder

            case "custom explicit butcher tableau"
               if ~isfield(options,'bucher_tableau_b')
                  TRYMPC2.usererror('the b-vector of a butcher tableau must be defined to use a custom butcher tableau. Try f.ex. "ERK4" for a predefined RK method.')
               end
               if ~isfield(options,'bucher_tableau_A')
                  TRYMPC2.usererror('the A-matrix of a butcher tableau must be defined to use a custom butcher tableau. Try f.ex. "ERK4" for a predefined RK method.')
               end

               C.order = 0;
               C.plicity = "explicit";

               % Butcher Tableau:
               C.BT_b = options.bucher_tableau_b;
               C.BT_A = options.bucher_tableau_A;

               ERK_builder

            case "custom implicit butcher tableau"
               if ~isfield(options,'bucher_tableau_b')
                  TRYMPC2.usererror('the b-vector of a butcher tableau must be defined to use a custom butcher tableau. Try f.ex. "IRK4" for a predefined RK method.')
               end
               if ~isfield(options,'bucher_tableau_A')
                  TRYMPC2.usererror('the A-matrix of a butcher tableau must be defined to use a custom butcher tableau. Try f.ex. "IRK4" for a predefined RK method.')
               end

               C.order = 0;
               C.plicity = "implicit";

               % Butcher Tableau:
               C.BT_b = options.bucher_tableau_b;
               C.BT_A = options.bucher_tableau_A;

               IRK_builder

            case "custom collocation"

               % save collocation settings:
               C.collocation.d = options.collocation_polynomial_order;
               C.collocation.polynomial_type = options.collocation_polynomial_type;
               
               C.plicity = "implicit";
               
               % specify order
               switch C.collocation.polynomial_type
                  case "legendre"
                     C.order = C.collocation.d*2;
                  case "radau"
                     C.order = C.collocation.d*2 - 1;
                     warning('DEVELOPER ERROR: ops, I am unsure if I have built the collocation scheme specifically for Legendre, of if simply choosing Radau point will procude the correct Radau collocation shceme. - Trym')
                  otherwise
                     error('DEVELOPER ERROR: an invalid colloocation polynomial type was selected')
               end


               %%%%%%%% MAKE LGARANGE POLYNMIAL COEFFIECIENTS
               d = C.collocation.d; % order of integration and order of polynomial
               C.collocation.tau = [0 casadi.collocation_points(d, char(C.collocation.polynomial_type))]; % choose collocaiton points (based on either gauss-legendre or gauss-radau quadrature)
               C.collocation.Li = [];
               C.collocation.dLi = [];
               % Construct Lagrange interpolation polynomials:
               for j=1:d+1
                  coeff = 1;
                  for r=1:d+1
                     if r ~= j
                        coeff = conv(coeff, [1, -C.collocation.tau(r)]);
                        coeff = coeff / (C.collocation.tau(j)-C.collocation.tau(r));
                     end
                  end
                  C.collocation.Li(j,:) = coeff;
                  C.collocation.dLi(j,:) = polyder(coeff);
               end
               %%%%%%%% ... END OF MAKING LAGRANGE POLYNOMIAL COEFFICIENTS

               collocation_builder
         end
         
         %% Create Integrator Functions



         %%% Prepare args:
         C.discretizer_args.state = C.model.args.state;
         C.discretizer_args.input = C.model.args.input;
         C.discretizer_args.param = C.model.args.param;
         C.discretizer_args.aux_var = C.aux_var.vec;
         C.discretizer_args.Dt = Dt;
         arg = C.discretizer_args;
         infields = fieldnames(arg);

         % Discrete dynamics
         arg.out = next_state;
         C.Next = casadi.Function('F_next_state',arg,infields,{'out'}); % discrete dynamics

         % Auxiliary functions (=0)
         if C.plicity == "implicit"
            arg.out = aux_expr.vec;
            C.Algebraics = casadi.Function('F_algebraics',arg,infields,{'out'}); % Expressions that should be zero. (both algebraic equations, and auxiliary integration equations)
         end


         disp(['done.  ',sec2str(toc(def_time)),' Name: "',char(C.Name),'"'])
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END OF CONSTRUCTOR



         %% RK Builders

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% RUNGE-KUTTA:
         function ERK_builder


            % Buthcer Tableau:
            b = reshape(C.BT_b,[],1);
            A = C.BT_A;
            stages = length(b);

            % Prepare variables to simplify creation of integrator:
            in_vars = C.model.args;

            % Integrate system using RK method described by Butcher Tableau:
            for inc = 1:C.n_increments
               f = C.model.dynamics.call(in_vars).out;
               for s = 2:stages
                  in_vars.state = next_state;
                  for a = 1:s-1
                     in_vars.state = in_vars.state + A(s,a)*f(:,a)*dt;
                  end
                  f = [f C.model.dynamics.call(in_vars).out]; %#ok<AGROW>
               end
               next_state = next_state + (f*b).*dt;
            end
         end


         function IRK_builder

            % Buthcer Tableau:
            b = reshape(C.BT_b,[],1);
            A = C.BT_A;
            samp = length(b); % the number of sample points of the RK scheme (a.k.a. "stages", but that can be confused with stages as in each discrete time point on the prediction horizon)

            % First loop over and create relevant variables:
            for inc = 1:C.n_increments
               for s = 1:samp
                  C.aux_var.str.(['inc_',num2str(inc)]).state.(['samp_',num2str(s)]) = casadi.SX.sym(['dynamics_inc',num2str(inc),'_samp',num2str(s)],[C.model.dim.state,1]);
               end
            end

            % Prepare variables to simplify creation of integrator:
            in_vars = C.model.args;

            % Integrate system using RK method described by Butcher Tableau:
            for inc = 1:C.n_increments
               for s = 1:samp
                  in_vars.state = next_state;
                  for a = 1:samp
                     in_vars.state = in_vars.state + A(s,a)*C.aux_var.str.(['inc_',num2str(inc)]).state.(['samp_',num2str(a)])*dt;
                  end

                  aux_expr.str.(['inc_',num2str(inc)]).(['samp_',num2str(s)]).dynamics = C.aux_var.str.(['inc_',num2str(inc)]).state.(['samp_',num2str(s)]) - C.model.dynamics.call(in_vars).out;

               end

               next_state = next_state + (reshape(structor.subvec(C.aux_var,C.aux_var.str.(['inc_',num2str(inc)]).state),C.model.dim.state,samp)*b) .* dt;
            end

         end




         %% COLLOCATION
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%% COLLOCATION
         function collocation_builder


            % Create symbolic variables for collocation coefficients (analogous to state variables at the collocation points)
            for inc = 1:C.n_increments
               for c = 1:C.collocation.d
                  C.aux_var.str.(['inc_',num2str(inc)]).state.(['samp_',num2str(c)]) = casadi.SX.sym(['collocation_state_inc',num2str(inc),'_samp',num2str(c)],[C.model.dim.state,1]);
               end
            end


            %%%%%%%%%%% Create casadi-Function to represent the collocation polynomial p (and dp) to use on each interval:

            %%% Use collocation variables of first step as to define casadi-function:

            % Prepare args for casadi function:
            tau = casadi.SX.sym('tau');
            args = C.model.args; % Start with the standard (state,algeb,input,param) that the dynamics require
            args.tau = tau; % we add this here, so that it is easy to evaluate the polynomial at any tau later. (not par of auxiliary variables...)
            args.aux_var = structor.subvec(C.aux_var,C.aux_var.str.inc_1);

            % Create matrix whose columns are the collocation variables
            col_var_matrix = [args.state reshape(args.aux_var,C.model.dim.state,C.collocation.d)];
            

            % Create collocation polynomial p
            tau_powers = transpose(tau.^(C.collocation.d:-1:0));
            p.expr = col_var_matrix * C.collocation.Li *  tau_powers;

            % Create casadi function
            infields = fieldnames(args);
            args.out = p.expr;
            p.F = casadi.Function('F_collocation_p',args,infields,{'out'});


            % Create collocation polynomial derivative; dp
            tau_powers = transpose(tau.^(C.collocation.d-1:-1:0));
            dp.expr = col_var_matrix * C.collocation.dLi * tau_powers;

            % Create casadi function
            args.out = dp.expr;
            dp.F = casadi.Function('F_collocation_dp',args,infields,{'out'});


            %%%%%%%%%% Create casadi-Functions for p and dp at each step
            args = struct;            % used for calling the polynomials p,dp
            args_dyn = C.model.args;  % used for calling the dynamics

            for inc = 1:C.n_increments
               args.aux_var = structor.subvec(C.aux_var,C.aux_var.str.(['inc_',num2str(inc)])); % get auxiliary variables for current increment
               args.state = args_dyn.state; % initial condition for increment
               for c = 1:C.collocation.d
                  
                  % Update the tau (evaluation point) and colloctaion
                  % variable (state) used to evaluate dynamcis.
                  args.tau = C.collocation.tau(c+1);
                  args_dyn.state = C.aux_var.str.(['inc_',num2str(inc)]).state.(['samp_',num2str(c)]);

                  % collocation constraint
                  aux_expr.str.(['inc_',num2str(inc)]).(['samp_',num2str(c)]).dynamics = C.model.dynamics.call(args_dyn).out*dt - dp.F.call(args).out;
               end

               % Update the state variable to be the next state:
               args.tau = 1; % evaluate the polynomial at tau=1 to find next state
               args_dyn.state = p.F.call(args).out; % update the state to be the next state
            end

            % Finally, the final/('next') state variable is:
            next_state = args_dyn.state;

         end


      end
   end



   methods
      function args = example_args_Next(C)
         args = C.model.example_args;
         args.Dt = 1;
         args.aux_var = rand(C.aux_var.len,1);
      end
   end
end