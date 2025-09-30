classdef trympcMODEL
   

   properties(SetAccess=immutable)

      %%% Admin:
      Name (1,1) string


      %%% Main properties:

      % variables
      state  % "x"
      input  % "u"
      param  % "p"

      % expressions
      dynamics      %    x_dot = f(x,u,p)
      output        %   y = output(x,u,p)

   end

   % Advanced properties
   properties(SetAccess = immutable)
      Jac % contains jacobians of all expressions wrt all variables
      dim % struct containing the dimensions of variables
   end

   % Numeric values
   properties(SetAccess=public)
      numeric_model (1,:) trympcNUMERIC_MODEL % structor containing the numeric values for the parameters
   end
   
   properties(Hidden,SetAccess=immutable)
      ID (1,1) string
      names

      % Auxiliary properties
      variable_types = ["state", "input"];
      variable_types_par = ["state", "input", "param"];
      expression_types

      % Struct holding the all variables, to be used for function calls
      args % in-arguments to the dynamics and algebraics

      % Fast versions: the same as dynamics/algebraics, but arguments are not in a struct, but positional instead. (faster and sometimes more convenient)
      dynamics_alt
      output_alt
   end






   methods
      function C = trympcMODEL(model_name,names,model)

         arguments
            model_name
            names.state (1,:) string
            names.input (1,:) string
            names.param (1,:) string = [];
            model.dynamics (1,1) function_handle
            model.output (1,1) function_handle
         end
         fprintf('Creating Model...  ');
         def_time = tic;
         
         % % Call superclass constructor first
         % C@TRYMPC2_overclass(model_name);
         C.Name = model_name;
         C.ID = TRYMPC2.generate_id;
         C.names = names;
         
         % Characterize 'types':
         C.expression_types = string(fieldnames(model))';

         
         % Ensure proper initialization:
         if ~isfield(names,'state')
            TRYMPC2.usererror('must provide a list of state names.')
         end
         if ~isfield(names,'input')
            TRYMPC2.usererror('must provide a list of input names.')
         end


         % Set model type
         if ~isfield(model,'dynamics')
            TRYMPC2.usererror('must define dynamics: function_handle with stucts as arguments. F.ex: "dynamics = @(s,a,i,p) s.x1^2 + 2*i.u1*p.mass" \newline The same goes for algebraics if relevant.')
         end



         % Define variables:
         for type = C.variable_types_par
            C.(type) = structor;
            C.dim.(type) = length(names.(type));
            for name = names.(type)
               C.(type).str.(name) = casadi.SX.sym(char(type + "_" + name));
            end
         end





         %%% Create dynamics / ouput
         for expr_type = C.expression_types


            try
               expr = model.(expr_type)(C.state,C.input,C.param);
            catch ME
               warning(['USER ERROR: the ',char(expr_type),' that was provided fails when evaluating it with casadi-variables (note that the arguments are structors, and one must use either s.str for the struct part, or s.vec for the vector part).'])
               disp(ME)
            end

            if expr_type == "dynamics"
               if size(expr,1) ~= C.dim.state || size(expr,2) ~= 1
                  TRYMPC2.usererror('The dimensions of the dynamics must be (n_state,1)')
               end
            elseif expr_type == "output"
               C.dim.output = size(expr,1);
               if size(expr,2) ~= 1
                  TRYMPC2.usererror('The dimensions of the output must be (n_output,1) (i.e. a column vector)')
               end
            end
            

            % Input arguments
            for type = C.variable_types_par
               C.args.(type) = C.(type).vec;
            end
            infields = fieldnames(C.args);

            % Output argument
            args = C.args;
            args.out = expr;

            % Create Casadi-function:
            C.(expr_type) = casadi.Function(char("F_"+expr_type),args,infields,{'out'});

            % create jacobians
            for type = C.variable_types_par
               args.out = jacobian(expr,C.(type).vec);
               C.Jac.(expr_type).(type) = casadi.Function(char("J_"+expr_type+"_wrt_"+type),args,infields,{'out'});
            end

         end



         % Generate fast versions:
         args = C.args;
         C.dynamics_alt   = casadi.Function('F_dynamics_fast',  {args.state,args.input,args.param},{C.dynamics.call(  args).out});
         if ismember(C.expression_types,"output")
            C.output_alt   = casadi.Function('F_output_fast',  {args.state,args.input,args.param},{C.output.call(  args).out});
         end

         disp(['done.  ',sec2str(toc(def_time)),' Name: "',char(C.Name),'"'])
      end
   end



   % Helpful methods
   methods
      function args = example_args(C)
         args.state = rand(C.dim.state,1);
         args.input = rand(C.dim.input,1);
         args.param = rand(C.dim.param,1);
      end
   end






   % Create help struct
   methods(Hidden,Static)
      function Help(field)
         arguments
            field (1,1) string
         end
         switch field
            case "state"
               disp("state - is the (CasADi) symbolic state vector of your system. It is a 'structor', which means that is acts both as a stuct and a vector. \newline ---- Try 'myMODEL.state.str' and 'myMODEL.state.vec'.")
            case "input"
               disp("input - is the (CasADi) symbolic input vector of your system. It is a 'structor' (see help('state') for more).")
            case "param"
               disp("param - is the (CasADi) symbolic parameter vector of your system. It is a 'structor' (see help('state') for more).")
            case "dynamics"
               disp("dynamics - is the system dynamic model. I.e. state_dot = dynamics(state,algeb,input,param). It is a CasADi-function.")
            case "Jac"
               disp("Jac - contains the jacobian of dynamics/algebraics w.r.t. the various variable-vectors. Try for example 'myMODEL.Jac.dynamics.state' or 'myMODEL.Jac.algebraics.input'")
            case "dim"
               disp("dim - contains the dimension of the variable vectors.")
         end
      end
   end
end

