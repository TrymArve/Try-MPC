classdef trympcNUMERIC_MODEL
   properties(SetAccess=immutable)
      % Give the model instance a name
      Name (1,1) string 
   end

   properties
      % Specific parameter values for a model
      parameters structor = structor;

      % A suggested initial state for a model (for the given parameters)
      initial_state structor % structor

      % A reference to track (function_handle: @(t))
      ref struct = struct % fieldsnames should be members of ["state","input","output"]
   end


   properties
      % A stable equilibrium of a model (for the given parameters)
      stable_equilibrium_state % structor
      stable_equilibrium_input % structor
      stable_equilibrium_output % structor

      % An unstable equilibrium of a model (for the given parameters)
      unstable_equilibrium_state % structor
      unstable_equilibrium_input % structor
      unstable_equilibrium_output % structor
   end


   methods
      function C = trympcNUMERIC_MODEL(Name,param,options)
         arguments
            Name (1,1) string

            % model parameters 
            param structor

            % A suggested initial state for a model (for the given parameters)
            options.initial_state structor


            % references  struct with fields "state", "input" and/or "output" ( @(t) a function of time )
            options.ref struct



            % A stable equilibrium of a model (for the given parameters)
            options.stable_equilibrium_state structor
            options.stable_equilibrium_input structor
            options.stable_equilibrium_output structor

            % An unstable equilibrium of a model (for the given parameters)
            options.unstable_equilibrium_state structor
            options.unstable_equilibrium_input structor
            options.unstable_equilibrium_ouput structor
         end

         C.Name = Name;
         C.parameters = param;

         Props = ["initial_state",...
                  "stable_equilibrium_state",...
                  "stable_equilibrium_input",...
                  "unstable_equilibrium_state",...
                  "unstable_equilibrium_input",...
                  "ref",...
                  ];


         for prop = Props
            if isfield(options,prop)
               C.(prop) = options.(prop);
            end
         end

      end
   end
end