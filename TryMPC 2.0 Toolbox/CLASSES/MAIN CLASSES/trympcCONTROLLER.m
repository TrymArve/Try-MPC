classdef (Abstract) trympcCONTROLLER < handle

   properties
      Name (1,1) string
      ID (1,1) string
   end

   methods
      function C = trympcCONTROLLER(Name)
         C.Name = Name;
         C.ID = TRYMPC2.generate_id;
      end
   end


   methods
      function u = control(C,time,state)
         arguments(Input)
            C
            time (1,1) double
            state (:,1) double
         end
         arguments(Output)
            u function_handle % on form @(t,x) to match the ode45 solver
         end
         u = C.internal_control(time,state);
      end
   end

   methods(Abstract,Access=protected)
      u = internal_control(time,state)
   end
end