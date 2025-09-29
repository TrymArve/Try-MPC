classdef TRYMPC2


   %%% BASIC PROPERTIES
   properties(SetAccess=immutable)
      Name (1,1) string
   end
   properties(Hidden)
      ID (1,1) string
   end

   

   %%% Contructor:
   methods
      function C = TRYMPC2(class_instance_name)
         arguments
            class_instance_name (1,1) string
         end

         C.Name = class_instance_name; % Apply name:
         C.ID = C.generate_id; % Create (probably) unique ID for class intance:

      end
   end


   %%% Static Helper functions:
   methods(Access=public, Static)

      %%% generate (probably) unique ID:
      function id = generate_id
         characters = ['A':'Z', 'a':'z', '0':'9'];
         id = string(characters(randi(numel(characters), [1, 6])));
      end


      %%% Convert to seconds:
      function out = minutes(m)
         out = m*60;
      end
      function out = hours(h)
         out = h*3600;
      end
      function out = days(d)
         out = d*24*3600;
      end


      %%% Error messaging:
      function usererror(text)
         arguments
            text (1,:) char = ''
         end
         error(['USER ERROR: ',char(text),' \newline \tab --- Try ... for help.'])
      end
   end
end
