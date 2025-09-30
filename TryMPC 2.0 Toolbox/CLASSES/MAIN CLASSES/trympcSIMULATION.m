classdef trympcSIMULATION



   properties(SetAccess = immutable)
      
      %%%%% SIMULATION RESULTS
      states
      inputs
      times

      %%%%% INFO ABOUT SIMULATION
      % Controller
      controller_name = []
      controller_ID = []
      % timing
      duration
      start_time
      control_delay
      % solver
      ode_solver
      ode_options


   end





   methods(Static)
      function C = trympcSIMULATION(sim_model,initial_state,duration,sampling_time,options)
         arguments
            %%% Fundamental Arguments:
            sim_model (1,1) funciton_handle % on form @(t,x,u) (used by "odeXX" as @(t,x) sim_model(t,x,K(t)))
            initial_state (:,1) {double,string} 
            duration (1,1) double {mustBePositive}
            sampling_time (1,1) double {mustBePositive}

            %%% Main Options:
            options.controller {mustBeA(options.controller,["trympcCONTROLLER","double","function_handle"])}
            %{
The input must either be a controller class inherited from
trympcCONTROLLER --which produces function_handles-- or a function handle
directly. 
The function handle should be on the form @(t), i.e. only a function of
time, and produce a double (dim_input,1)
Alternatively, one may provide a double vector directly, in which the
control signal is consant.
            %}
            
            %%% Simulation Settings:
            options.start_time (1,1) double = 0;
            options.ode_solver (1,1) string {mustBeMember(options.ode_solver,["ode45","ode23","ode113","ode78","ode89","ode15s","ode23s","ode23t","ode23tb","ode15i"])} = "ode45";
            options.ode_options (1,1) struct = odeset('RelTol',10^(-7),'AbsTol',10^(-7));
            options.control_delay (1,1) logical = false; 
            %{
CONTROL DELAY:
If control delay is set to true, then the a control delay of a full
sampling time is added to the simulation. This funcitons as a worst case
control delay, since the real control delay should be at most that large.
If the true contorl delay is smaller, then the performance should only be
better (presumably...). 
Note that; a real controller implementation can always be made to act as
the worst case simulation by adding artificial delay.
In code, this control delay is simulated by simply reporting the previous
measurement to the controller rather than the current measurement (which is reported next sample instead).
            %}
         end


         if isfield(options,'controller')
            if isa(options.controller,'double')
               make_controller = @(~,~) (@(~) options.controller); % returns a function_handle, that only returns the constant input
            elseif isa(options.controller,'function_handle')
               make_controller = @(~,~) options.controller; % just return the function handle that was privided.
            elseif isa(options.controller,'trymocCONTROLLER')
               make_controller = @options.controller.control; % turn "make_controller" into an alias for the "control" function of the controller.
            end
         end
         % Prepare misc:
         end_time = options.start_time + duration;
         sampling_time = max(sampling_time,duration);
         C.states = dictionary;
         C.inputs = dictionary;
         C.times = dictionary;
         counter = 1;

         % Prepare variables:
         time = options.start_time;
         state = {initial_state, initial_state}; % [current_state, previous_state]
         
         

         % Simulate:
         while time < end_time

            K = make_controller(time,state); % produces a function handle @(t)
            [time_sim,state_sim] = options.ode_solver(@(t,s) sim_model(t,s,K(t)), time+[0,sampling_time], state{1+C.control_delay}, options.ode_options);

            % store simulation
            C.states(counter) = state_sim';
            C.times(counter) = time_sim';
            C.inputs(counter) = K(time_sim');


            state{2} = state{1};        % store previous solution
            state{1} = state_sim(1,:)'; % update current solution

            counter = counter + 1;
         end





         %%%%% SET CLASS PROPERTIES:
         % Controller
         C.controller_name = options.controller.Name;
         C.controller_ID = optins.controller.ID;
         % timing
         C.duration = duration;
         C.start_time = options.start_time;
         C.control_delay = options.control_delay;
         % solver
         C.ode_solver = options.ode_solver;
         C.ode_options = options.ode_options;


      end
   end





end