classdef trympcSIMULATION



   properties(SetAccess = immutable)
      
      %%%%% SIMULATION RESULTS
      States
      Inputs
      Times

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

   properties(Dependent)
      states
      inputs
      times
   end





   methods
      function C = trympcSIMULATION(sim_model,initial_state,duration,sampling_time,controller,options)
         arguments
            %%% Fundamental Arguments:
            sim_model (1,1) function_handle % on form @(t,x,u) (used by "odeXX" as @(t,x) sim_model(t,x,K(t)))
            initial_state (:,1) {double,string} 
            duration (1,1) double {mustBePositive}
            sampling_time (1,1) double {mustBePositive}
            controller (1,1) {mustBeA(controller,["trympcCONTROLLER","double","function_handle"])}
            %{
The input must either be a controller class inherited from
trympcCONTROLLER --which produces function_handles-- or a function handle
directly. 
The function handle should be on the form @(time,state), i.e. only a function of
time and state, and produce a double (dim_input,1)
Alternatively, one may provide a 'double'-vector directly, in which the
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


         if isa(controller,'double')
            make_controller = @(~,~) (@(~,~) controller); % returns a function_handle, that only returns the constant input
            C.controller_name = "constant";
         elseif isa(controller,'function_handle')
            make_controller = @(~,~) controller; % just return the function handle that was privided.
            C.controller_name = "a function_handle";
         elseif isa(controller,'trympcCONTROLLER')
            make_controller = @controller.control; % turn "make_controller" into an alias for the "control" function of the controller.
            C.controller_name = controller.Name;
            C.controller_ID = controller.ID;
         end


         %%%%% SET CLASS PROPERTIES:
         % timing
         C.duration = duration;
         C.start_time = options.start_time;
         C.control_delay = options.control_delay;
         % solver
         C.ode_solver = options.ode_solver;
         C.ode_options = options.ode_options;


         % Prepare misc:
         odeXX = str2func(C.ode_solver);
         end_time = C.start_time + duration;
         sampling_time = min(sampling_time,duration);
         counter = 1;

         % Prepare variables:
         time = C.start_time;
         state = {initial_state, initial_state}; % [current_state, previous_state]
         
         

         % Simulate:
         while time < end_time

            K = make_controller(time,state{1+C.control_delay}); % produces a function handle @(t,s)
            [time_sim,state_sim] = odeXX(@(t,s) sim_model(t,s,K(t,s)), time+[0,sampling_time], state{1}, C.ode_options);

            % store simulation
            C.States{counter} = state_sim';
            C.Times{counter} = time_sim';
            C.Inputs{counter} = K(time_sim');

            
            % UPDATE:
            state{2} = state{1};        % store previous solution
            state{1} = state_sim(1,:)'; % update current solution
            time = time + sampling_time; % increment time
            counter = counter + 1; % increment counter
         end



      end


      function traj = get.states(C)
         traj = horzcat(C.States{:});
      end
      function traj = get.inputs(C)
         traj = horzcat(C.Inputs{:});
      end
      function traj = get.times(C)
         traj = horzcat(C.Times{:});
      end
   end





end