%% Cart-Pendulum Tutorial 3 - Simulating NMPC controller

%{
This tutorial is meant to be run one section at a time; (ctrl + enter) or
(ctrl + shift + enter).

Before running a section, read through the comments.
After each section is run, take the time to inspect what happened. Look at
plots, and use the command window to see what the class instance has
stored.

Note: to guarantee that everything runs as intended, only rund the section
in order from start to finish, without anything else in between.
%} 


% NB: It is recommended that you go through tutorial 2 first, to get
% familiar with the optimization framework.


%% Define problem:

fresh % start fresh

C = TRYMPC('Tester Instance',...
    state = ["x" "th" "dx" "dth"], ... % provide state variable names
    input = "ux", ...                  % provide input variable names
    param = ["L" "g" "mx" "mth"]);      % provide parameter names (keep them symbolic, to avoid hard-coding them)

% Dynamics:
u_pendulum = 0;
dynamics = @(s,a,i,p) [s.dx;
                       s.dth;
         	           -(p.L*i.ux + u_pendulum*cos(s.th) + p.L^2*s.dth^2*p.mth*sin(s.th) - p.L*p.g*p.mth*cos(s.th)*sin(s.th))/(p.L*(- p.mth*cos(s.th)^2 + p.mx + p.mth));
                       -(p.mx*u_pendulum+ p.mth*u_pendulum + p.L*p.mth*i.ux*cos(s.th) - p.L*p.g*p.mth^2*sin(s.th) + p.L^2*s.dth^2*p.mth^2*cos(s.th)*sin(s.th) - p.L*p.g*p.mx*p.mth*sin(s.th))/(p.L^2*p.mth*(- p.mth*cos(s.th)^2 + p.mx + p.mth))];

C.def_dynamics(dynamics)                                  % define dynamics
C.def_objective(quadratic=["Q","R","dR","Q_terminal"])    % define objective
C.def_integrator("ERK4","n_increments",2)                 % define integrator
C.def_stage_constraints("upper_bounds","ux")              % define stage constraints
C.def_horizon(50)                                         % define horizon (also build DOP on that horizon length)

% Define LaTex Display names for nice plots:
C.plotting.display_names.state.x = "$x$";
C.plotting.display_names.state.th = "$\theta$";
C.plotting.display_names.state.dx = "$\dot{x}$";
C.plotting.display_names.state.dth = "$\dot{\theta}$";
C.plotting.display_names.input.ux = "$u_{x}$";


%% Provide numeric values of your problem
% Parameters
C.parameters.str.L = 1;
C.parameters.str.g = 9.81;
C.parameters.str.mx = 5;
C.parameters.str.mth = 3;
% Initial state (cart-offset by 1)
C.initial_state.str.x   = 1;
C.initial_state.str.th  = 0;
C.initial_state.str.dx  = 0;
C.initial_state.str.dth = 0;

% % Fail for this initial condition:
% C.initial_state.str.x   = 1.0019*0;
% C.initial_state.str.th  = 0.0045;
% C.initial_state.str.dx  = -0.0099*0;
% C.initial_state.str.dth = 0.0039;

% quadratic cost penalities
C.quadratic_cost.Q = [1 10 0.1 0.1]; % state penalty
C.quadratic_cost.R = 1;              % input penalty
C.quadratic_cost.dR = 1;             % change in input penalty
C.quadratic_cost.Q_terminal = C.quadratic_cost.Q*350; % terminal state penalty
% upper bound on input
C.bounds.upper.ux = 0.4; % limit at which the algorithm fails... (works: 0.430415446809)

% Define the horizon length in terms of how far into the future you want to
% predict (given in seconds):
C.set_T(7); % predict a specific number of seconds into the future



%% Let's inspect the solution:
% C.clear_archive; close all
C.set_SQP_settings("tolerance_lagrangian",1000)
C.solve("sqp",max_iterations=10);
C.display_optimization;

%{
We have now solved a single optimization problem, representing the MPC
prediction at the first time instance. In an MPC setting, this problem is
solved once for every time step (every sample period / control cycle).
%}

%% Now we try to simulate the system, using the problem formulation we created above:

duration = 20; % duration of simulation

% This may take some time (up to a few minutes, depending on your computer)
C.simulate(duration,controller_type="NMPC_ipopt",sampling_time=0.1)
C.display_simulation("mark_samples",false);


%{
The reason it takes some time is
   1) NMPC is computationally demanding, and each optimization problem is
   solved fully at each time step.

   2) The implementation stores a lot of data during simulation, in order
   to analyze, inspect, and debug the simulation and optimization problems
   that were solved.

   3) Not relevant to this simulation, but nice to know:
   The simulator is quite general, and allows lots of different methods
   and settings. Even using external controllers. In particular, you can
   add "augmenting controllers", which are controllers u = F(state), which
   are assumed to be part of the systems dynamics, such that the new
   augmented system is essentially just simulated without control. This
   is the reason for the "archiving input signals" process after the
   simulation, which is also quite heavy. This process computes all the
   input signals that appeared during simulation, such that the controller
   signals can be analyzed. 
   Note; the "ode***()" MATLAB ode-solvers use error correction, thus
   variable step size. Therefore, it is hard to know exactly what
   input signals were used throughout the integration. There are several
   ways of keeping track of these signals, and the one used in this
   class simply does a new pass over the simulated state trajectory and
   re-computes the controller signals.
   There also the ability to add input disturbance in the form of agumented
   dynamics (function handle), for which all of the above is also true.
   Therefore when using input disturbance, the extra time is doubled.


Depending on when you are going through this tutorial, there might be an
option to run a "light" (or some other name) simulation, that does not store
as much data, and does not bother to reconstruct input signals, and thus
runs faster. Although, this will come at the cost of deminished ability to
analyze and debug after the fact.
%}

%% Or plot with sample markers:

% Marks every time-point and correcponding state value that was measured at
% that time. 
C.display_simulation("mark_samples",true);

%% A lot of information about the simulation is now available:

% the simulation is archived in the same way as the optimization problems,
% and each new simulation is appended to the end to the cell array:
C.archive.simulations{end}


%% Even the individual optimization problems within each simulation are logged:

C.archive.simulations{end}.optimizations

%% Let's look at the 3. optimizaiton problem that was solved:

C.archive.simulations{end}.optimizations{3}

%% The predicted state trajectory at that time was:

% the predicted state trajectory at the third sample-time: t = 0.3s
C.archive.simulations{end}.optimizations{3}.decision.str.state

%% We can easily plot this optimization problem using the familiar "display_optimization" method

C.display_optimization("optimizations",C.archive.simulations{end}.optimizations,... % specify the "optimizaitons" cell array where we find the problem to plot
                       optimization_number=50); % specify that the third problem should be plotted
                       


%% We can even the prediction at multiple sample-times simultaneously using the familiar multiplot method:

C.display_optimization("optimizations",C.archive.simulations{end}.optimizations,...
                       optimization_number=[3 25 50 75 100 150 200],... % specify the problems at sample-times 3, 25, 50, 75, 100, 150, and 200
                       multiplot="ontop",...    % specify to plot predictions in the same axes
                       color_match="solution",... % specify that colors should match across solution, rather than within axes.
                       linestyle="-"); % specify the linestyle








%% We now attempt to simulate using less requent sampling:

% first we define a halting condition. That is, if this condition becomes
% positive, the simulation stops/halts. This can be used to detect that a certain
% maneuver failed. For example; if the pendulum falls over during the
% offset correction maneuver.
my_halt_condition = @(s,i,p) abs(s(2)) - pi/100; % unlike the other function_handles, these arguments are vectors instead of structs.

% Simulate with the slower controller:
C.simulate(duration,controller_type="NMPC_ipopt",sampling_time=0.145,...
           halting_condition = my_halt_condition); % apply haulting condition

C.display_simulation("mark_samples",false);

% as we can see, the slower NMPC controller is unable to complete the
% maneuver, and the simulation stops already at around 2.5 seconds because
% the halting condition was met. That is, the pendulum fell to far.


%% What did the predictions look like?

[tiles,Layout] = C.display_optimization("optimizations",C.archive.simulations{end}.optimizations,...
                       optimization_number=[1 5 10 14],... % specify the problems at sample-times 1, 5, 10, 15, and 17
                       multiplot="ontop",...    % specify to plot predictions in the same axes
                       color_match="plot"); % specify that colors should match across plots/axes

%{
Notice that we can retrieve the axes handles as an output of the
"display_optimization" function. All axes handles are stored in the struct
"tiles", and can be edited as you wish. Also the "tiledLayout" object that
they are contained in is provided via the second output, which we save to the
"Layout" variable. 
%}

%% Look, we can edit the plots!

tiles.state.x.YLabel.String = "$x$ (cart pos.)";



%% By passing these handles on to another display-funciton, we can plot more stuff in the same plots!

% Let's add the simulation, to better understand what happened:
C.display_simulation(tiles=tiles,...    % pass on the "tiles" struct with all axes handles from above
                     colors="black",... 
                     linestyle="--",...
                     legend="simulation")


%% Furthermore, we can provide a reference to the NMPC controller.

%{
The default reference is the origin, but a reference trajectory can be
provded manually. Use the "C.set_ref()" method to specify the reference for
the state, algebraic states, or input. 
The reference must be given as a setpoint-vector, a series of
column-vectors consisting of the time-stamp as the top element, the
reference vector as the remaining elements, or a "function_handle" that
takes a time or row-vector of times as argument and returns a matrix whose
columns are the corresponding reference vectors.
%}


% Let's define a set-point, such that the cart only has to move slightly;
% to 0.8
C.set_ref("state",[0.8 0 0 0]')


%% What happens if we try to simulate with this reference instead?

C.simulate(duration,controller_type="NMPC_ipopt",sampling_time=0.1);
eq_sim = length(C.archive.simulations);
C.display_simulation("mark_samples",false);

% Look at that! It converged to the new reference instead of the origin!







%% Let's Simulate using the built-in SQP algorithm:

C.set_ref("state",[0.8 0 0 0]')
my_halt_condition = @(s,i,p) abs(s(2)) - pi/100; 

% Simulate with the SQP-NMPC controller:
C.simulate(duration,controller_type="NMPC_sqp",sampling_time=0.1,...
           halting_condition = my_halt_condition); % apply haulting condition

C.display_simulation("mark_samples",false);


%% Let's see its predictions along the way:

tiles = C.display_optimization("optimizations",C.archive.simulations{end}.optimizations,...
                       "multiplot","ontop",...
                       "optimization_number",[1:20:100 150 200]);

%% And again add the simulation itself
C.display_simulation(tiles=tiles,mark_samples=false,colors="red",opacity=1,linewidth=2.5,linestyle="--");







%% Lastly, we will generate an open-loop optimal offset-correction maneuver, and try to track that trajectory


%%%% We start over, and make a high resolution, long horizon problem:

fresh % start fresh

C = TRYMPC('Tester Instance',...
    state = ["x" "th" "dx" "dth"], ... % provide state variable names
    input = "ux", ...                  % provide input variable names
    param = ["L" "g" "mx" "mth"]);      % provide parameter names (keep them symbolic, to avoid hard-coding them)

% Dynamics:
u_pendulum = 0;
dynamics = @(s,a,i,p) [s.dx;
                       s.dth;
         	           -(p.L*i.ux + u_pendulum*cos(s.th) + p.L^2*s.dth^2*p.mth*sin(s.th) - p.L*p.g*p.mth*cos(s.th)*sin(s.th))/(p.L*(- p.mth*cos(s.th)^2 + p.mx + p.mth));
                       -(p.mx*u_pendulum+ p.mth*u_pendulum + p.L*p.mth*i.ux*cos(s.th) - p.L*p.g*p.mth^2*sin(s.th) + p.L^2*s.dth^2*p.mth^2*cos(s.th)*sin(s.th) - p.L*p.g*p.mx*p.mth*sin(s.th))/(p.L^2*p.mth*(- p.mth*cos(s.th)^2 + p.mx + p.mth))];

C.def_dynamics(dynamics)                                  % define dynamics
C.def_objective(quadratic=["Q","R","dR","Q_terminal"])    % define objective
C.def_integrator("IRK4 (L-stable)","n_increments",3)      % define integrator
C.def_stage_constraints("upper_bounds","ux")              % define stage constraints
C.def_horizon(100)                                        % define horizon (also build DOP on that horizon length)

% Define LaTex Display names for nice plots:
C.plotting.display_names.state.x = "$x$";
C.plotting.display_names.state.th = "$\theta$";
C.plotting.display_names.state.dx = "$\dot{x}$";
C.plotting.display_names.state.dth = "$\dot{\theta}$";
C.plotting.display_names.input.ux = "$u_{x}$";

% Parameters
C.parameters.str.L = 1;
C.parameters.str.g = 9.81;
C.parameters.str.mx = 5;
C.parameters.str.mth = 3;

% Initial state (cart-offset by 1)
C.initial_state.str.x   = 1;
C.initial_state.str.th  = 0;
C.initial_state.str.dx  = 0;
C.initial_state.str.dth = 0;

% quadratic cost penalities
C.quadratic_cost.Q = [1 10 0.1 0.1]; % state penalty
C.quadratic_cost.R = 1;              % input penalty
C.quadratic_cost.dR = 1;             % change in input penalty
C.quadratic_cost.Q_terminal = C.quadratic_cost.Q*350; % terminal state penalty
% upper bound on input
C.bounds.upper.ux = 0.4; 

C.set_T(30); % predict a specific number of seconds into the future


%% Generate an open-loop optimal maneuver:

C.solve("sqp");
C.display_optimization;

%% Now use the solution to create a reference:

state_trajectory = C.archive.optimizations{end}.decision.str.state;
input_trajectory = C.archive.optimizations{end}.decision.str.input;
time = C.archive.optimizations{end}.time;

C.set_ref("state",[time; state_trajectory],...            % set the state reference
          "input",[time; input_trajectory([1:end end])])  % set the input reference


%% And then we create the problem we will use for NMPC trakcing controller

C.clear_problem % clear the current problem, so we can make the NMPC problem

% After clearing the problem, the following must be re-defined:
C.def_objective(quadratic=["Q","R","dR","Q_terminal"])    % define objective
C.def_integrator("ERK4","n_increments",1)                 % define integrator
C.def_stage_constraints("upper_bounds","ux")              % define stage constraints
C.def_horizon(5)                                          % define horizon (also build DOP on that horizon length)

C.set_T(2); % choose a short prediction horizon


%% Now simulate the NMPC trakcing problem:

my_halt_condition = @(s,i,p) abs(s(2)) - pi/100; 

C.simulate(27,controller_type="NMPC_sqp",...
           sampling_time=0.1,...
           halting_condition=my_halt_condition)
tiles = C.display_simulation;


%% Let's look at some of the predictions it made along the way:

C.display_optimization(tiles=tiles,...  % provide the axes handles of the simulation plot
                       optimizations=C.archive.simulations{end}.optimizations,... % pass the optimizations cell array of the relevant simulation
                       optimization_number=1:20:270); % plot predictions made at these sample times
                       

%% Also plot the reference ontop

Ref = structor;
Ref.str.state = state_trajectory;
Ref.str.input = input_trajectory([1:end end]);

C.display_trajectory(Ref,time,tiles=tiles,colors="black",linewidth=1.5);

% We can now see that the solution seems to oscillate around the reference

%% Now we try to increase the resolution, integrator increments and prediction horizon to improve performance

C.clear_problem % clear the current problem

% After clearing the problem, the following must be re-defined:
C.def_objective(quadratic=["Q","R","dR","Q_terminal"])    % define objective
C.def_integrator("ERK4","n_increments",2)                 % define integrator
C.def_stage_constraints("upper_bounds","ux")              % define stage constraints
C.def_horizon(10)                                         % define horizon (also build DOP on that horizon length)

C.set_T(3); % choose a short prediction horizon

%% Now simulate again

my_halt_condition = @(s,i,p) abs(s(2)) - pi/100; 

C.simulate(27,controller_type="NMPC_sqp",...
           sampling_time=0.1,...
           halting_condition=my_halt_condition)


%% plot the reference first:
tiles_2 = C.display_trajectory(Ref,time,colors="black",linewidth=1.5);

%% Then plot the simulation result on top
C.display_simulation(tiles=tiles_2);


% It seems to follow fairly well!

