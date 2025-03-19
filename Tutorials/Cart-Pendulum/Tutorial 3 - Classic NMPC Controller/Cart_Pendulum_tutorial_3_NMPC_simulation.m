%% Cart-Pendulum Tutorial 2 - Building dynamic optimization problems and optimizing


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
C.def_horizon(50)                                        % define horizon (also build DOP on that horizon length)

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
% quadratic cost penalities
C.quadratic_cost.Q = [1 10 0.1 0.1]; % state penalty
C.quadratic_cost.R = 1;              % input penalty
C.quadratic_cost.dR = 1;             % change in input penalty
C.quadratic_cost.Q_terminal = C.quadratic_cost.Q*350; % terminal state penalty
% upper bound on input
C.bounds.upper.ux = 0.4;

% Define the horizon length in terms of how far into the future you want to
% predict (given in seconds):
C.set_T(7); % predict 30 seconds into the future



%% Let's inspect thte solution:

C.solve("ipopt");
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

%% We can then plot this



%% We can attempt to simulate using less requent sampling:

% first we define a halting condition. That is, if this condition becomes
% positive, the simulation stops/halts. This can be used to detect that a certain
% maneuver failed. For example; if the pendulum falls over during the
% offset correction maneuver.
my_halt_condition = @(s,i,p) abs(s(2)) - pi/10; % unlike the other function_handles, these arguments are vectors instead of structs.

% Simulate with the slower controller:
C.simulate(duration,controller_type="NMPC_ipopt",sampling_time=0.15,...
           halting_condition = my_halt_condition); % apply haulting condition

C.display_simulation("mark_samples",false);

% as we can see, the slower NMPC controller is unable to comlete the
% maneuver, and the simulation stops already at around 2.5 seconds because
% the halting condition was met.





%% We can also provide a non-zero reference for the NMPC controller:










%% Lets Simulate using the built-in SQP algorithm:

my_halt_condition = @(s,i,p) abs(s(2)) - pi/10; % unlike the other function_handles, these arguments are vectors instead of structs.

% Simulate with the SQP-NMPC controller:
C.simulate(duration,controller_type="NMPC_sqp",sampling_time=0.1,...
           halting_condition = my_halt_condition); % apply haulting condition

C.display_simulation("mark_samples",false);