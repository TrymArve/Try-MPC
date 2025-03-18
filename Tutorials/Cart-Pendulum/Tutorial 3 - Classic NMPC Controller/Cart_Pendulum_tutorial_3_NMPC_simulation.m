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


%% Defien problem:

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
init_state = structor;
init_state.str.x   = 1;
init_state.str.th  = 0;
init_state.str.dx  = 0;
init_state.str.dth = 0;
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

C.solve("ipopt",initial_state=init_state.vec);
C.display_optimization;

%{
We have now solved a single optimization problem, representing the MPC
prediction at the first time instance. In an MPC setting, this problem is
solved once for every 
%}

%% 