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


% NB: It is recommended that you go through tutorial 1 first.


%% Preview:


%{
First off, notice that you don't really need to know everything that is
explained here. Much of this is to give a basic understanding of the
sturcture of how the probimization problem is defined, built , stored and
solved.

Mostly you only need to know how to define the problem, and the TRYMPC
handles everything for you. Therefore to not deter you, we give you a
preview of what you actually need to know yourself, in order to start
optimizing.
%}


%%% A quick demo of how to optimize:
fresh % start fresh

C = TRYMPC('Tester Instance',...
    state = ["x" "th" "dx" "dth"], ... % provide state variable names
    input = "ux", ...                  % provide input variable names
    param = ["L" "g" "mx" "mth"]);      % provide parameter names (keep them symbolic, to avoid hard-coding them)

%%% Define dynamics:
u_pendulum = 0;
dynamics = @(s,a,i,p) [s.dx;
                       s.dth;
         	           -(p.L*i.ux + u_pendulum*cos(s.th) + p.L^2*s.dth^2*p.mth*sin(s.th) - p.L*p.g*p.mth*cos(s.th)*sin(s.th))/(p.L*(- p.mth*cos(s.th)^2 + p.mx + p.mth));
                       -(p.mx*u_pendulum+ p.mth*u_pendulum + p.L*p.mth*i.ux*cos(s.th) - p.L*p.g*p.mth^2*sin(s.th) + p.L^2*s.dth^2*p.mth^2*cos(s.th)*sin(s.th) - p.L*p.g*p.mx*p.mth*sin(s.th))/(p.L^2*p.mth*(- p.mth*cos(s.th)^2 + p.mx + p.mth))];

C.def_dynamics(dynamics)                                  % define dynamics
C.def_objective(quadratic=["Q","R","dR","Q_terminal"])    % define objective
C.def_integrator("ERK4","n_increments",2)                 % define integrator
C.def_stage_constraints("upper_bounds","ux")              % define stage constraints
C.def_horizon(100)                                        % define horizon (also build DOP on that horizon length)

%%% Provide numeric values of your problem
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
C.set_T(30); % predict 15 seconds into the future

% Solve optimizaiton problem using ipopt:
C.solve("ipopt",initial_state=init_state.vec);

% Display solution:
C.display_optimization;



%%%%%%%%%%%%%%%%%%%%%% END of demo








%% Acutal beginning of Tutorial: Prepare a class instance and define dynamics, as in tutorial 1:

fresh % start fresh

C = TRYMPC('Tester Instance',...
    state = ["x" "th" "dx" "dth"], ... % provide state variable names
    input = "ux", ...                  % provide input variable names
    param = ["L" "g" "mx" "mth"]);      % provide parameter names (keep them symbolic, to avoid hard-coding them)


%%% Define dynamics:
u_pendulum = 0;
dynamics = @(s,a,i,p) [s.dx;
                       s.dth;
         	           -(p.L*i.ux + u_pendulum*cos(s.th) + p.L^2*s.dth^2*p.mth*sin(s.th) - p.L*p.g*p.mth*cos(s.th)*sin(s.th))/(p.L*(- p.mth*cos(s.th)^2 + p.mx + p.mth));
                       -(p.mx*u_pendulum+ p.mth*u_pendulum + p.L*p.mth*i.ux*cos(s.th) - p.L*p.g*p.mth^2*sin(s.th) + p.L^2*s.dth^2*p.mth^2*cos(s.th)*sin(s.th) - p.L*p.g*p.mx*p.mth*sin(s.th))/(p.L^2*p.mth*(- p.mth*cos(s.th)^2 + p.mx + p.mth))];

C.def_dynamics(dynamics)




%% We start by defining the objective

% We choose a quadratic cost in the states ("Q" is the state penalty vector)
C.def_objective(quadratic="Q") 

%% Have a look in the "cas" property, to find the objective-related expressions:

% A casadi.SX variable was automatically generated to represent the state
% penalty weights:
C.cas.objective.weights.Q 

% Also, a casadi.Function was generated to evaluate the total stage cost;
% (which in this case only consists of the quadratic state penalty Q)
C.cas.objective.cost.stage_total.F


%% We then define the integrator:

% define an explicit euler inetgrator, with 3 integration steps ("increments")
C.def_integrator("Explicit Euler","n_increments",3)

%{
 Note: the integration "increments" are different from the number of state
 variables on the horizon of your dynamic NLP. F.ex. in this case, we have
 3 integration steps ("increments") between state(t_1) and state(t_2)
%}

% Have a look at the newly created properties:
C.cas.integrator

% Notably, we now have a casadi.Function that takes in a state and control
% variable, and output the next state (at a time Dt into the future):
C.cas.integrator.next.F

%% Now we may use the integrator to define a series of state and input variables on a horizon that obeys our model:

% horizon length in terms of number of stages to evaluate:
N = 50;

%{
Note: A "stage" is a point in time, at which we represent our system state and control input by
variables, and evaluate the constraints and objective.

In this case, we chose to have N stages into the future, which results in
N+1 total stages when including the initial stage (at time = 0).
%}

% and simply build the dynamic optimization problem (DOP)
C.def_horizon(N)

%% Now Look at all the new expressions!

% All decision variables (a.k.a. "optimization variables" or "primal variables")
C.cas.horizon.decision.str.state
C.cas.horizon.decision.str.input

%% Note that we can easily obtain the decision vector by:
C.cas.horizon.decision.vec

%% Also we obtained dynamic constraints (symbolic):
C.cas.horizon.constraints.equality.str % The dynamic constraints that ensure integration from stage k to k+1.

%% F.ex. look at the dynamics constraint that ensure integration from stage 3 to 4:
C.cas.horizon.constraints.equality.str.stage_3.dynamic.continuation

%{
Note that there can be other constraints at a stage, such as bounds, thus
the dynamic constraints are sotred into a "dynamic" field.

Furthermore, the dynamic constraints can consist of more than just
continuation constraints, thus the continuation constraints are sorted as
well.
("continuation constraints" are the constraints that forces the "shooting gaps" to close)
%}

%% Lastly, we emphasize that the enitre symbolic equality constrainst vector can be obtained as
C.cas.horizon.constraints.equality.vec


%% Now look at the "problem" property:

% This property contains the dynamic optimization problem:
C.cas.problem
C.cas.problem.decision.str  % The decision structor
C.cas.problem.objective     % The objective function (.expr for symbolic expression, and .F for casadi.Function)


%% This is all we need to start optimizing the cart pendulum!

% Start by defining an offset correction problem (Initial state: the pendulum is upwards, but the cart is offset by 1)
init_state = structor;
init_state.str.x   = 1;
init_state.str.th  = 0;
init_state.str.dx  = 0;
init_state.str.dth = 0;

% Define the numerical values for parameters
C.parameters.str.L = 1;
C.parameters.str.g = 9.81;
C.parameters.str.mx = 5;
C.parameters.str.mth = 3;

% Define the quadratic state penalty weight Q:
C.quadratic_cost.Q = [1 10 0.1 0.1];

% And define the the prediction horizon lenth in seconds
C.set_T(5)

% And simply optimize:
C.solve("ipopt",initial_state=init_state.vec);

%% Now look at the solution:

% Define LaTex Display names for nice plots:
C.plotting.display_names.state.x = "$x$";
C.plotting.display_names.state.th = "$\theta$";
C.plotting.display_names.state.dx = "$\dot{x}$";
C.plotting.display_names.state.dth = "$\dot{\theta}$";
C.plotting.display_names.input.ux = "$u_{x}$";

% and plot solution:
C.display_optimization;



%% We notice that the control signal is very large, perhaps we should put some constrainst on it:

% we can easily define input bounds and rebuild the problem:

% start by clearing the problem:
C.clear_problem

% Redefine the objective and integrator:
C.def_objective(quadratic="Q") 
C.def_integrator("Explicit Euler","n_increments",3)

% Now, define stage constraints (in this case; "bounds" on the input signal)
C.def_stage_constraints("lower_bounds","ux","upper_bounds","ux")

% And provide a numerical value to the bounds
C.bounds.lower.ux = -5;
C.bounds.upper.ux =  5;

% rebuild problem:
C.def_horizon(N);

% And simply optimize:
C.solve("ipopt",initial_state=init_state.vec);
C.display_optimization; % plot solution

%{
Notice that all the numerical values and display-names have not been
deleted, and we could simply re-optimize immediately.
%}

