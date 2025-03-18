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
C.set_T(30); % predict 30 seconds into the future

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

%% Now we can plot the solutions in the same figure to compare them

C.display_optimization("optimization_number",[1 2],...
                       "multiplot","separate") % Specify that we want the plots to appear in different axes (subfigures) 


%% Perhaps we should loosen the input bound a little:

% Try an upper bound of 10 instead:
C.bounds.upper.ux = 10;

% And simply optimize again:
C.solve("ipopt",initial_state=init_state.vec);
C.display_optimization; % plot solution

%% Now we can easily comare all three solutions

C.display_optimization("optimization_number",1:length(C.archive.optimizations),...
                       "multiplot","ontop") % Let's specify that they should appear on top of each other

%% Lookin gat the plots, we see that the trajectories are very coarse, so let's try finer integrations:

N = 150;

% rebuild problem:
C.def_horizon(N);

% And simply optimize again:
C.solve("ipopt",initial_state=init_state.vec);


%% Let's Seee how this changed the trajectories:

C.display_optimization("optimization_number",length(C.archive.optimizations)+(-1:0),...
                       "multiplot","ontop") % Let's specify that they should appear on top of each other

%% notice that the number of iterations and solve time did not change significantly:

disp('n. iterations:')
disp(['   -- (N =  50): ',num2str(C.archive.optimizations{end-1}.n_iterations)])
disp(['   -- (N = 150): ',num2str(C.archive.optimizations{end  }.n_iterations)])
disp('solve times:')
disp(['   -- (N =  50): ',sec2str(C.archive.optimizations{end-1}.solve_time.total)])
disp(['   -- (N = 150): ',sec2str(C.archive.optimizations{end  }.solve_time.total)])

%% Close all plots (if you have not already)
close all




%%%%%%%%%%%%%%%%%%%%%% We will not go more into the various settings, and possibilities:

%% Create a new problem with more complicated objective, integrator and constraints:

C.clear_problem

%% We go through them one by one, starting with objective:

%%%%%%%%%%%%%%%%%  Objective:
%{
 Specify:
    - "Q"  for state deviation (from reference, which is zero at the moment),
   (- "Z"  for algebraic state deviation. This is not relevant here. )
    - "R"  for input deviation,
    - "dQ" for penalty to change in state from one stage to the next ( Delta_state_k = state_(k) - state_(k-1) )
   (- "dZ" for penalty to algebraic state change. This is not relevant here. )
    - "dR" for penalty to input changes
    - "Q_terminal" for adding quadratic penalty to terminal state (deviation from terminal reference)
   (- "Z_terminal" for adding quadratic penalty to terminal algebraic state (deviation from terminal reference) )

Furthermore, one may manually add objective terms, but limited to stage-wise costs:
    - "stage_cost":    for adding an externally provided (not necessarily quadratic) stage cost
    - "terminal cost": for adding an externally provided (not necessarily quadratic) terminal cost
%}

% Let's go with quadratic terms on state, input, input-change, and terminal
% state.
quad_terms = ["Q","R","dR","Q_terminal"];

% and we also would like to add a term representing the sum of cart
% position (x) and pendulum angle (th):
non_quadratic_stage_cost = C.cas.state.x + C.cas.state.th;

% Now apply to TRYMPC:
C.def_objective(quadratic=quad_terms, stage_cost=non_quadratic_stage_cost);


%% Now we look at the integrator:


%{
We can easily define any RK method by proving the butcher tableau: 
    F.ex. Explicit RK4:
        - bucher_tableau_b = [1/6 1/3 1/3 1/6];
        - bucher_tableau_A = [ 0   0   0   0 ;
                              1/2  0   0   0 ;
                               0  1/2  0   0 ;
                               0   0   1   0 ];
        - note that the "c" vector automatically calculated as the sum of
          the rows of A

     To define an explicit RK integrator, simply use:
                " C.def_integrator("explicit butcher tableau","bucher_tableau_A",A,"bucher_tableau_b",b) "
     where A and b make up the butcher tableau.

     To define an implicit RK integrator, just use the "implicit butcher tableau" argument instead:
                " C.def_integrator("implicit butcher tableau","bucher_tableau_A",A,"bucher_tableau_b",b) "

     For examples of RK methods, see: https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods

Predifined Integrators:
  A few predefined integrators are available:
      - Explicit Euler         (multiple shooting)  (Forward Euler)
      - Implicit Euler         (fully simultaneous) (Backward Euler)
      - ERK4                   (multiple shooting)  (the "classic" RK method)
      - ERK4 (simultaneous)    (fully simultaneous) (A fully simultaneous version of the classic RK4 method)
      - etc...

        
%}

%%%%%%%%%%%%%%%%% Define Integrator:
% Let's use a collocation based integrator
C.def_integrator("collocation",...                       % Select "collocation" as the integration scheme
                 "collocation_polynomial_order",4,...    % Choose polynomial order 4 (often denoted "d")
                 "n_increments",3,...                    % Number of integration steps per stage
                 collocation_polynomial_type="legendre") % Choose polynomial type




%% We will add some more constraints this time:

%{

Similarly to the "def_objective", the "def_constraints" has built in the
ability to add bounds directly, simply by spcifying the argument
"lower_bounds" or "upper_bounds".

You may also add general stage-wise constraints by passing a struct of
expressions to either the "equality" or "inequality" argument

Similarly we may add terminal constraints via "def_terminal_constraint",
which works the same way, but only acts on the terminal state.

%}

%%%%%%%%%%%%%%%%% Stage Constraints:
% Manually provide a stage-wise constraint on the total kinetic energy of
% the system to be smaller than 0.12 J: (note that inequalityies are on the form: h(z) >= 0)
max_kinetic_energy = 0.12;
inequality.kinetic_energy = @(s,a,i,p) -(0.5*p.mx*s.dx^2 + 0.5*p.mth*(s.dth*p.L)^2)  + max_kinetic_energy;


% We also want to add bounds to the pendulum angle (th), and the control input (ux):
C.def_stage_constraints("lower_bounds",["th","ux"],"upper_bounds",["th","ux"],inequality=inequality)


%%%%%%%%%%%%%%%%% Terminal Constraints:
% Force the cart to be at zero, with the pendulum at rest in the upright
% position: (note that the cart may be moving sideways)
terminal_equality.cart_at_zero = @(s,a,i,p) s.x;
terminal_equality.pendulum_at_zero = @(s,a,i,p) s.th;
terminal_equality.pendulum_at_still = @(s,a,i,p) s.dth;
% Define terminal constrainst, but also add bounds on the terminal speed of
% the cart:
C.def_terminal_constraint("equality",terminal_equality,...  % Append the terminal equality constraints
                          "upper_bounds","dx",...           % Add upper bound on the terminal cart speed
                          "lower_bounds","dx")              % Add lower bound on the terminal cart speed



%% Now build the problem

N = 50;
C.def_horizon(N)


%% We can also inspect the problem that has been created:

% The problem is stored in:
C.cas.problem

% ... and some stats are stored in
C.display.problem

%% Display information about the problem

print_level = 4;  % Choose a level from 1 to 4
C.display_problem(print_level) 


%% Also we can inspect the various sparsity patterns of our problem!

% Check sparsity of the Hessian of the objective and Lagrangian (Hessian w.r.t. the decision variables (primal))
C.display_sparsity(["objective","Lagrangian"])


%% Constraint Jacobian

% And check out the sparsity of the Jacobian of the constraint vectors  (w.r.t. the decision variables (primal)):
C.display_sparsity(["equality","inequality"])


%% KKT matrix!

% Lastly we can inspect the KKT matrix!
C.display_sparsity("KKT")


%% We can rebuild our problem, but select "sorted" for the primal-dual vector:

% Selecting "sorted" makes the primal-dual vector sorted by the relevant
% time stage:
C.def_horizon(N,"primaldual","sorted")

% The KKT matrix then becomes diagonal !
C.display_sparsity("KKT")







%%%%%%%%%%%%%%%%%%%%%% It's time to solve our new problem !

%% Set numerical values:

% Parameters:
C.parameters.str.L = 1;
C.parameters.str.g = 9.81;
C.parameters.str.mx = 5;
C.parameters.str.mth = 3;

% Initial state:
init_state = structor;
init_state.str.x   = 1;
init_state.str.th  = 0;
init_state.str.dx  = 0;
init_state.str.dth = 0;

% Cost coefficients:
C.quadratic_cost.Q = [1 10 0.1 0.1];  % state deviation
C.quadratic_cost.R = 1;               % input deviation
C.quadratic_cost.dR = 1;              % input change
C.quadratic_cost.Q_terminal = C.quadratic_cost.Q*350; % terminal state deviation

% Bounds on pendulum angle:
C.bounds.lower.th = -0.5;
C.bounds.upper.th =  0.5;

% Bounds on input:
C.bounds.lower.ux = -1;
C.bounds.upper.ux =  1;

% Terminal bounds on cart speed:
C.terminal_bounds.lower.dx = -0.1;
C.terminal_bounds.upper.dx = 0.1;

% Set the horizon length
C.set_T(20);


%% Before solving, we specify some ipopt options:

% options for casadi's "nlpsol": (which can pass options to ipopt)
nlpsol_options.ipopt.print_level = 0;
nlpsol_options.ipopt.max_iter = 100;
nlpsol_options.ipopt.tol = 1e-9;
nlpsol_options.ipopt.acceptable_tol = 1e-6;
nlpsol_options.ipopt.compl_inf_tol = 1e-9;
nlpsol_options.ipopt.constr_viol_tol = 1e-9;
nlpsol_options.ipopt.dual_inf_tol = 1e-9;


%% Solve:

sol = C.solve("ipopt","initial_state",init_state.vec,"ipopt_nlpsol_options",nlpsol_options)
[tiles,Layout] = C.display_optimization;

%% Note that the "display_optimization" funciton outputs "tiles" and "Layout" which allows us to edit the figure:

% Each axes object can be accessed through the a handle found sorted into
% structs as follows:
tiles.state.x    % handle for the "x"   axes
tiles.state.dx   % handle for the "dx"  axes
tiles.state.th   % handle for the "th"  axes
tiles.state.dth  % handle for the "dth" axes
tiles.input.ux   % handle for the "ux" axes

% Also, the "tiledLayout" is available through the output "Layout":
Layout


%% Now, we plot the kinetic energy at each stage, and verify that it is smaller than 0.12 J

% One option is to retrieve the variables form the solution vector:
dx  = C.archive.optimizations{end}.decision.str.state(3,:); % state 3 is "dx"
dth = C.archive.optimizations{end}.decision.str.state(4,:); % state 4 is "dth"
mx  = C.parameters.str.mx;   % get parameter "mx"
mth = C.parameters.str.mth;  % get parameter "mth"
L   = C.parameters.str.L;    % get parameter "L"

% and compute the kinetic enegery manually:
K = 0.5*mx*dx.^2 + 0.5*mth*(dth*L).^2; % compute the kinetic energy

% and plot it manually:
figure(Name='Kinetic Energy')
ax = axes;
title(ax,'Kinetic energy of the optimal Cart Pendulum offset-correction maneuver')
hold(ax,"on"); grid(ax,"on")
xlabel('stage'); ylabel('Kinetic Energy')
plot(ax,[0,N],[1 1]*max_kinetic_energy,Color='k',LineStyle='--',LineWidth=1.1,DisplayName='Max. Kinetic Energy')
plot(ax,0:N,K,'Marker','+',Color=GetColorCode('y',0.9),MarkerSize=8,LineStyle='-',LineWidth=1.3,DisplayName='Kinetic Energy')
legend(ax)

%% However!, the TRYMPC class supports automatically plotting the general stage-constraint values across stages:

% Plot the "kinetic_energy" stage constraints that was defined earlier:
C.display_constraints("inequality","kinetic_energy");

%{
Note that this constraint value is not the kinetic energy itself, rather it
is the negative of the kinetic energy plus the "max_kinetic_energy". This
is because the constraint itself is defined that way in order to formulate
the bound on the form; h() >= 0.
%}
