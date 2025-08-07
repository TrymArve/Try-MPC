%% Full test
fresh



%%%%%%%%%%%%%%%%% Initiate class with variable names:
C = TRYMPC('Tester Instance',...
    state = ["x" "th" "dx" "dth"], ...
    algeb = "z",...
    input = "ux", ...
    param = ["L" "g" "mx" "mth"]);




%%%%%%%%%%%%%%%%%  Define dynamics:
s = C.cas.state;
a = C.cas.algeb;
i = C.cas.input;
p = C.cas.param;
%{
 --- This is equivalent to doing this:
s.x = sym('x');
s.th = sym('th');
s.dx = sym('dx');
s.dth = sym('dth');
a.z = sym('z');
i.ux = sym('ux');
p.L = sym('L');
p.g = sym('g');
p.mx = sym('mx');
p.mth = sym('mth');
---------------------
%}

% Cart-Pendulum dynamics:
dynM = [       p.L^2*p.mth, -p.L*p.mth*cos(s.th);
     -p.L*p.mth*cos(s.th),       p.mx + p.mth];
dynC = [               0, 0;
     p.L*s.dth*p.mth*sin(s.th), 0];
dynG = [-p.L*p.g*p.mth*sin(s.th);
              0];
dynamics = [ s.dx;
             s.dth;
             dynM\(-dynC*[s.dx; s.dth] - dynG + [i.ux; 0])];
% ... and some algebraics (random example)
algebraics = a.z + s.x;

% Apply to Trympc:
C.def_dynamics(dynamics,algebraics)





%%%%%%%%%%%%%%%%%  Objective:

% Some weird non-quadratic stage cost is fine: (note the state variables represent deviation from reference, not absolute value. So this means: (x-ref) + (th-th_ref)^2*(dx-dx_ref))
stage_cost = s.x + s.th^2*s.dx;

% Define objective:
C.def_objective(stage_cost=stage_cost,quadratic=["Q","Z","R","dQ","dZ","dR"]);

% Then determine some quadratic cost weights (deviation from reference)
C.quadratic_cost.Q = [1 10 0.1 0.1];
C.quadratic_cost.Z = 2;
C.quadratic_cost.R = 1;

% Furthermore, one may add penalty to change in variables from one stage to
% the next:
C.quadratic_cost.dQ = [0 0 1 0];
C.quadratic_cost.dZ = 2;
C.quadratic_cost.dR = 1;


%%%%%%%%%%%%%%%%% Define Integrator:

%%% Some examples:

% C.def_integrator("IRK4","n_increments",4)
% C.def_integrator("ERK4","n_increments",4)
% C.def_integrator("Explicit Euler","n_increments",4)
% C.def_integrator("Implicit Euler","n_increments",4)
C.def_integrator("collocation", ...
                 "n_increments",5, ...
                 "collocation_polynomial_order",5, ...
                 "collocation_polynomial_type","legendre")

%%% Inspect result:
% present(C.integrator)



%%%%%%%%%%%%%%%%% Constraints:

% Define stage-wise constraints:
equality.fixed_th = s.th;
equality.constant_th = s.dth;
inequality.lower_bound_x = 4 - s.x;
inequality.upper_bound_x = s.x - 3;

% Apply:
C.def_stage_constraints(equality=equality,inequality=inequality)



%%%%%%%%%%%%%%%%% Horizon:

%{
In this context, "horizon" refers to the variables on each stage of the
discretized horizon. "Horizon length" means the number of stages
(discrete instances), or potentially the amount of time into the future.

Here we build all the horizon variables (state/algeb/input trajectories),
and simultaneously construct the dynamic constraints based on the
interator, the algebraic constraints and general stage-wise constraints
for each stage on the horizon.
Furthermore, all objective terms (all terms for all stages) over the horizon are created both
individually and as a sum across the horizon.

By defining the horizon, you also automatically call the internal
"def_problem" function, which constructs the final problem formulation.
Note that this problem formulation is somewhat distinc from the horizon and
its state constraints and costs. This is because, dependent on the
optimization algorithm (solver), the formulation might be different
depending on the interface. Therefore the problem fomrmulation is a specifi
formulation that is used by the internal solvers. Other formulations will
automatically be build for whatever solve you try to use.
%}


horizon_length = 50;
Dt = 1;
C.def_horizon(horizon_length,Dt=Dt)
C.display_problem(4)
C.display_sparsity(["objective","inequality","equality"])


%% Test without inputs not algebraic states

fresh



%%%%%%%%%%%%%%%%% Initiate class with variable names:
C = TRYMPC('Tester Instance',...
    state = ["x" "th" "dx" "dth"], ...
    param = ["L" "g" "mx" "mth"]);




%%%%%%%%%%%%%%%%%  Define dynamics:
s = C.cas.state;
p = C.cas.param;

% Cart-Pendulum dynamics:
dynM = [       p.L^2*p.mth, -p.L*p.mth*cos(s.th);
     -p.L*p.mth*cos(s.th),       p.mx + p.mth];
dynC = [               0, 0;
     p.L*s.dth*p.mth*sin(s.th), 0];
dynG = [-p.L*p.g*p.mth*sin(s.th);
              0];
dynamics = [ s.dx;
             s.dth;
             dynM\(-dynC*[s.dx; s.dth] - dynG + [0; 0])];

% Apply to Trympc:
C.def_dynamics(dynamics)





%%%%%%%%%%%%%%%%%  Objective:

% Some weird non-quadratic stage cost is fine: (note the state variables represent deviation from reference, not absolute value. So this means: (x-ref) + (th-th_ref)^2*(dx-dx_ref))
stage_cost = s.x + s.th^2*s.dx;

% Define objective:
C.def_objective(stage_cost=stage_cost,quadratic=["Q","dQ"]);

% Then determine some quadratic cost weights (deviation from reference)
C.quadratic_cost.Q = [1 10 0.1 0.1];

% Furthermore, one may add penalty to change in variables from one stage to
% the next:
C.quadratic_cost.dQ = [0 0 1 0];


%%%%%%%%%%%%%%%%% Define Integrator:

%%% Some examples:

% C.def_integrator("IRK4","n_increments",4)
% C.def_integrator("ERK4","n_increments",4)
C.def_integrator("Explicit Euler","n_increments",4)
% C.def_integrator("Implicit Euler","n_increments",4)
% C.def_integrator("collocation", ...
%                  "n_increments",5, ...
%                  "collocation_polynomial_order",5, ...
%                  "collocation_polynomial_type","legendre")

%%% Inspect result:
% present(C.integrator)



%%%%%%%%%%%%%%%%% Constraints:

% % Define stage-wise constraints:
% equality.fixed_th = s.th;
% equality.constant_th = s.dth;
% inequality.lower_bound_x = 4 - s.x;
% inequality.upper_bound_x = s.x - 3;
% 
% % Apply:
% C.def_stage_constraints(equality,inequality)



%%%%%%%%%%%%%%%%% Horizon:

%{
In this context, "horizon" refers to the variables on each stage of the
discretized horizon. "Horizon length" means the number of stages
(discrete instances), or potentially the amount of time into the future.

Here we build all the horizon variables (state/algeb/input trajectories),
and simultaneously construct the dynamic constraints based on the
interator, the algebraic constraints and general stage-wise constraints
for each stage on the horizon.
Furthermore, all objective terms (all terms for all stages) over the horizon are created both
individually and as a sum across the horizon.

By defining the horizon, you also automatically call the internal
"def_problem" function, which constructs the final problem formulation.
Note that this problem formulation is somewhat distinc from the horizon and
its state constraints and costs. This is because, dependent on the
optimization algorithm (solver), the formulation might be different
depending on the interface. Therefore the problem fomrmulation is a specifi
formulation that is used by the internal solvers. Other formulations will
automatically be build for whatever solve you try to use.
%}


horizon_length = 500;
Dt = 1;
C.def_horizon(horizon_length,Dt=Dt)
C.display_problem(4)
C.display_sparsity("equality")




%% Genuine Test for cart pendulum
fresh


%%%%%%%%%%%%%%%%% Initiate class with variable names:
C = TRYMPC('Tester Instance',...
    state = ["x" "th" "dx" "dth"], ...
    input = "ux", ...
    param = ["L" "g" "mx" "mth"]);

%%%%%%%%%%%%%%%%%  Define dynamics:
s = C.cas.state;
% a = C.cas.algeb;
i = C.cas.input;
p = C.cas.param;

u_pendulum = 0;

% % cart acceleration
% ddx = -(p.L*i.ux + u_pendulum*cos(s.th) + p.L^2*s.dth^2*p.mth*sin(s.th) - p.L*p.g*p.mth*cos(s.th)*sin(s.th))/(p.L*(- p.mth*cos(s.th)^2 + p.mx + p.mth));
% 
% % pendulum angular acceleration
% ddth = -(p.mx*u_pendulum+ p.mth*u_pendulum + p.L*p.mth*i.ux*cos(s.th) - p.L*p.g*p.mth^2*sin(s.th) + p.L^2*s.dth^2*p.mth^2*cos(s.th)*sin(s.th) - p.L*p.g*p.mx*p.mth*sin(s.th))/(p.L^2*p.mth*(- p.mth*cos(s.th)^2 + p.mx + p.mth));
% 
% % dynamics (dq = f(q,u))
% dynamics = [s.dx;
%             s.dth;
%          	ddx;
%             ddth];

% cart acceleration
ddx = @(s,a,i,p) -(p.L*i.ux + u_pendulum*cos(s.th) + p.L^2*s.dth^2*p.mth*sin(s.th) - p.L*p.g*p.mth*cos(s.th)*sin(s.th))/(p.L*(- p.mth*cos(s.th)^2 + p.mx + p.mth));

% pendulum angular acceleration
ddth = @(s,a,i,p) -(p.mx*u_pendulum+ p.mth*u_pendulum + p.L*p.mth*i.ux*cos(s.th) - p.L*p.g*p.mth^2*sin(s.th) + p.L^2*s.dth^2*p.mth^2*cos(s.th)*sin(s.th) - p.L*p.g*p.mx*p.mth*sin(s.th))/(p.L^2*p.mth*(- p.mth*cos(s.th)^2 + p.mx + p.mth));

% dynamics (dq = f(q,u))
dynamics = @(s,a,i,p) [s.dx;
                       s.dth;
         	           ddx(s,a,i,p);
                       ddth(s,a,i,p)];


% Apply to Trympc:
C.def_dynamics(dynamics)

%%%%%%%%%%%%%%% Define LaTex Display names:
C.plotting.display_names.state.x = "$x$";
C.plotting.display_names.state.th = "$\theta$";
C.plotting.display_names.state.dx = "$\dot{x}$";
C.plotting.display_names.state.dth = "$\dot{\theta}$";
C.plotting.display_names.input.ux = "$u_{x}$";

%% Simulate (PD - controller)
close all
C.clear_archive

C.parameters.str.L = 1;
C.parameters.str.g = 9.81;
C.parameters.str.mx = 5;
C.parameters.str.mth = 3;

init_state = structor;
init_state.str.x   = 0.1;
init_state.str.th  = 0.2;
init_state.str.dx  = 0;
init_state.str.dth = 0;

% PD controller: L = 1
x_kp =  -1;
x_kd =  -1;
th_kp = 100;
th_kd = 10;

controller = @(~,~,x) x(1)*x_kp + x(2)*th_kp + x(3)*x_kd + x(4)*th_kd;

disp('------------ SIMULATING: cart-pendulum - PD controlled')
C.simulate(50,init_state.vec,...   % Basic settings (required)
   sampling_time=inf,...           % Choose control rate (optinal)
   simulator="ode15s",...          % Choose simulator (optinal)
   ode_options = odeset('RelTol',10^(-7),'AbsTol',10^(-7)),... % choose ode settings (optinal)
   controller_type="custom",...    % Choose controller type
   controller_custom=controller,...% "custom" specific: apply controller handle
   control_delay_scaler=0)         % scale control delay (set to 0 to turn off control delay)


C.display_simulation;


%% Simulate (PD - controller with zero-order hold (piecewise constant with given sampling time))
close all
C.clear_archive

C.parameters.str.L = 5;
C.parameters.str.g = 9.81;
C.parameters.str.mx = 5;
C.parameters.str.mth = 3;

init_state = structor;
init_state.str.x   = 0.01;
init_state.str.th  = 0.0;
init_state.str.dx  = 0;
init_state.str.dth = 0;

% PD controller: L = 5
x_kp =  -50;
x_kd =  -10;
th_kp = 1500;
th_kd = 150;

controller = @(C,t,x) PD1(C,t,x,x_kp,x_kd,th_kp,th_kd);

disp('------------ SIMULATING: cart-pendulum - PD controlled w/ sampling time')
C.simulate(10,init_state.vec,...   % Basic settings (required)
   sampling_time=0.1,...           % Choose control rate (optinal)
   simulator="ode15s",...          % Choose simulator (optinal)
   ode_options = odeset('RelTol',10^(-7),'AbsTol',10^(-7)),... % choose ode settings (optinal)
   controller_type="custom",...    % Choose controller type
   controller_custom=controller,...% "custom" specific: apply controller handle
   control_delay_scaler=0)         % scale control delay (set to 0 to turn off control delay)

C.display_simulation;


%% Simulate (PD - controller with zero-order hold and control delay)
close all
C.clear_archive

C.parameters.str.L = 5;
C.parameters.str.g = 9.81;
C.parameters.str.mx = 5;
C.parameters.str.mth = 3;

init_state = structor;
init_state.str.x   = 0.01;
init_state.str.th  = pi+0.2;
init_state.str.dx  = 0;
init_state.str.dth = 0;

% PD controller: L = 5
x_kp =  0;
x_kd =  0;
th_kp = 100;
th_kd = -1000;

controller = @(C,t,x) PD2(C,t,x,x_kp,x_kd,th_kp,th_kd);

disp('------------ SIMULATING: cart-pendulum - PD controlled w/ sampling time, w/ control delay')
C.simulate(1,init_state.vec,...   % Basic settings (required)
   sampling_time=0.2,...           % Choose control rate (optinal)
   simulator="ode15s",...          % Choose simulator (optinal)
   ode_options = odeset('RelTol',10^(-7),'AbsTol',10^(-7)),... % choose ode settings (optinal)
   controller_type="custom",...    % Choose controller type
   controller_custom=controller,...% "custom" specific: apply controller handle
   control_delay_scaler=1000)         % scale control delay (set to 0 to turn off control delay)

C.display_simulation;



%% Simulate (PD - controller, comparisons)
close all
C.clear_archive

C.parameters.str.L = 1;
C.parameters.str.g = 9.81;
C.parameters.str.mx = 5;
C.parameters.str.mth = 3;

init_state = structor;
init_state.str.x   = 0.1;
init_state.str.th  = 0.2;
init_state.str.dx  = 0;
init_state.str.dth = 0;

% PD controller:
x_kp =  -1;
x_kd =  -1;
th_kp = 100;
th_kd = 10;

controller = @(C,t,x) PD1(C,t,x,x_kp,x_kd,th_kp,th_kd);

disp('------------ SIMULATING: cart-pendulum - PD controlled, comparisons')
disp('---- Sim 1...')
C.simulate(50,init_state.vec,...   % Basic settings (required)
   sampling_time=inf,...           % Choose control rate (optinal)
   simulator="ode15s",...          % Choose simulator (optinal)
   ode_options = odeset('RelTol',10^(-7),'AbsTol',10^(-7)),... % choose ode settings (optinal)
   controller_type="custom",...    % Choose controller type
   controller_custom=controller,...% "custom" specific: apply controller handle
   control_delay_scaler=0)         % scale control delay (set to 0 to turn off control delay)
C.display_simulation("simulation_number",1);
%% zero-order hold is significant
disp('---- Sim 2...')
C.simulate(50,init_state.vec,...   % Basic settings (required)
   sampling_time=0.1,...           % Choose control rate (optinal)
   simulator="ode15s",...          % Choose simulator (optinal)
   ode_options = odeset('RelTol',10^(-7),'AbsTol',10^(-7)),... % choose ode settings (optinal)
   controller_type="custom",...    % Choose controller type
   controller_custom=controller,...% "custom" specific: apply controller handle
   control_delay_scaler=0)         % scale control delay (set to 0 to turn off control delay)
C.display_simulation("simulation_number",2,"mark_samples",false);
%% Control delay is insignificant for a PD
disp('---- Sim 3...')
C.simulate(50,init_state.vec,...   % Basic settings (required)
   sampling_time=0.1,...           % Choose control rate (optinal)
   simulator="ode15s",...          % Choose simulator (optinal)
   ode_options = odeset('RelTol',10^(-7),'AbsTol',10^(-7)),... % choose ode settings (optinal)
   controller_type="custom",...    % Choose controller type
   controller_custom=controller,...% "custom" specific: apply controller handle
   control_delay_scaler=1500)         % scale control delay (set to 0 to turn off control delay)
C.display_simulation("simulation_number",3,"mark_samples",false);




%% Define NMPC problem

C.clear_problem

%%%%%%%%%%%%%%%%%  Objective:
% Define objective:
C.def_objective(quadratic=["Q","R","dR","Q_terminal"]);


%%%%%%%%%%%%%%%%% Define Integrator:
C.def_integrator("Explicit Euler","n_increments",1)
% C.def_integrator("Implicit Euler","n_increments",1)
% C.def_integrator("ERK4","n_increments",1);
% C.def_integrator("collocation","collocation_polynomial_order",4,  "n_increments",3,   collocation_polynomial_type="legendre")

%%%%%%%%%%%%%%%%% Constraints:
% % Stage
% inequality.sum_of_pos = s.x + s.th + 10;
% C.def_stage_constraints("lower_bounds",["x","ux"],"upper_bounds","ux",inequality=inequality)
C.def_stage_constraints("lower_bounds","ux","upper_bounds","ux")

% % Terminal
% terminal_equality.cart_at_zero = s.x;
% terminal_equality.pendulum_at_zero = s.th;
% terminal_equality.pendulum_at_still = s.dth;
% C.def_terminal_constraint("equality",terminal_equality,"upper_bounds","dx","lower_bounds","dx")

%%%%%%%%%%%%%%%%% horizon:
horizon_length = 150;
C.def_horizon(horizon_length,"primaldual","stacked")

C.display_problem(4)
C.display_sparsity("KKT")

%% Set parameters and initial state
C.parameters.str.L = 1;
C.parameters.str.g = 9.81;
C.parameters.str.mx = 5;
C.parameters.str.mth = 3;

init_state = structor;
init_state.str.x   = 1;
init_state.str.th  = 0;
init_state.str.dx  = 0;
init_state.str.dth = 0;

C.quadratic_cost.Q = [1 10 0.1 0.1];
C.quadratic_cost.R = 1;
C.quadratic_cost.dR = 1;
C.quadratic_cost.Q_terminal = C.quadratic_cost.Q*350;

% C.bounds.lower.x = -10;
C.bounds.lower.ux = -1;
C.bounds.upper.ux = 0.4;
% 
% C.terminal_bounds.lower.dx = -0.01;
% C.terminal_bounds.upper.dx = 0.01;

C.set_T(15);
%% Open-Loop optimization (ipopt)
   ipopt_nlpsol_options.ipopt.print_level = 0;
   ipopt_nlpsol_options.ipopt.max_iter = 100;
   ipopt_nlpsol_options.ipopt.tol = 1e-9;
   ipopt_nlpsol_options.ipopt.acceptable_tol = 1e-6;
   ipopt_nlpsol_options.ipopt.compl_inf_tol = 1e-9;
   ipopt_nlpsol_options.ipopt.constr_viol_tol = 1e-9;
   ipopt_nlpsol_options.ipopt.dual_inf_tol = 1e-4;
sol = C.solve("ipopt","initial_state",init_state.vec,"ipopt_nlpsol_options",ipopt_nlpsol_options);
index_ipopt = length(C.archive.optimizations);
C.display_optimization;

%% Open-Loop optimization (SQP)
C.set_SQP_settings("max_N_iterations",20,"tolerance_lagrangian",5,"tolerance_equality",10^-6)
sol = C.solve("sqp","initial_state",init_state.vec);
index_sqp = length(C.archive.optimizations);
C.display_optimization;

%% Plot together (sqp vs. ipopt)
C.display_optimization("optimization_number",[index_sqp index_ipopt])
%% Plot together (manual)
C.display_optimization("optimization_number",[1 2])

%% NMPC controlled (ipopt)

disp('------------ SIMULATING: cart-pendulum - NMPC (ipopt)')
C.simulate(5,init_state.vec,...   % Basic settings (required) (duration, initial-state)
   sampling_time=0.05,...             % Choose control rate (optinal)
   simulator="ode15s",...          % Choose simulator (optinal)
   ode_options = odeset('RelTol',10^(-7),'AbsTol',10^(-7)),... % choose ode settings (optional)
   controller_type="NMPC_ipopt",...    % Choose controller type
   control_delay_scaler=0,...   % scale control delay (set to 0 to turn off control delay)
   speed=false)
C.display_simulation(mark_samples=false);



%% testing convergence...

QPIPoptions = optimoptions("quadprog","Algorithm","interior-point-convex","Display","none",...
                           MaxIterations=2000,...
                           ConstraintTolerance = 1e-15,...
                           OptimalityTolerance=1e-15,...
                           StepTolerance=1e-15);
C.set_SQP_settings("quadprog_IP_options",QPIPoptions,...
                   "backtracking_rate", 0.5,...
                   "backtracking_min_step_size",1e-20,...
                   "tolerance_lagrangian", 1,...
                   "tolerance_equality",   4.5e-07,...
                   "tolerance_inequality", 1e-5)
C.def_integrator("ERK4","n_increments",3)
C.def_horizon(30)
C.solve("sqp");
C.archive.optimizations{end}.solver_specific.step_size;
% C.display_optimization;
 




%% controllers

function u = PD1(C,~,x,x_kp,x_kd,th_kp,th_kd)
tic;
u = x(1)*x_kp + x(2)*th_kp + x(3)*x_kd + x(4)*th_kd;
cd = toc;

if C.simulation.control_delay_scaler
   C.set_control_delay(cd);
   disp(['PD controller called: u = ',num2str(u)])
end
end

function u = PD2(C,~,x,x_kp,x_kd,th_kp,th_kd)
tic;
u = x(1)*x_kp + (x(2)-pi)*th_kp + x(3)*x_kd + x(4)*th_kd;
cd = toc;

C.set_control_delay(cd);

disp(['PD controller called: u = ',num2str(u),' - control delay: ',sec2str(cd)])
end