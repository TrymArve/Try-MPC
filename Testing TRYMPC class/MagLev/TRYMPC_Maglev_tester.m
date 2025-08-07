fresh

%%% Initiate class with variable names:
C = TRYMPC('Tester Instance',...
    state = ["r1" "r2" "r3" "phi" "theta" "v1" "v2" "v3" "omega_x" "omega_y"], ...
    input = ["u1" "u2" "u3" "u4"],...
    param = "dummy_parameter");%, ...
    %param = ["m_lev" "M_lev" "I_lev" "g" "mu0" "r_base_array" "m_base" "r_base" "n_w"]);

params = MagLev_parameters;
    
%%% Define Dynamics:

state = @(s) [s.r1 s.r2 s.r3 s.phi s.theta s.v1 s.v2 s.v3 s.omega_x s.omega_y]';
input = @(i) [i.u1 i.u2 i.u3 i.u4]';

% dynamics (dq = f(q,u))
dynamics = @(s,a,i,p) f(state(s),input(i),params);

% Apply to Trympc:
C.def_dynamics(dynamics)

%% Configure:

C.plotting.display_names.state.r1      = "$x$";
C.plotting.display_names.state.r2      = "$y$";
C.plotting.display_names.state.r3      = "$z$";
C.plotting.display_names.state.phi     = "$\phi$";
C.plotting.display_names.state.theta   = "$\theta$";
C.plotting.display_names.state.v1      = "$\dot{x}$";
C.plotting.display_names.state.v2      = "$\dot{y}$";
C.plotting.display_names.state.v3      = "$\dot{z}$";
C.plotting.display_names.state.omega_x = "$\dot{\phi}$";
C.plotting.display_names.state.omega_y = "$\dot{\theta}$";
C.plotting.display_names.input.u1      = "$u_1$";
C.plotting.display_names.input.u2      = "$u_2$";
C.plotting.display_names.input.u3      = "$u_3$";
C.plotting.display_names.input.u4      = "$u_4$";

C.parameters.str.dummy_parameter = 0; % ignore this

%%% Find equilibrium point (Hans' code)
index = @(A,i) A(i);
fz = @(z) index(f([0,0,z,zeros(1,7)]',[0,0,0,0]',params),8);  % state is now 10x1
zeq =  fzero(fz,0.1);
xeq = [0,0,zeq,zeros(1,7)]';
ueq = [0,0,0,0]';

%%% Define initial state:
equilibrium = structor;
equilibrium.str.r1      = 0      + 0;
equilibrium.str.r2      = 0      + 0;
equilibrium.str.r3      = zeq    + 0;
equilibrium.str.phi     = 0      + 0;
equilibrium.str.theta   = 0      + 0;
equilibrium.str.v1      = 0      + 0;
equilibrium.str.v2      = 0      + 0;
equilibrium.str.v3      = 0      + 0;
equilibrium.str.omega_x = 0      + 0;
equilibrium.str.omega_y = 0      + 0;
C.initial_state = equilibrium.copy;

duration = 2;
%% Simulate

C.simulate(duration,simulator="ode15s")
C.display_simulation("state",["r1", "r2", "r3", "phi", "theta"],"input","all");

%% LQR Controller

% Linearize model
xlp = xeq;
ulp = ueq;
[A_lqr,B_lqr,~] = linearizeModel(@f, @h, xlp, ulp, params);

% Define LQR controller
Q = diag([ ...
    1e1,1e1,1e1, ...     % position
    1e1,1e1, ...         % orientation (phi, theta)
    1e1,1e1,1e5, ...     % velocity
    1e1,1e1              % angular velocity
]);
R = 1e2*eye(params.n_u);

% Create controller
K = lqr(ss(A_lqr, B_lqr, eye(10), []), Q, R);
u = @(C,t,x) -K*x;


%% Simulate with minor offset:

C.initial_state.str.r1      = 0      - 0.001;
C.initial_state.str.r2      = 0      + 0.001;
C.initial_state.str.r3      = zeq    + 0.02;
C.initial_state.str.phi     = 0      + pi/10;
C.initial_state.str.theta   = 0      + 0;
C.initial_state.str.v1      = 0      + 0;
C.initial_state.str.v2      = 0      + 0;
C.initial_state.str.v3      = 0      + 0;
C.initial_state.str.omega_x = 0      + 0;
C.initial_state.str.omega_y = 0      + 0;

C.simulate(duration,controller_type="custom",controller_custom=u,simulator="ode15s")
C.display_simulation("state",["r1", "r2", "r3", "phi", "theta"],"input","all");



%% Define optimization problem:

C.clear_problem

%%% Objective:
% C.def_objective(quadratic=["Q","R","dQ","dR","Q_terminal"]);
C.def_objective(quadratic=["Q","R","Q_terminal"]);

%%% Define Integrator:
C.def_integrator("ERK4","n_increments",1)
% C.def_integrator("Explicit Euler","n_increments",1)
% C.def_integrator("collocation", ...
%                  "n_increments",2, ...
%                  "collocation_polynomial_order",2, ...
%                  "collocation_polynomial_type","legendre")


%%% Constraints:
C.def_stage_constraints("lower_bounds",["u1" "u2" "u3" "u4"],"upper_bounds",["u1" "u2" "u3" "u4"])


%%% Horizon:
horizon_length = 10; % n. of points on horizon
C.def_horizon(horizon_length)

%%% Display problem properties:
C.display_problem(3)
% C.display_sparsity("equality")

%% Configure optimization problem:

%%% Quadratic Objective:

% %%% stabilierer innen for |u|<1
% Q_penalty = structor;
% Q_penalty.str.r1      = 10;
% Q_penalty.str.r2      = 10;
% Q_penalty.str.r3      = 10000;
% Q_penalty.str.phi     = 1;
% Q_penalty.str.theta   = 1;
% Q_penalty.str.v1      = 100;
% Q_penalty.str.v2      = 100;
% Q_penalty.str.v3      = 1000;
% Q_penalty.str.omega_x = 1;
% Q_penalty.str.omega_y = 1;
% C.quadratic_cost.Q = Q_penalty.vec;
% C.quadratic_cost.R = ones(1,4)*1000;
% C.quadratic_cost.Q_terminal = C.quadratic_cost.Q*100000;
% % initial state
% offset = structor;
% offset.str.r1      = 0.00006; % 0.001;
% offset.str.r2      = 0.0001; % 0.001;
% offset.str.r3      = 0.01; % 0.02;
% offset.str.phi     = 0.001; % pi/10;
% offset.str.theta   = 0;
% offset.str.v1      = 0;
% offset.str.v2      = 0;
% offset.str.v3      = 0;
% offset.str.omega_x = 0;
% offset.str.omega_y = 0;
% % change in offset:
% offset.str.r1      = 0.00011;
% offset.str.r2      = 0.0002;
% offset.str.r3      = 0.012; 
% offset.str.phi     = 0.0007; 


Q_penalty = structor;
Q_penalty.str.r1      = 10;
Q_penalty.str.r2      = 10;
Q_penalty.str.r3      = 10000;
Q_penalty.str.phi     = 1;
Q_penalty.str.theta   = 1;
Q_penalty.str.v1      = 100;
Q_penalty.str.v2      = 100;
Q_penalty.str.v3      = 1e3;
Q_penalty.str.omega_x = 1;
Q_penalty.str.omega_y = 1;
C.quadratic_cost.Q = Q_penalty.vec;
% dQ_penalty = structor;
% dQ_penalty.str.r1      = 1;
% dQ_penalty.str.r2      = 1;
% dQ_penalty.str.r3      = 0;
% dQ_penalty.str.phi     = 1;
% dQ_penalty.str.theta   = 1;
% dQ_penalty.str.v1      = 1;
% dQ_penalty.str.v2      = 1;
% dQ_penalty.str.v3      = 10000;
% dQ_penalty.str.omega_x = 1;
% dQ_penalty.str.omega_y = 1;
% C.quadratic_cost.dQ = dQ_penalty.vec;

C.quadratic_cost.R = ones(1,4)*1e4;
% C.quadratic_cost.dR = ones(1,4);

C.quadratic_cost.Q_terminal = C.quadratic_cost.Q*1e5;
% C.quadratic_cost.Q_terminal.r1 = 1e7;
% C.quadratic_cost.Q_terminal.r2 = 1e7;
% %%% state bounds
control_bound = 20;
% % Lower
% C.bounds.lower.r1      = -0.005;
% C.bounds.lower.r2      = -0.001;
% C.bounds.lower.r3      = 0;
% C.bounds.lower.phi     = -pi/4;
% C.bounds.lower.theta   = -pi/4;
% C.bounds.lower.v1      = -0.1;
% C.bounds.lower.v2      = -0.1;
% C.bounds.lower.v3      = -0.1;
% C.bounds.lower.omega_x = -100;
% C.bounds.lower.omega_y = -100;
C.bounds.lower.u1 = -control_bound;
C.bounds.lower.u2 = -control_bound;
C.bounds.lower.u3 = -control_bound;
C.bounds.lower.u4 = -control_bound;

% % Upper
% C.bounds.upper.r1      = 0.005;
% C.bounds.upper.r2      = 0.001;
% C.bounds.upper.r3      = zeq + 0.1;
% C.bounds.upper.phi     = pi/4;
% C.bounds.upper.theta   = pi/4;
% C.bounds.upper.v1      = 0.1;
% C.bounds.upper.v2      = 0.1;
% C.bounds.upper.v3      = 0.1;
% C.bounds.upper.omega_x = 100;
% C.bounds.upper.omega_y = 100;
C.bounds.upper.u1 = control_bound;
C.bounds.upper.u2 = control_bound;
C.bounds.upper.u3 = control_bound;
C.bounds.upper.u4 = control_bound;

% Reference:
C.set_ref("state",equilibrium.vec)

% initial state
offset = structor;
offset.str.r1      = 0.00006; % 0.001;
offset.str.r2      = 0.0001; % 0.001;
offset.str.r3      = 0.01; % 0.02;
offset.str.phi     = 0.001; % pi/10;
offset.str.theta   = 0;
offset.str.v1      = 0;
offset.str.v2      = 0;
offset.str.v3      = 0;
offset.str.omega_x = 0;
offset.str.omega_y = 0;
C.initial_state = equilibrium.retrieve(equilibrium.vec + offset.vec);

% C.initial_state = equilibrium.copy;
initial_guess = C.cas.problem.decision.zeros;
initial_guess.str.state = equilibrium.vec.*ones(1,C.horizon.N+1);

T = 0.5;
C.set_T(T);

%% SQP settings:
quadprog_options = optimoptions("quadprog","Algorithm","interior-point-convex","ConstraintTolerance",1e-8,"Display","none","OptimalityTolerance",1e-5,"StepTolerance",1e-5);
C.set_SQP_settings("linesearch_method","none","tolerance_equality",1e-4,"QP_solver","quadprog_IP","max_N_iterations",1,"quadprog_IP_options",quadprog_options);

%% Optimize
ipoptions = struct('print_time',0);
ipoptions.ipopt.print_level         = 0;
ipoptions.ipopt.tol                 = 1e-8;
ipoptions.ipopt.acceptable_tol      = 1e-6;
ipoptions.ipopt.compl_inf_tol       = 1e-4;
ipoptions.ipopt.constr_viol_tol     = 1e-4;
ipoptions.ipopt.dual_inf_tol        = 1e-4;

C.solve("ipopt","max_iterations",5000,initial_guess_primal=initial_guess.vec,ipopt_nlpsol_options=ipoptions);
C.display_optimization("state",["r1", "r2", "r3", "phi", "theta"],"input","all");
ind_ipopt_sol = length(C.archive.optimizations);

%% Save solution
Solution = C.opt.decision.vec;
initial_guess = C.opt.decision.retrieve(Solution);

%% slightly change initial state
offset.str.r1      = 0.00011;
offset.str.r2      = 0.0002;
offset.str.r3      = 0.012; 
offset.str.phi     = 0.0007; 
C.initial_state = equilibrium.retrieve(equilibrium.vec + offset.vec);

%% Solve with SQP

C.solve("sqp","max_iterations",1,initial_guess_primal=initial_guess.vec);
C.display_optimization("state",["r1", "r2", "r3", "phi", "theta"],"input","all");
ind_sqp_sol = length(C.archive.optimizations);



%% Simulate with open-loop optimal controller:

% Define OL optimal controller
u_OL = C.opt.decision.str.input;
u_OL = @(C,t,x) u_OL(:,find((t+C.horizon.Dt*0.1)<C.opt.time,1)-1);

% Halt condition:
% halt = @(s,i,p) s(3) < 0.005 || norm(i,inf) > 10 || norm(s(1:3),inf) > 0.1 || norm(s([C.names.ind.state.phi, C.names.ind.state.theta]),inf) > 0.5;
halt = @(s,i,p) s(3) < 0.005 || norm(s(1:2),inf) > 10e-3  || norm(i,inf) > control_bound;

C.simulate(T,C.initial_state.vec,"controller_type","custom","controller_custom",u_OL,simulator="ode15s",sampling_time=C.horizon.Dt,halting_condition=halt)

tiles = C.display_simulation("state",["r1", "r2", "r3", "phi", "theta"],"input","all",mark_samples=true);

%% Add prediction to plot:

C.display_optimization(tiles=tiles,optimization_number=length(C.archive.optimizations),linestyle="--",opacity=0.5);


%% Clear Archive:
C.archive
C.clear_archive



%% Closed-Loop NMPC simulation

C.simulate(2,C.initial_state.vec,sampling_time=0.005,controller_type="RTI",halting_condition=halt,initial_guess_primal=initial_guess.vec,simulator="ode15s")
sim_tiles = C.display_simulation("state",["r1", "r2", "r3", "phi", "theta"],"input","all",mark_samples=false);

% NMPC info:
sim_opts = C.archive.simulations{end}.optimizations;
opts = [sim_opts{:}];
sol_spec = [opts.solver_specific];
QP_sol = [sol_spec.QP_sol];QP_sol = [QP_sol{:}];
exit = [QP_sol.EXITFLAG];
niter = [opts.n_iterations];
soltimes = [opts(:).solve_time];
total = [soltimes(:).total];
disp(['    n. iter: min(',num2str(min(niter)),') | mean(',num2str(mean(niter)),') | max(',num2str(max(niter)),')'])
disp(['solve times: min(',num2str(min(total)),') | mean(',num2str(mean(total)),') | max(',num2str(max(total)),')'])
disp([' QPs converged: ',num2str(100*sum(exit)/length(exit)),'%'])
%% Add open-loop predictions
C.display_optimization(tiles=sim_tiles,optimizations=sim_opts,optimization_number=1:5:length(sim_opts),linestyle="-",multiplot="ontop",opacity=0.5,colors=["black"]);

%% 
while 1

   u = C.solve(x)


end