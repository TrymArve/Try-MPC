%%% Testing simulator with cart-pendulum system

fresh


M_PC = trympcmodels("Pendulum_on_Cart");

model_number = 1;
parameters = M_PC.numeric_model(model_number).parameters.vec;
initial_state = M_PC.numeric_model(model_number).initial_state.vec;
tf = 10;

% uncontrolled:
u0 = zeros(M_PC.dim.input,1);




%% Dynamics function comparison

slow = [];
fast = [];
for i = 1:1000
   args = M_PC.example_args;
   tic; num_dyn      = M_PC.dynamics.call(args).out; slow(end+1) = toc;
   tic; num_dyn_fast = M_PC.dynamics_alt(args.state,args.input,args.param); fast(end+1) = toc;
end
disp(['------ number of tests: ', num2str(i)])
disp("dyn     : " + sec2str(max(slow)) +" "+ sec2str(mean(slow)) +" "+ sec2str(median(slow)) +" "+ sec2str(min(slow)))
disp("dyn_fast: " + sec2str(max(fast)) +" "+ sec2str(mean(fast)) +" "+ sec2str(median(fast)) +" "+ sec2str(min(fast)))



%% Simulate

high_tol = 1e-10;
low_rel_tol = 1e-5;
low_abs_tol = 1e-5;
low_tol = 1e-6;

% get a true solution:
opts = odeset(RelTol=high_tol,AbsTol=high_tol);
tic
[t_sim,s_sim] = ode45(@(t,s) full(M_PC.dynamics_alt(s,u0,parameters)), [0, tf], initial_state,opts);
high_time = toc;

plot(t_sim,s_sim,Color='g'); hold on
plot([0 10],[0 0],'k--',HandleVisibility='off')
plot([0 10],[1 1]*pi,'k--',HandleVisibility='off')
plot([0 10],[1 1]*pi*2,'k--',HandleVisibility='off')

% reduce the abs tolerance for speed, and compare:
opts = odeset(RelTol=high_tol,AbsTol=low_abs_tol);
tic
[t_sim,s_sim] = ode45(@(t,s) full(M_PC.dynamics_alt(s,u0,parameters)), [0, tf], initial_state,opts);
high_rel_time = toc;
plot(t_sim,s_sim,Color='b');

% reduce the rel tolerance for speed, and compare:
opts = odeset(RelTol=low_rel_tol,AbsTol=high_tol);
tic
[t_sim,s_sim] = ode45(@(t,s) full(M_PC.dynamics_alt(s,u0,parameters)), [0, tf], initial_state,opts);
high_abs_time = toc;
plot(t_sim,s_sim,Color='r');

% reduce the rel tolerance for speed, and compare:
opts = odeset(RelTol=low_tol,AbsTol=low_tol);
tic
[t_sim,s_sim] = ode45(@(t,s) full(M_PC.dynamics_alt(s,u0,parameters)), [0, tf], initial_state,opts);
low_time = toc;
plot(t_sim,s_sim,Color='m');

disp("High both: (green)" + sec2str(high_time))
disp("High rel : (blue) " + sec2str(high_rel_time))
disp("High abs : (red)  " + sec2str(high_abs_time))
disp("Low both : (magen)" + sec2str(low_time))


%% Discretizer (Integrator)

D_ERK4 = trympcDISCRETIZER("ERK",M_PC,"ERK4");


%% Propogate the discrete dynamics:

Dt = 0.005;
time = 0:Dt:tf;
state = initial_state;
args.Dt = Dt;
args.input = u0;
args.param = parameters;
tic
for i = time(2:end)
   args.state = state(:,end);
   state(:,end+1) = full(D_ERK4.Next.call(args).out);
end
ERK4_time = toc;

plot(time(1:size(state,2)),state,'k--')
disp("ERK4     : (black)" + sec2str(ERK4_time))

%% Create NLP (cart pendulum)
fresh
M_PC = trympcmodels("Pendulum_on_Cart");
num_mod = M_PC.numeric_model(1);

D_IRK = trympcDISCRETIZER("IRK",M_PC,"Crank-Nicolson (Implicit)","n_increments",2);

lb_state = structor;
lb_state.str.cart_position = -5;
lb_state.str.pendulum_angle = -pi*2;
lb_state.str.cart_speed = -5;
lb_state.str.pendulum_speed = -pi;

ub_state = structor;
ub_state.str.cart_pos = 5;
ub_state.str.pendulum_angle = pi*2;
ub_state.str.cart_speed = 5;
ub_state.str.pendulum_speed = pi;

input_bounds = [-1 1]*1;

% DOP   = trympcDOP("PC Multiple Shooting",D_IRK,...
%                   N_horizon=5,...
%                   quad_cost=["state","input"],...
%                   bounds_state=[lb_state.vec ub_state.vec],...
%                   bounds_input=input_bounds);

DOP   = trympcDOP("PC Multiple Shooting",D_IRK,...
                  N_horizon=100,...
                  quad_cost=["state","input"],bounds_input=input_bounds);


%% Solve
num_mod.parameters.str.pendulum_length = 2;
num_mod.initial_state.str.pendulum_angle = 0.01;

quad.state =  10*ones(4,1);
quad.input = 1;

T = 10;

[solution,solver] = DOP.solve(num_mod,quad=quad,T_horizon=T);
DOP.show(solution.decision);


%% Simulate result

nmpc = NMPC("Classic NMPC of Pendulum",DOP);

nmpc.numeric_model = num_mod;
nmpc.toggle_archive = true;
nmpc.quad = quad;
nmp.T_horizon = 10;


