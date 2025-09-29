%% Create System:

fresh
Create_TSEx_Qm % Define the "TSEx" polymer CSTR reactor model

%% Parameter Values (original):
C.parameters.str.rhoC_p      = 360;
C.parameters.str.rhoC_pc     = 966.3;
C.parameters.str.f_eff       = 0.6;
C.parameters.str.hA          = 70;
C.parameters.str.I_f         = 0.5888;      % Don't change
C.parameters.str.A_d         = 5.59*10^(13);
C.parameters.str.A_p         = 1.06*10^(7); % Don't change
C.parameters.str.A_t         = 1.25*10^(9); % Don't change
C.parameters.str.E_d         = 14897;       % Don't change
C.parameters.str.E_p         = 3557;
C.parameters.str.E_t         = 843;
C.parameters.str.M_f         = 8.6981;
%C.parameters.str.Q_m         = 0.105;
C.parameters.str.Q_i         = 0.040001583018170; %0.03; % added this
C.parameters.str.T_cf        = 295;
C.parameters.str.T_f         = 330;
C.parameters.str.V           = 1000; % original: 3000
C.parameters.str.V_c         = 1000; % original: 3312.4
C.parameters.str.nDHr        = 16700;      % Don't change
C.parameters.str.M_m         = 104.14;
C.parameters.str.gamma       = 0.2;
C.parameters.str.lambda      = 0.6;

%% Parameter Values (stiff):
C.parameters.str.rhoC_p      = 360;
C.parameters.str.rhoC_pc     = 966.3;
C.parameters.str.f_eff       = 0.6;
C.parameters.str.hA          = 70;
C.parameters.str.I_f         = 0.5888;      % Don't change
C.parameters.str.A_d         = 5.59*10^13; %8*10^16;        % OLD 5.59*10^(13) % Crucial to make stiff !
C.parameters.str.A_p         = 1.06*10^(7); % Don't change
C.parameters.str.A_t         = 1.25*10^(9); % Don't change
C.parameters.str.E_d         = 14897;       % Don't change
C.parameters.str.E_p         = 4500;        % OLD 3557
C.parameters.str.E_t         = 843;
C.parameters.str.M_f         = 8.6981;
%C.parameters.str.Q_m         = 0.105;
C.parameters.str.Q_i         = 0.040001583018170; %0.03; % added this
C.parameters.str.T_cf        = 295;
C.parameters.str.T_f         = 330;
C.parameters.str.V           = 1000;              % OLD 3000
C.parameters.str.V_c         = 500;              % OLD 3312.4
C.parameters.str.nDHr        = 16700;      % Don't change
C.parameters.str.M_m         = 104.14;
C.parameters.str.gamma       = 0.2;
C.parameters.str.lambda      = 0.6;


%% Parameter Error: (difference from optimization model to simulation model)
param_error = C.parameters.zeros;
param_error.str.rhoC_p      = 10; % Connected with Q_m
param_error.str.rhoC_pc     = 10; % Connected with Q_m
% param_error.str.f_eff       = 0.01;
param_error.str.hA          = 1; % Connected with Q_m
% param_error.str.I_f         = 0.0001;      % Don't change
% param_error.str.A_d         = 1*10^(16);
% param_error.str.A_p         = 1.0*10^(10); % Don't change
% param_error.str.A_t         = 1.0*10^(12); % Don't change
% param_error.str.E_d         = 10;       % Don't change
% param_error.str.E_p         = 1;
% param_error.str.E_t         = 1;
% param_error.str.M_f         = 0.0001; % Connected with Q_m
% param_error.str.Q_m         =  0.001;
% param_error.str.Q_i         = 0.005; % added this
% param_error.str.T_cf        = 1;
% param_error.str.T_f         = 1;
param_error.str.V           = 3; % Connected with Q_m
param_error.str.V_c         = 3.4; % Connected with Q_m
% param_error.str.nDHr        = 10;      % Don't change
% param_error.str.M_m         = 0.14;
% param_error.str.gamma       = 0.0;
% param_error.str.lambda      = 0.0;

param_error_scaler = 1e2; %use -5 for super stiff system
%% Alternative 2: Volume
param_error = C.parameters.zeros;
param_error.str.rhoC_p      = -100; % Connected with Q_m
param_error.str.rhoC_pc     = -200; % Connected with Q_m
param_error.str.hA          = -20; % Connected with Q_m
param_error.str.V           = -2000; % Connected with Q_m
param_error.str.V_c         =  0; % Connected with Q_m

param_error_scaler = 0e-3; %use -5 for super stiff system
start_varying = 1/10; % fraction of duration at which to start the parameter varying

%% Duration 300 hours
duration = hours(300);
%% Duration 150 hours
duration = hours(150);
%% Duration 50 hours
duration = hours(50);
%% Duration 20 hours
duration = hours(20);
%% Duration 1 hour
duration = hours(1);

%% define randomness:
deviation = 1;
randy = @() 2*(rand-0.5)*deviation + 1;

%% Varying Parameters: step
varying_parameters = @(t) C.parameters.vec + randy()*((t>start_varying*duration) + 2*(t>(start_varying + (1-start_varying)/2)*duration))*param_error.vec*param_error_scaler; % Make a varying parameter function handle

%% Varying Parameters: gradual (linear)
varying_parameters = @(t) C.parameters.vec + randy()*((max(t,start_varying*duration)-start_varying*duration)/(duration*(1-start_varying)))*param_error.vec*param_error_scaler; % Make a varying parameter function handle

%% Varying Parameters: gradual (sqrt)
varying_parameters = @(t) C.parameters.vec + randy()*sqrt((max(t,start_varying*duration)-start_varying*duration)/(duration*(1-start_varying)))*param_error.vec*param_error_scaler; % Make a varying parameter function handle



%% Steady States

% (stable) 
stable_eq = structor;
stable_eq.str.x_I    = 0.000678299473481    *100;
stable_eq.str.x_M    = 0.024868324646906    *100;
stable_eq.str.x_T    =  3.309412492542791    *100;
stable_eq.str.x_Tc   = 3.179149361140690    *100;
stable_eq.str.x_D0   =   0.000005624677816    *100;
stable_eq.str.x_D1   =  0.187149596787991    *100;

stable_input = structor;
stable_input.str.Q_c    = 0.041180244769836;
stable_input.str.Q_m    = 0.104999313835524;

% (unstable)
unstable_eq = structor;
unstable_eq.str.x_I  =   0.063289699425013;
unstable_eq.str.x_M  =   2.121424920685393;
unstable_eq.str.x_T  = 345.0000337734414;
unstable_eq.str.x_Tc = 326.8782732745260;
unstable_eq.str.x_D0 =   0.003285841672448;
unstable_eq.str.x_D1 =  56.770090728413585;

unstable_input = structor;
unstable_input.str.Q_c    = 0.041180244769836;
unstable_input.str.Q_m    = 0.104999313835524;

%% Alternative Stable Steady State
load('stable_eq_saved.mat')
load('unstable_eq_saved.mat')
% (stable) 
stable_eq = stable_eq_saved(1);
stable_input = stable_eq_saved(2);

% (unstable) 
unstable_eq = unstable_eq_saved(1);
unstable_input = unstable_eq_saved(2);

%% Reference (stable)
C.set_ref("state",stable_eq.vec,"input",stable_input.vec)
%% Reference (unstable)
C.set_ref("state",unstable_eq.vec,"input",unstable_input.vec)
%% Stable Equilibrium as initial condition
initial_state = stable_eq;
steady_controller = stable_input;
%% Unstable Equilibrium as initial condition
initial_state = unstable_eq;
steady_controller = unstable_input;
%% Modify initial state:
initial_state_scaler = 1;
initial_state.str.x_I    = initial_state_scaler*initial_state.str.x_I;
initial_state.str.x_M    = initial_state_scaler*initial_state.str.x_M;
initial_state.str.x_T    = initial_state_scaler*initial_state.str.x_T;
initial_state.str.x_Tc   = initial_state_scaler*initial_state.str.x_Tc;
initial_state.str.x_D0   = initial_state_scaler*initial_state.str.x_D0;
initial_state.str.x_D1   = initial_state_scaler*initial_state.str.x_D1;

%% Simulate (uncontrolled)

C.simulate(duration,initial_state.vec,simulator="ode15s",controller_type="constant",controller_constant=steady_controller.vec,simulation_parameters=varying_parameters)
C.display_simulation("time_order","hours",reference=["state","input"]);
% C.archive.simulations{end}.sim.state(:,end) % use to extract steady state for new parameters.

%% In order to define stable steady state exactly
sim_end = C.archive.simulations{end}.sim.state(:,end);
sim_inputs = C.archive.simulations{end}.sim.input_effective(:,end);
% (stable) 
stable_eq_saved = [structor; structor];
stable_eq_saved(1).str.x_I    = sim_end(1);
stable_eq_saved(1).str.x_M    = sim_end(2);
stable_eq_saved(1).str.x_T    = sim_end(3);
stable_eq_saved(1).str.x_Tc   = sim_end(4);
stable_eq_saved(1).str.x_D0   = sim_end(5);
stable_eq_saved(1).str.x_D1   = sim_end(6);
stable_eq_saved(2).str.Q_c   = sim_inputs(1);
stable_eq_saved(2).str.Q_m   = sim_inputs(2);
% save('stable_eq_saved',"stable_eq_saved")



%% Get sensitivity matrix and extract element

 

%% Create Optimization Problem
N = 50;
C.clear_problem
C.def_objective("quadratic",["Q","R","Q_terminal"]);
C.def_integrator("ERK4","n_increments",6);
C.def_stage_constraints("lower_bounds",["x_I","x_M","x_T","x_Tc","x_D0","x_D1" "Q_c","Q_m"]);
C.def_horizon(N);

%% Configure:

%%% Objective:

% state:
Q = structor;
Q.str.x_I    = 1;
Q.str.x_M    = 0; 
Q.str.x_T    = 0;
Q.str.x_Tc   = 0;
Q.str.x_D0   = 1;
Q.str.x_D1   = 1;
C.quadratic_cost.Q = Q.vec;
C.quadratic_cost.Q_terminal = Q.vec*100;

% input:
R = structor;
R.str.Q_c    = 1;
R.str.Q_m    = 1;
C.quadratic_cost.R = R.vec*1e-2;

%%% Bounds:
C.bounds.lower.x_I    = 0;
C.bounds.lower.x_M    = 0; 
C.bounds.lower.x_T    = 0;
C.bounds.lower.x_Tc   = 0;
C.bounds.lower.x_D0   = 0;
C.bounds.lower.x_D1   = 0;
C.bounds.lower.Q_c    = 0;
C.bounds.lower.Q_m    = 0;

% Horizon Length:
C.set_T(hours(30));


%% Stable Equilibrium as initial condition
initial_state = stable_eq;
%% Unstable Equilibrium as initial condition
initial_state = unstable_eq;

%% Generate Initial Guess Based on Last Simulation:
initial_guess.primal = C.sim2guess; % interpolates an initial guess based on the last simulation

%% Open-Loop
C.solve("ipopt","display_result",true,"initial_state",initial_state.vec,initial_guess_primal=initial_guess.primal.vec);
C.display_optimization(time_order="hours",reference=["state","input"]);

%% Store as initial guess:
initial_guess = C.opt2guess;
% initial_dual_eq = C.archive.optimizations{end}.dual_eq; % cannot do this if the constraints are different
% initial_dual_in = C.archive.optimizations{end}.dual_in;

%% Open-Loop (with good initial guess)
C.solve("ipopt","display_result",true,"initial_state",initial_state.vec,initial_guess_primal=initial_guess.primal.vec,initial_guess_dual_eq=initial_guess.dual_eq,initial_guess_dual_in=initial_guess.dual_in);
C.display_optimization(time_order="hours",reference=["state","input"]);

%% Open-Loop (with good initial guess) SQP-algoritm
C.set_SQP_settings("tolerance_lagrangian",100,"tolerance_equality",1e-10,"backtracking_min_step_size",1e-4)
C.solve("sqp","display_result",true,"initial_state",initial_state.vec,initial_guess_primal=initial_guess.primal.vec,initial_guess_dual_eq=initial_guess.dual_eq,initial_guess_dual_in=initial_guess.dual_in);
C.display_optimization(time_order="hours",reference=["state","input"]);

%% SIMULATION

%% Stable Equilibrium as initial condition
initial_state = stable_eq;

%% Simulate NMPC (ipopt)
dur = hours(60);
samp_time = minutes(60);
C.simulate(dur,initial_state.vec,"controller_type","NMPC_ipopt","sampling_time",samp_time,initial_guess_primal=initial_guess.primal.vec,initial_guess_dual_eq=initial_guess.dual_eq,initial_guess_dual_in=initial_guess.dual_in,simulation_parameters=varying_parameters,use_sim_param_for_NMPC=true);
sim_ipopt = length(C.archive.simulations);
sim_number = sim_ipopt;
%% Simulate NMPC (SQP)
C.set_SQP_settings("max_N_iterations",1000,"tolerance_lagrangian",1e1)

dur = hours(20);
samp_time = minutes(60);
C.simulate(dur,initial_state.vec,"controller_type","NMPC_sqp","sampling_time",samp_time,initial_guess_primal=initial_guess.primal.vec,initial_guess_dual_eq=initial_guess.dual_eq,initial_guess_dual_in=initial_guess.dual_in,simulation_parameters=varying_parameters);
sim_RTI = length(C.archive.simulations);
sim_number = sim_RTI;

%% Plot
tiles = C.display_simulation(simulation_number=sim_number,time_order="hours",reference=["state","input"]);

% And stable reference
stable = structor;
stable.str.state = stable_eq.vec.*[1 1];
stable.str.input = stable_input.vec*[1 1];
time = C.archive.simulations{sim_number}.start_time + [0 C.archive.simulations{sim_number}.duration];
C.display_trajectory(stable,time,tiles=tiles,colors="black",linewidth=1.5,linestyle=":",time_order="hours");
%% Plot predictions along the way:
N_optimizations = length(C.archive.simulations{sim_number}.optimizations);
optnum_to_disp = round(linspace(1,N_optimizations,5));
C.display_optimization(tiles=tiles,optimizations=C.archive.simulations{sim_number}.optimizations,optimization_number=optnum_to_disp,time_order="hours",linestyle="--");

%% Extract end state
end_state_NMPC = C.archive.simulations{sim_number}.sim.state(:,end); % use to extract steady state for new parameters.
end_input_NMPC = C.archive.simulations{end}.optimizations{end}.decision.str.input(:,end); % use to extract steady input for new parameters.
%% End state as initial condition
initial_state = initial_state.retrieve(end_state_NMPC);

%% Simulate (uncontrolled)
C.simulate(duration,end_state_NMPC,simulator="ode15s",controller_type="constant",controller_constant=end_input_NMPC,simulation_parameters=varying_parameters)
C.display_simulation("time_order","hours",reference=["state","input"]);
% C.archive.simulations{end}.sim.state(:,end) % use to extract steady state for new parameters.




%% If a unstable steady state is found via closed loop simulation, save with this:
sim_end = C.archive.simulations{end}.sim.state(:,end);
sim_inputs = C.archive.simulations{end}.sim.input_effective(:,end);
% (unstable) 
unstable_eq_saved = [structor; structor];
unstable_eq_saved(1).str.x_I    = sim_end(1);
unstable_eq_saved(1).str.x_M    = sim_end(2);
unstable_eq_saved(1).str.x_T    = sim_end(3);
unstable_eq_saved(1).str.x_Tc   = sim_end(4);
unstable_eq_saved(1).str.x_D0   = sim_end(5);
unstable_eq_saved(1).str.x_D1   = sim_end(6);
unstable_eq_saved(2).str.Q_c   = sim_inputs(1);
unstable_eq_saved(2).str.Q_m   = sim_inputs(2);
save('unstable_eq_saved',"unstable_eq_saved")


%% Check closed loop NMPC parameters:
dp = [];
for i = 1:length(C.archive.simulations{end}.optimizations)-1
   dp = [dp C.archive.simulations{end}.optimizations{i+1}.parameters.vec - C.archive.simulations{end}.optimizations{i}.parameters.vec];
end
disp(dp)