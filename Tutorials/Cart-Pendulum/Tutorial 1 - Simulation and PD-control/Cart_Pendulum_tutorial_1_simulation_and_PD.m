%% Cart-Pendulum Tutorial 1 - Defining dynamics and simulating


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



%% Initiate class with variable names:
fresh % start fresh

C = TRYMPC('Tester Instance',...
    state = ["x" "th" "dx" "dth"], ... % provide state variable names
    input = "ux", ...                  % provide input variable names
    param = ["L" "g" "mx" "mth"])      % provide parameter names (keep them symbolic, to avoid hard-coding them)


%% Define dynamics:
u_pendulum = 0; % no input on the pendulum

% cart acceleration
ddx = @(s,a,i,p) -(p.L*i.ux + u_pendulum*cos(s.th) + p.L^2*s.dth^2*p.mth*sin(s.th) - p.L*p.g*p.mth*cos(s.th)*sin(s.th))/(p.L*(- p.mth*cos(s.th)^2 + p.mx + p.mth));

% pendulum angular acceleration
ddth = @(s,a,i,p) -(p.mx*u_pendulum+ p.mth*u_pendulum + p.L*p.mth*i.ux*cos(s.th) - p.L*p.g*p.mth^2*sin(s.th) + p.L^2*s.dth^2*p.mth^2*cos(s.th)*sin(s.th) - p.L*p.g*p.mx*p.mth*sin(s.th))/(p.L^2*p.mth*(- p.mth*cos(s.th)^2 + p.mx + p.mth));

% dynamics (dq = f(q,u))
dynamics = @(s,a,i,p) [s.dx;
                       s.dth;
         	           ddx(s,a,i,p);
                       ddth(s,a,i,p)];


% Apply dynamics to TRYMPC instance:
C.def_dynamics(dynamics)


%% Define LaTex Display names for nice plots:
C.plotting.display_names.state.x = "$x$";
C.plotting.display_names.state.th = "$\theta$";
C.plotting.display_names.state.dx = "$\dot{x}$";
C.plotting.display_names.state.dth = "$\dot{\theta}$";
C.plotting.display_names.input.ux = "$u_{x}$";


%% Prepare numerical values for parameters

%{
Note that these values can be changed at any time in order to update the
system model, even during simulation.

The datatype for the parameters is a specialized class "structor", which can
act both as a "struct" and as a "vector". This is not important here, but
is very convenient in this class, and is used frequently.
%}

C.parameters.str.L = 1.1;
C.parameters.str.g = 9.81;
C.parameters.str.mx = 5;
C.parameters.str.mth = 3;

%% Choose an initial state

%{
here we also use a structor, though this is unnecessary here, we just need
a vector. Using "init_state.vec" gives the values as a vector.

One may also simply do this:
init_state = [0.1; 0.2; 0; 0];
%}

init_state = structor;
init_state.str.x   = 0.2;
init_state.str.th  = 0.2;
init_state.str.dx  = 0;
init_state.str.dth = 0;


%% Try simulating:

duration = 10; % simulate for 10 seconds

C.simulate(duration,init_state.vec)

%% We now have simulation stored in the archive, have a look:

C.archive.simulations{1}

%% We can also plot the results of the simulation in a neat figure:

C.display_simulation


%% Try a different simulation:

duration = 15; % simulate for 15 seconds

% Why not increase gravity?
C.parameters.str.g = 11.6;

% Choose different initial values for x and th:
init_state.str.x   = 0.1;
init_state.str.th  = 0.4;

C.simulate(duration,init_state.vec, simulator="ode113") % specify "ode113" for fun... (The default simulator is "ode45", and can be changed)

%% The new simulation is now stored as well !

C.archive.simulations{2}


%% Let plot the new simulation too, but only the states

C.display_simulation("state","all")


%% Forgot the first simulation results? You can still plot it. Let's only plot "x" and "th" this time...

C.display_simulation("simulation_number",1,"state",["x","th"])


%% Let's define a simple PD controller

x_kp =  -1;
x_kd =  -1;
th_kp = 100;
th_kd = 10;

controller_handle = @(~,~,x) x(1)*x_kp + x(2)*th_kp + x(3)*x_kd + x(4)*th_kd;


%% Now, lets simulate the PD-controlled system

C.parameters.str.L = 1.1;  % change the pendulum length
C.parameters.str.g = 9.81; % reset gravity

% some other initial conditions:
init_state.str.th    = 0.2;
init_state.str.dth   = 0.1;

C.simulate(50,init_state.vec,...   
   simulator="ode15s",...                    % Choose the "ode15s" simulator for fun...
   controller_type="custom",...              % Choose controller type ("custom" means that you provide an external controller)
   controller_custom=controller_handle)      % "custom"-specific: apply controller handle here

C.display_simulation; % display solution


%% Zero-order hold (piecewise constant with given sampling time)
%{
Until now, we have simply augmented the system dynamics with a PD
controller.
In reality, we must measure the state of the system, then compute a control
input. This input is typically kept constant until the updated input is
available. This is called zero-order hold, or simply piece-wise constant.

Let try simulating this, slightly more realistic, scenario.
%}

C.parameters.str.L = 5; % Let's update the pendulum length again

% Choose an initial state that is not too far from 0
init_state.str.x   = 0.01;
init_state.str.th  = 0;
init_state.str.dx  = 0;
init_state.str.dth = 0;

% Update the PD controller:
x_kp =  -100;
x_kd =  -30;
th_kp = 1500;
th_kd = 300;
controller = @(C,t,x) [x_kp th_kp x_kd th_kd]*x;

% Choose a sampling time (time of measuring state, and thus updating the control input)
samp_time = 0.07;

duration = 15; % choose a duration of 15s

% Let's configue the ode-solver (default: "ode45") tolerances just in case:
ode_opts = odeset('RelTol',10^(-7),'AbsTol',10^(-7));

% Simulate
C.simulate(duration,init_state.vec,...   
   controller_type   = "custom",...    
   controller_custom = controller,...
   sampling_time     = samp_time,...     % Set the sampling time (this automatically causes a zero-order hold controller)
   ode_options       = ode_opts);        % provide external ode options

% save the simulation index, so we can plot this again later in the tutorial
index_without_control_delay = length(C.archive.simulations);

%% Now lets plot the results:

C.display_simulation("gridstyle","vertical",... % I want the plots to be stacked vertically
                     "state",["x","th"],...     % Let's only plot the states: "x" and "th"
                     "mark_samples",true);      % Let's add a marker on the sampled measurements


%% And the input signal: (notice that it's piecewise constant)

C.display_simulation("input","all");     % Only plot the input


%% Now we will try to be even more realistic!

%{
Let's account for the time it takes to compute the control signal after
measuring the state (sampling).

This is called control delay, and causes the state to change accoring to
the previous input signal for a short duration, while the new signal is
being computed.

This casues an effect where, when finally applying the control signal, the
system is not longer in the state that control signal was based on. This
means that the control signal is now 'wrong', in some sense. 
This effect can cause instability.

Perticularly for NMPC, the control delay is often
very large, since NMPC is computaionally demanding.
%}

close all

% Change to a non-anonymous function, such that we can time exactly how long it takes to
% compute the control signal:
controller = @(C,t,x) PD(C,t,x,[x_kp th_kp x_kd th_kd]);

% To get the point accross, we will pretend that the control delay is
% larger than is really is:
scale_cd = 1e10; % scale control delay by 1e10. 
%{
Depending on your computer, the control delay will vary since the
computation time varies dpeneding on computing power of your computer.
Also, if the delays are very small, and the system is sensitive to the control delays, random fluctuations in computation time will have
larger impact on the performance.

Therefore, depending on your computer speed, you may want to play with the
control delay scaler.

In this tutorial, we scale the control delay by a stupidly large amount
(1e10) just to ensure that everyone sees the same delays.
The control delay is capped at the sampling time, thus when scaling the
control delay so much, we guarantee that everyone simply sees the sampling
time as the control delay. 
This means that the control signal always lags one sample cycle behind.

NB: this will cause a bunch of warning to appear. Ignore them!
%}


% Simulate
C.simulate(duration,init_state.vec,...   
   controller_type      = "custom",...    
   controller_custom    = controller,...
   sampling_time        = samp_time,...     
   ode_options          = odeset('RelTol',10^(-7),'AbsTol',10^(-7)),...
   control_delay_scaler = scale_cd); % Set the control delay scaler to any value in greater than zero to account for control delay. (F.ex. setting it to 0.7 scales the control delay down to 70%)

% save the simulation index
index_with_control_delay = length(C.archive.simulations);


% Display previous and new result on top of each other and compare the effect of control delay:
C.display_simulation("simulation_number",[index_without_control_delay,index_with_control_delay],"state",["x","th"],"input","ux");





%% Functions:

function u = PD(C,~,x,K)

   % Measure how much time is required to compute the control signal:
   tic;
   u = K*x;
   control_delay = toc;

   % Tell the TRYMPC class instance how long you spend computing the
   % control signal:
   C.set_control_delay(control_delay); 
end