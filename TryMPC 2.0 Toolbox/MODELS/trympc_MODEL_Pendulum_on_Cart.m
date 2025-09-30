


state_names = ["cart_pos" "pendulum_angle" "cart_speed" "pendulum_speed"]; % state variable names
input_names = "cart_force";                                               % input variable names
param_names = ["pendulum_length" "gravity" "mass_cart" "mass_pendulum"];  % parameter names (keep them symbolic, to avoid hard-coding them)

%%% Define dynamics:
u_pendulum = 0;
dynamics = @(s,i,p) [s.str.cart_speed;
                       s.str.pendulum_speed;
                       -(p.str.pendulum_length*i.str.cart_force + u_pendulum*cos(s.str.pendulum_angle) + p.str.pendulum_length^2*s.str.pendulum_speed^2*p.str.mass_pendulum*sin(s.str.pendulum_angle) - p.str.pendulum_length*p.str.gravity*p.str.mass_pendulum*cos(s.str.pendulum_angle)*sin(s.str.pendulum_angle)) / ...
                        (p.str.pendulum_length*(- p.str.mass_pendulum*cos(s.str.pendulum_angle)^2 + p.str.mass_cart + p.str.mass_pendulum));
                       -(p.str.mass_cart*u_pendulum + p.str.mass_pendulum*u_pendulum + p.str.pendulum_length*p.str.mass_pendulum*i.str.cart_force*cos(s.str.pendulum_angle) - p.str.pendulum_length*p.str.gravity*p.str.mass_pendulum^2*sin(s.str.pendulum_angle) + ...
                         p.str.pendulum_length^2*s.str.pendulum_speed^2*p.str.mass_pendulum^2*cos(s.str.pendulum_angle)*sin(s.str.pendulum_angle) - p.str.pendulum_length*p.str.gravity*p.str.mass_cart*p.str.mass_pendulum*sin(s.str.pendulum_angle)) / ...
                        (p.str.pendulum_length^2*p.str.mass_pendulum*(- p.str.mass_pendulum*cos(s.str.pendulum_angle)^2 + p.str.mass_cart + p.str.mass_pendulum))];

% %%% Define dynamics:
% u_pendulum = 0;
% dynamics = @(s,a,i,p) [s.str.dx;
%                        s.str.dth;
%                        -(p.str.L*i.str.ux + u_pendulum*cos(s.str.th) + p.str.L^2*s.str.dth^2*p.str.mth*sin(s.str.th) - p.str.L*p.str.g*p.str.mth*cos(s.str.th)*sin(s.str.th))/(p.str.L*(- p.str.mth*cos(s.str.th)^2 + p.str.mx + p.str.mth));
%                        -(p.str.mx*u_pendulum + p.str.mth*u_pendulum + p.str.L*p.str.mth*i.str.ux*cos(s.str.th) - p.str.L*p.str.g*p.str.mth^2*sin(s.str.th) + p.str.L^2*s.str.dth^2*p.str.mth^2*cos(s.str.th)*sin(s.str.th) - p.str.L*p.str.g*p.str.mx*p.str.mth*sin(s.str.th))/(p.str.L^2*p.str.mth*(- p.str.mth*cos(s.str.th)^2 + p.str.mx + p.str.mth))];

M = trympcMODEL("Pendulum on Cart (ODE)",...
   "state",state_names,...
   "input",input_names,...
   "param",param_names,...
   "dynamics",dynamics);

param = structor;
param.str.pendulum_length = 0.5;
param.str.gravity         = 9.81;
param.str.mass_cart       = 2;
param.str.mass_pendulum       = 2;

init_state = structor;
init_state.str.cart_pos       = 0;
init_state.str.pendulum_angle = 0.01;
init_state.str.cart_speed     = 0;
init_state.str.pendulum_speed = 0;

uneq_state = structor;
uneq_state.str.cart_pos       = 0;
uneq_state.str.pendulum_angle = 0;
uneq_state.str.cart_speed     = 0;
uneq_state.str.pendulum_speed = 0;

ref.state = @(t) zeros(4,1).*t;
ref.input = @(t) 0.*t;

M.numeric_model = trympcNUMERIC_MODEL("Pendulum Swingup",param,...
     "initial_state",init_state,...
     "unstable_equilibrium_state",uneq_state,...
     "ref",ref);