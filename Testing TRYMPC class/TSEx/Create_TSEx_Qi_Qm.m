%%% Set up stuff:
state_names = ["x_I","x_M","x_T","x_Tc","x_D0","x_D1"];
input_names = ["Q_c","Q_i","Q_m"];
param_names = [
"rhoC_p";
"rhoC_pc";
"f_eff";
"hA";
"I_f";
"A_d";
"A_p";
"A_t";
"E_d";
"E_p";
"E_t";
"M_f";
%"Q_m";
"T_cf";
"T_f";
"V";
"V_c" ;
"nDHr";
"M_m"; 
"gamma";
"lambda"];

minutes = @(m) m*60;
hours = @(h) h*3600;
days  = @(d) d*24*3600;


%%% Initiate class with variable names:
C = TRYMPC('TSEx (Qm control)',...
    state = state_names, ...
    input = input_names,...
    param = param_names);
    
%%% Define Dynamics:
a_m = @(s,a,i,p)      (1/(1-p.lambda))*i.Q_m/p.V;
a_gamma = @(s,a,i,p)  (1-p.gamma)/(1-p.lambda);
k_d = @(s,a,i,p)       p.A_d*exp(-p.E_d/s.x_T);
xi = @(s,a,i,p)        p.A_p*(2*p.f_eff*p.A_d/p.A_t)^(1/2)*exp((p.E_t-p.E_d-2*p.E_p)/(2*s.x_T));
xMI = @(s,a,i,p)       s.x_M*s.x_I^(1/2);

state = @(s) [s.x_I, s.x_M, s.x_T, s.x_Tc, s.x_D0, s.x_D1]';

A_lin = @(s,a,i,p) -a_m(s,a,i,p)*diag([1 1 1 0 1 1])...
                 + [ 0 0 0 0 0 0;
                     0 0 0 0 0 0;
                     0 0 [-1 1].*( p.hA/(p.rhoC_p*p.V) ) 0 0;
                     0 0 [1 -1].*( p.hA/(p.rhoC_pc*p.V_c) ) 0 0;
                     0 0 0 0 0 0;
                     0 0 0 0 0 0];

f_nonlin = @(s,a,i,p)...
           [ -k_d(s,a,i,p)*s.x_I;
             -xi(s,a,i,p)*xMI(s,a,i,p);
             p.nDHr/(p.rhoC_p)*xi(s,a,i,p)*xMI(s,a,i,p);
             0;
             p.f_eff*k_d(s,a,i,p)*s.x_I;
             p.M_m*xi(s,a,i,p)*xMI(s,a,i,p)];

C_const = @(s,a,i,p)...
          [0;
           i.Q_m/p.V*p.M_f;
           a_m(s,a,i,p)*p.T_f;
           0; 0; 0];

B_control = @(s,a,i,p)...
            [1/p.V*(p.I_f-a_gamma(s,a,i,p)*s.x_I)*i.Q_i;
             -s.x_M/p.V*a_gamma(s,a,i,p)*i.Q_i;
             1/p.V*a_gamma(s,a,i,p)*(p.T_f-s.x_T)*i.Q_i;
             1/p.V_c*(p.T_cf-s.x_Tc)*i.Q_c;
             -s.x_D0/p.V*a_gamma(s,a,i,p)*i.Q_i;
             -s.x_D1/p.V*a_gamma(s,a,i,p)*i.Q_i];

% dynamics (dq = f(q,u)):
dynamics = @(s,a,i,p) A_lin(s,a,i,p)*state(s) + f_nonlin(s,a,i,p) + C_const(s,a,i,p) + B_control(s,a,i,p);

% Apply to Trympc:
C.def_dynamics(dynamics)


%%% Displaying:
% states:
C.plotting.display_names.state.x_I     = "$x_I$";
C.plotting.display_names.state.x_M     = "$x_M$";
C.plotting.display_names.state.x_T     = "$x_T$";
C.plotting.display_names.state.x_Tc    = "$x_{Tc}$";
C.plotting.display_names.state.x_D0    = "$x_{D_0}$";
C.plotting.display_names.state.x_D1    = "$x_{D_1}$";
C.plotting.color.state.x_I  = GetColorCode("g");
C.plotting.color.state.x_M  = GetColorCode("p");
C.plotting.color.state.x_T  = GetColorCode("b");
C.plotting.color.state.x_Tc = GetColorCode("c");
C.plotting.color.state.x_D0 = GetColorCode("o");
C.plotting.color.state.x_D1 = GetColorCode("r");
% inputs:
C.plotting.display_names.input.Q_c    = "$Q_{c}$";
C.plotting.display_names.input.Q_i    = "$Q_{i}$";
C.plotting.display_names.input.Q_m    = "$Q_{m}$";
C.plotting.color.input.Q_c    = GetColorCode("c");
C.plotting.color.input.Q_i    = GetColorCode("g");
C.plotting.color.input.Q_m    = GetColorCode("p");