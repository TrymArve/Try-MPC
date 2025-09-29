% Defines state, input and parameter names:
% state_names = ["x_I","x_M","x_T","x_Tc","x_D0","x_D1"];
state_names = ["I", "M", "T", "Tc", "D0", "D1"];
input_names = ["Qc","Qm"];
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
"Q_i";
"T_cf";
"T_f";
"V";
"V_c" ;
"nDHr";
"M_m"; 
"gamma";
"lambda"];



%%% Define Dynamics:
a_m = @(s,a,i,p)      (1/(1-p.str.lambda))*i.str.Qm/p.str.V;
a_gamma = @(s,a,i,p)  (1-p.str.gamma)/(1-p.str.lambda);
k_d = @(s,a,i,p)       p.str.A_d*exp(-p.str.E_d/s.str.T);
xi = @(s,a,i,p)        p.str.A_p*(2*p.str.f_eff*p.str.A_d/p.str.A_t)^(1/2)*exp((p.str.E_t-p.str.E_d-2*p.str.E_p)/(2*s.str.T));
MI = @(s,a,i,p)        s.str.M*s.str.I^(1/2);

state = @(s) [s.str.I, s.str.M, s.str.T, s.str.Tc, s.str.D0, s.str.D1]';

A_lin = @(s,a,i,p) -a_m(s,a,i,p)*diag([1 1 1 0 1 1])...
                 + [ 0 0 0 0 0 0;
                     0 0 0 0 0 0;
                     0 0 [-1 1].*( p.str.hA/(p.str.rhoC_p*p.str.V) ) 0 0;
                     0 0 [1 -1].*( p.str.hA/(p.str.rhoC_pc*p.str.V_c) ) 0 0;
                     0 0 0 0 0 0;
                     0 0 0 0 0 0];

f_nonlin = @(s,a,i,p)...
           [ -k_d(s,a,i,p)*s.str.I;
             -xi(s,a,i,p)*MI(s,a,i,p);
             p.str.nDHr/(p.str.rhoC_p)*xi(s,a,i,p)*MI(s,a,i,p);
             0;
             p.str.f_eff*k_d(s,a,i,p)*s.str.I;
             p.str.M_m*xi(s,a,i,p)*MI(s,a,i,p)];

C_const = @(s,a,i,p)...
          [0;
           i.str.Qm/p.str.V*p.str.M_f;
           a_m(s,a,i,p)*p.str.T_f;
           0; 0; 0];

B_control = @(s,a,i,p)...
            [1/p.str.V*(p.str.I_f-a_gamma(s,a,i,p)*s.str.I)*p.str.Q_i;
             -s.str.M/p.str.V*a_gamma(s,a,i,p)*p.str.Q_i;
             1/p.str.V*a_gamma(s,a,i,p)*(p.str.T_f-s.str.T)*p.str.Q_i;
             1/p.str.V_c*(p.str.T_cf-s.str.Tc)*i.str.Qc;
             -s.str.D0/p.str.V*a_gamma(s,a,i,p)*p.str.Q_i;
             -s.str.D1/p.str.V*a_gamma(s,a,i,p)*p.str.Q_i];
% dynamics (dq = f(q,u)):
CSTR_dynamics = @(s,a,i,p) A_lin(s,a,i,p)*state(s) + f_nonlin(s,a,i,p) + C_const(s,a,i,p) + B_control(s,a,i,p);


%%% Initiate class with variable names:
M = trympcMODEL('CSTR_ex1',...
    state = state_names, ...
    input = input_names,...
    param = param_names,...
    dynamics = CSTR_dynamics);