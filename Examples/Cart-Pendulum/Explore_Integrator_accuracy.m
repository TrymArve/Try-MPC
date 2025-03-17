%% Example: Cart-Pendulum - Explore accuracy of integrators

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

C.def_dynamics(dynamics)                     % define dynamics
C.def_objective(quadratic=["Q","R","dR","Q_terminal"])    % define objective
C.def_stage_constraints("upper_bounds","ux","lower_bounds","ux") % define stage constraints
% C.def_terminal_constraint("equality",C.cas.state) % must terminate at the origin
N = 30; % Choose a high resolution grid for the high accuracy solutions


%% Provide numeric values of your problem

% Parameters
C.parameters.str.L = 1;
C.parameters.str.g = 9.81;
C.parameters.str.mx = 5;
C.parameters.str.mth = 3;

% Initial state (cart-offset by 1)
C.initial_state.str.x   = 1;
C.initial_state.str.th  = 0;
C.initial_state.str.dx  = 0;
C.initial_state.str.dth = 0;

% quadratic cost penalities
C.quadratic_cost.Q = [1 10 0.1 0.1]; % state penalty
C.quadratic_cost.R = 1;              % input penalty
C.quadratic_cost.dR = 1;             % change in input penalty
C.quadratic_cost.Q_terminal = C.quadratic_cost.Q*1; % terminal state penalty

% upper bound on input
C.bounds.lower.ux = -0.4;
C.bounds.upper.ux =  0.4;

% Define the horizon length in terms of how far into the future you want to
% predict (given in seconds):
C.set_T(20); % predict T seconds into the future

% Define LaTex Display names:
C.plotting.display_names.state.x = "$x$";
C.plotting.display_names.state.th = "$\theta$";
C.plotting.display_names.state.dx = "$\dot{x}$";
C.plotting.display_names.state.dth = "$\dot{\theta}$";
C.plotting.display_names.input.ux = "$u_{x}$";

%% First define the true optimal continuous solution: (we pretend that this is a continuous solution)
N_high = 200;
inc = 3;
C.def_integrator("ERK4","n_increments",inc) % ERK with many steps!
C.def_horizon(N_high)
C.solve("ipopt");
index_truth = length(C.archive.optimizations);
name_truth = "Continuous sol.";


%% Generate High Accuracy solutions to function as the "true" solution (with respect to the given discretization/transcription)

inc = 10;
C.def_integrator("ERK4","n_increments",inc) % ERK with many steps!
C.def_horizon(N)
C.solve("ipopt");
index_high_ERK = length(C.archive.optimizations);
name_high_ERK = "Discrete sol.";

%% Gauss-Legendre 6. order (a different integration scheme to ensure that they all give the same answer, reassuring us that it actually the true solution)
inc = 5;
C.def_integrator("Gauss-Legendre (6. order)","n_increments",inc) % ERK with many steps!
C.def_horizon(N)
C.solve("ipopt");
index_high_GL6 = length(C.archive.optimizations);
name_high_GL6 = append("Gauss-Legendre 6. (inc:",string(inc),", N = ",string(N),")");

%% Collocation (and one more integrator)
inc = 1;
d = 3;
C.def_integrator("collocation","n_increments",inc,"collocation_polynomial_order",d) % ERK with many steps!
C.def_horizon(N)
C.solve("ipopt");
index_high_collocation = length(C.archive.optimizations);
name_high_collocation = append("Collocation (d:",string(d),",inc:",string(inc),", N = ",string(N),")");


%% Plot thogether

% They are all equal, demonstrating that the solutions are accurate!
Ind = [index_truth,index_high_ERK,index_high_GL6,index_high_collocation];
Names = [append(name_truth," (erk4, inc:",string(inc),", N = ",string(N_high),")"),...
         append(name_high_ERK," (erk4, inc:",string(inc),", N = ",string(N),")"),name_high_GL6,name_high_collocation];
Styles = ["-","-.","--","-"];
C.display_optimization("optimization_number",Ind,"legend",Names,"linestyle",Styles,"multiplot","ontop");


%{
Note that all three trajectories on the horizon are somewhat inaccurate
because of the discretization of time into stages (transcription).
Therefore the true optimal trajectory is smooth and does not perfectly
align with the diescrete trajectories.

By integrating with three different high accuracy integrators, and seeing
that thay all produce the same result, we are reassured that the solutions
are in fact very accurate, and can be considered the true solution.

Now we can move on to lower accuracy integrators, to see how well they
perform. 
%}














%%%%%%%%%%%%%%%%%%%% Compare Low accuracy solutions:

%% Now use a finer grid, such that it is faster to compute, but is less accurate

% Test a bunch of integrators:
[ind,names,style,color] = runall(C,"ipopt",300,N);


%% Plot all

   % Add True solution:
   Ind = [index_truth index_high_ERK ind];
   Names = [name_truth name_high_ERK names];
   Style = ["-" "-" style];
   Color = ["grey" "black" color];

C.display_optimization("optimization_number",Ind,"legend",Names,"linestyle",Style,"colors",Color,"multiplot","ontop","color_match","solution");


%% Stats
for i = 1:length(Ind)
   j = Ind(i);
   disp('--------------------------')
   disp(['   Integrator: ',char(Names(i))])
   disp(['       solver: ',char(C.archive.optimizations{j}.solver) ])
   disp(['   solve time: ',num2str(C.archive.optimizations{j}.solve_time.total)])
   disp(['n. iterations: ',num2str(C.archive.optimizations{j}.n_iterations)])
   disp(['return status: ',num2str(C.archive.optimizations{j}.return_status)])
   disp(['            N: ',num2str(C.archive.optimizations{j}.N)])
   disp( '  Convergence:'); present(C.archive.optimizations{j}.convergence)
end
disp('--------------------------')



%% Varying Max N. Iterations

% Choose an integrator to test:

% C.def_integrator("Explicit Euler","n_increments",1)
% C.def_integrator("Implicit Euler","n_increments",1)
% C.def_integrator("ERK4","n_increments",1)
% C.def_integrator("ERK4 (simultaneous)","n_increments",1)
% C.def_integrator("Implicit Midpoint","n_increments",1)
% C.def_integrator("Crank-Nicolson (Implicit)","n_increments",1)
% C.def_integrator("IRK4 (L-stable)","n_increments",1)
% C.def_integrator("Gauss-Legendre (4. order)","n_increments",1)
% C.def_integrator("Gauss-Legendre (6. order)","n_increments",1)
C.def_integrator("collocation","n_increments",1,"collocation_polynomial_order",2)


C.def_horizon(N)
Max_iter = [1 2 3 4 5 300];
ind_max = varying_iter(C,"ipopt",Max_iter);

C.display_optimization("optimization_number",ind_max,"legend",Max_iter);
















%% Functions

function ind = varying_iter(C,solver,max_iter)
   ind = [];
   for m_iter = max_iter
      C.solve(solver,max_iterations=m_iter);
      ind(end+1) = length(C.archive.optimizations);
   end
end


function [ind,names,style,color] = runall(C,solver,max_iter,N)

ind = [];
names = string([]);
style = string([]);
color = string([]);

   %
   C.def_integrator("Explicit Euler","n_increments",1)
   C.def_horizon(N)
   C.solve(solver,max_iterations=max_iter);
   ind(end+1) = length(C.archive.optimizations);
   names(end+1) = "Explicit Euler";
   style(end+1) = ":";
   color(end+1) = "yellow";
   %
   C.def_integrator("Implicit Euler","n_increments",1)
   C.def_horizon(N)
   C.solve(solver,max_iterations=max_iter);
   ind(end+1) = length(C.archive.optimizations);
   names(end+1) = "Implicit Euler";
   style(end+1) = ":";
   color(end+1) = "orange";
   %
   C.def_integrator("ERK4","n_increments",1)
   C.def_horizon(N)
   C.solve(solver,max_iterations=max_iter);
   ind(end+1) = length(C.archive.optimizations);
   names(end+1) = "ERK4";
   style(end+1) = "-.";
   color(end+1) = "green";
   %
   C.def_integrator("ERK4 (simultaneous)","n_increments",1)
   C.def_horizon(N)
   C.solve(solver,max_iterations=max_iter)
   ind(end+1) = length(C.archive.optimizations);
   names(end+1) = "ERK4 (sim)";
   style(end+1) = "-.";
   color(end+1) = "grey";
   %
   C.def_integrator("Implicit Midpoint","n_increments",1)
   C.def_horizon(N)
   C.solve(solver,max_iterations=max_iter);
   ind(end+1) = length(C.archive.optimizations);
   names(end+1) = "Implicit Midpoint";
   style(end+1) = "--";
   color(end+1) = "pink";
   %
   C.def_integrator("Crank-Nicolson (Implicit)","n_increments",1)
   C.def_horizon(N)
   C.solve(solver,max_iterations=max_iter);
   ind(end+1) = length(C.archive.optimizations);
   names(end+1) = "Crank-Nicolson";
   style(end+1) = "--";
   color(end+1) = "magenta";
   %
   C.def_integrator("IRK4 (L-stable)","n_increments",1)
   C.def_horizon(N)
   C.solve(solver,max_iterations=max_iter);
   ind(end+1) = length(C.archive.optimizations);
   names(end+1) = "IRK4 (L-stable)";
   style(end+1) = "--";
   color(end+1) = "cyan";
   %
   C.def_integrator("Gauss-Legendre (4. order)","n_increments",1)
   C.def_horizon(N)
   C.solve(solver,max_iterations=max_iter);
   ind(end+1) = length(C.archive.optimizations);
   names(end+1) = "Gauss-Leg. 4";
   style(end+1) = "-";
   color(end+1) = "red";
   %
   C.def_integrator("Gauss-Legendre (6. order)","n_increments",1)
   C.def_horizon(N)
   C.solve(solver,max_iterations=max_iter);
   ind(end+1) = length(C.archive.optimizations);
   names(end+1) = "Gauss-Leg. 6";
   style(end+1) = "-";
   color(end+1) = "purple";
   %
   pol_order = 2;
   C.def_integrator("collocation","n_increments",1,"collocation_polynomial_order",pol_order)
   C.def_horizon(N)
   C.solve(solver,max_iterations=max_iter);
   ind(end+1) = length(C.archive.optimizations);
   names(end+1) = append("Collocation (d:",string(pol_order),")");
   style(end+1) = "-";
   color(end+1) = "blue";

end