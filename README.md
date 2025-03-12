# Try-MPC
A MATLAB class for easily trying out NMPC control for your models (ODE/DAE) ! 

This class; "TRYMPC", is meant to be a super fast and easy way for anyone, beginners and experts, to try NMPC controllers for their systems.

It is written with respect to being very easy to use, analyze, debug, etc. All data is stored, such that the simulations are easily analyzed afterwards.

---
### (some) Capabilities:
 - Automatically generate dynamic optimization problem on a specified horizon, with any desired integrator (Euler, RK, collocation) and any approach (single shooting, multiple shooting, direct collocation).
 - Automatically generate casadi expressions and functions to evaluate all problem jacobians and hessians, and easily inspect sparsity patterns.
 - Solve open-loop optimization problems using your favorite solver.
 - Simulate your system using your favorite NMPC scheme.
   - also use any externally privded controller (PID, LQR, self-made NMPC, etc.)
 - Automatically generate nice plots of your solutions.
 - Store all data easily, and regenerate from .mat file.

---
### Some requirements:
- matlab must have CasADi on its path, such that using f.ex: "casadi.SX.sym(...)" works. CasADi: https://web.casadi.org/
- Some of my personal self-made function are used, so download my general matlab library and add to path: https://github.com/TrymArve/General-Matlab-Library
  - I will probably include the necessary functions directly into the class at some point, but just use the library for now.

---
### Tutorials
Try out the tutorials! Simply run each subsection in sequence and see what happens!
