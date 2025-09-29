%%% Create a CSTR model for the with the TRYMPC class

fresh

%% Model

M_CSTR = trympcmodels("CSTR_ex1");

%% Discretizer (Integrator)

% Some varieties:
D_ERK4 = trympcDISCRETIZER("ERK for CSTR",M_CSTR,"ERK4");
D_IRK4 = trympcDISCRETIZER("IRK (L-stable) for CSTR",M_CSTR,"IRK4 (L-stable)");
D_GL4  = trympcDISCRETIZER("Gauss-Legendre for CSTR",M_CSTR,"Gauss-Legendre (4. order)");
D_CLC4 = trympcDISCRETIZER("Collocation for CSTR",M_CSTR,"custom collocation",...
          "collocation_polynomial_order",2,...
          "collocation_polynomial_type","legendre");




%% 

