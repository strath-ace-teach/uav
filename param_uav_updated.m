function P  = param_uav_updated(alpha)
%{
ME528 Group Assignment: Aerosonde UAV parameters
 

The following returns a structure P which contains a large number of
parameters for the UAV including mass, physical dimensions, aerodynamic
coefficients. Due to the dependancy of the lift and drag coefficients to
the angle of attack of the UAV, some of the coefficients will change with
time.  A description of these variables and some of the background theory
is explained in the Project Description pdf on MyPlace. 

INPUTS: 
alpha = angle of attack (rad) of the UAV

OUTPUTS:
P = structure of all parameters for the Aerosonde UAV vehicle design and
operation
%}
 clear all
 clc
%% Physical parameters of airframe
mass = 13.5; %mass, kg
S_wing        = 0.55;  %wing area, m^2
b             = 2.8956;  %wing span, m
c             = 0.18994;  %, mean aerodynamic chord, m
e             = 0.9;  %Oswald efficiency
g = 9.80665; % gravitational acceleration, m/s^2
F_prop = 1.2*13.5*g;  %max thrust from propeller (** specific to this assignment)
AR = (b^2)/S_wing;  % aspect ratio

% moments of inertia (matrix J)
Jx   = 0.8244;  %kg-m^2
Jy   = 1.135; %kg-m^2
Jz   = 1.759; %kg-m^2
Jxz  = 0.1204; %kg-m^2

rho=1.225;
Va=20;
eve=(rho.*Va^2.*S_wing/2); 

Gamma = Jx*Jz - Jxz^2;
Gamma_1 = (Jxz*(Jx-Jy+Jz))/Gamma;
Gamma_2 = (Jz*(Jz-Jy)+Jxz^2)/Gamma;
Gamma_3 = Jz/Gamma;
Gamma_4 = Jxz/Gamma;
Gamma_5 = (Jz - Jx)/Jy;
Gamma_6 = Jxz/Jy;
Gamma_7 = ((Jx-Jy)*Jx + Jxz^2)/Gamma;
Gamma_8 = Jx/Gamma;

%% Aerodynamics as a function of angle of attack, alpha
% Use either Model A *or* Model B, just comment or uncomment as appropriate

M             = 50;
alpha_0        = 0.4712;
C_L_0         = 0.28;
C_L_alpha     = 3.45;
C_L_q         = 0.0;
C_L_delta_e   = -0.36;
C_D_0         = 0.03;
C_D_alpha     = 0.30;
C_D_p         = 0.0437;
C_D_q         = 0.0;
C_D_delta_e   = 0.0;

% MODEL A for CL, CD: Simplified aerodynamic model for small angle of attack
% C_L = C_L_0 + C_L_alpha*alpha;
% C_D = C_D_0 + C_D_alpha*alpha;

% MODEL B for CL, CD: Incorporates stall
sigma = ( 1 + exp(-M*(alpha - alpha_0)) +  exp(M*(alpha + alpha_0))) ...
    / ( (1+exp(-M*(alpha - alpha_0)))*(1+exp(M*(alpha + alpha_0))) );
C_L = (1-sigma)*(C_L_0+C_L_alpha*alpha) + sigma*(2*sign(alpha)*sin(alpha)^2*cos(alpha));
C_D = C_D_p + (C_L_0+C_L_alpha*alpha)^2/(pi*e*AR);

%% Aerodynamic coefficients 
C_prop        = 1.0;

% [m, n, l=ell]
C_m_0         = -0.02338;
C_m_alpha     = -0.38;
C_m_q         = -3.6;
C_m_delta_e   = -0.5;

C_ell_0       = 0.0;
C_ell_beta    = -0.12;
C_ell_p       = -0.26;
C_ell_r       = 0.14;
C_ell_delta_a = 0.08;
C_ell_delta_r = 0.105;

C_n_0         = 0.0;
C_n_beta      = 0.25;
C_n_p         = 0.022;
C_n_r         = -0.35;
C_n_delta_a   = 0.06;
C_n_delta_r   = -0.032;

% [p, q, r]
C_p_0 = Gamma_3*C_ell_0 + Gamma_4*C_n_0;
C_p_beta = Gamma_3*C_ell_beta + Gamma_4*C_n_beta;
C_p_p = Gamma_3*C_ell_p + Gamma_4*C_n_p;
C_p_r = Gamma_3*C_ell_r + Gamma_4*C_n_r;
C_p_delta_a = Gamma_3*C_ell_delta_a + Gamma_4*C_n_delta_a;
C_p_delta_r = Gamma_3*C_ell_delta_r + Gamma_4*C_n_delta_r;

C_r_0 = Gamma_4*C_ell_0 + Gamma_8*C_n_0;
C_r_beta = Gamma_4*C_ell_beta + Gamma_8*C_n_beta;
C_r_p = Gamma_4*C_ell_p + Gamma_8*C_n_p;
C_r_r = Gamma_4*C_ell_r + Gamma_8*C_n_r;
C_r_delta_a = Gamma_4*C_ell_delta_a + Gamma_8*C_n_delta_a;
C_r_delta_r = Gamma_4*C_ell_delta_r + Gamma_8*C_n_delta_r;

% [X, Y, Z]
C_X = -C_D*cos(alpha) + C_L*sin(alpha);
% C_X_0 = -C_D_0*cos(alpha) + C_L_0*sin(alpha);
% C_X_alpha = -C_D_alpha*cos(alpha) + C_L_alpha*sin(alpha);
C_X_q =  -C_D_q*cos(alpha) + C_L_q*sin(alpha);
C_X_delta_e =  -C_D_delta_e*cos(alpha) + C_L_delta_e*sin(alpha);

C_Z = -C_D*sin(alpha) - C_L*cos(alpha);
% C_Z_0= -C_D_0*sin(alpha) - C_L_0*cos(alpha);
% C_Z_alpha = -C_D_alpha*sin(alpha) + C_L_alpha*cos(alpha);
C_Z_q =  -C_D_q*sin(alpha) - C_L_q*cos(alpha);
C_Z_delta_e =  -C_D_delta_e*sin(alpha) - C_L_delta_e*cos(alpha);

C_Y_0         = 0.0;
C_Y_beta      = -0.98;
C_Y_p         = 0.0;
C_Y_r         = 0.0;
C_Y_delta_a   = 0.0;
C_Y_delta_r   = -0.17;
 

%% Physical parameters of airframe
P.mass = 13.5; %mass, kg
P.S_wing        = 0.55;  %wing area, m^2
P.b             = 2.8956;  %wing span, m
P.c             = 0.18994;  %, mean aerodynamic chord, m
P.e             = 0.9;  %Oswald efficiency
g = 9.80665; % gravitational acceleration, m/s^2
P.F_prop = 1.2*13.5*g;  %max thrust from propeller (** specific to this assignment)
P.AR = (b^2)/S_wing;  % aspect ratio

% moments of inertia (matrix J)
P.Jx   = 0.8244;  %kg-m^2
P.Jy   = 1.135; %kg-m^2
P.Jz   = 1.759; %kg-m^2
P.Jxz  = 0.1204; %kg-m^2

Gamma = Jx*Jz - Jxz^2;
P.Gamma_1 = (Jxz*(Jx-Jy+Jz))/Gamma;
P.Gamma_2 = (Jz*(Jz-Jy)+Jxz^2)/Gamma;
P.Gamma_3 = Jz/Gamma;
P.Gamma_4 = Jxz/Gamma;
P.Gamma_5 = (Jz - Jx)/Jy;
P.Gamma_6 = Jxz/Jy;
P.Gamma_7 = ((Jx-Jy)*Jx + Jxz^2)/Gamma;
P.Gamma_8 = Jx/Gamma;

%% Aerodynamics as a function of angle of attack, alpha
% Use either Model A *or* Model B, just comment or uncomment as appropriate

P.M             = 50;
P.alpha_0        = 0.4712;
P.C_L_0         = 0.28;
P.C_L_alpha     = 3.45;
P.C_L_q         = 0.0;
P.C_L_delta_e   = -0.36;
P.C_D_0         = 0.03;
P.C_D_alpha     = 0.30;
P.C_D_p         = 0.0437;
P.C_D_q         = 0.0;
P.C_D_delta_e   = 0.0;

% MODEL A for CL, CD: Simplified aerodynamic model for small angle of attack
% P.C_L = P.C_L_0 + P.C_L_alpha*alpha;
% P.C_D = P.C_D_0 + P.C_D_alpha*alpha;

% MODEL B for CL, CD: Incorporates stall
sigma = ( 1 + exp(-M*(alpha - alpha_0)) +  exp(M*(alpha + alpha_0))) ...
    / ( (1+exp(-M*(alpha - alpha_0)))*(1+exp(M*(alpha + alpha_0))) );
P.C_L = (1-sigma)*(C_L_0+C_L_alpha*alpha) + sigma*(2*sign(alpha)*sin(alpha)^2*cos(alpha));
P.C_D = C_D_p + (C_L_0+C_L_alpha*alpha)^2/(pi*e*AR);

%% Aerodynamic coefficients 
P.C_prop        = 1.0;

% [m, n, l=ell]
P.C_m_0         = -0.02338;
P.C_m_alpha     = -0.38;
P.C_m_q         = -3.6;
P.C_m_delta_e   = -0.5;

P.C_ell_0       = 0.0;
P.C_ell_beta    = -0.12;
P.C_ell_p       = -0.26;
P.C_ell_r       = 0.14;
P.C_ell_delta_a = 0.08;
P.C_ell_delta_r = 0.105;

P.C_n_0         = 0.0;
P.C_n_beta      = 0.25;
P.C_n_p         = 0.022;
P.C_n_r         = -0.35;
P.C_n_delta_a   = 0.06;
P.C_n_delta_r   = -0.032;

% [p, q, r]
P.C_p_0 = Gamma_3*C_ell_0 + Gamma_4*C_n_0;
P.C_p_beta = Gamma_3*C_ell_beta + Gamma_4*C_n_beta;
P.C_p_p = Gamma_3*C_ell_p + Gamma_4*C_n_p;
P.C_p_r = Gamma_3*C_ell_r + Gamma_4*C_n_r;
P.C_p_delta_a = Gamma_3*C_ell_delta_a + Gamma_4*C_n_delta_a;
P.C_p_delta_r = Gamma_3*C_ell_delta_r + Gamma_4*C_n_delta_r;

P.C_r_0 = Gamma_4*C_ell_0 + Gamma_8*C_n_0;
P.C_r_beta = Gamma_4*C_ell_beta + Gamma_8*C_n_beta;
P.C_r_p = Gamma_4*C_ell_p + Gamma_8*C_n_p;
P.C_r_r = Gamma_4*C_ell_r + Gamma_8*C_n_r;
P.C_r_delta_a = Gamma_4*C_ell_delta_a + Gamma_8*C_n_delta_a;
P.C_r_delta_r = Gamma_4*C_ell_delta_r + Gamma_8*C_n_delta_r;

% [X, Y, Z]
P.C_X = -C_D*cos(alpha) + C_L*sin(alpha);
% P.C_X_0 = -P.C_D_0*cos(alpha) + P.C_L_0*sin(alpha);
% P.C_X_alpha = -P.C_D_alpha*cos(alpha) + P.C_L_alpha*sin(alpha);
P.C_X_q =  -C_D_q*cos(alpha) + C_L_q*sin(alpha);
P.C_X_delta_e =  -C_D_delta_e*cos(alpha) + C_L_delta_e*sin(alpha);

P.C_Z = -C_D*sin(alpha) - C_L*cos(alpha);
% P.C_Z_0= -P.C_D_0*sin(alpha) - P.C_L_0*cos(alpha);
% P.C_Z_alpha = -P.C_D_alpha*sin(alpha) + P.C_L_alpha*cos(alpha);
P.C_Z_q =  -C_D_q*sin(alpha) - C_L_q*cos(alpha);
P.C_Z_delta_e =  -C_D_delta_e*sin(alpha) - C_L_delta_e*cos(alpha);

P.C_Y_0         = 0.0;
P.C_Y_beta      = -0.98;
P.C_Y_p         = 0.0;
P.C_Y_r         = 0.0;
P.C_Y_delta_a   = 0.0;
P.C_Y_delta_r   = -0.17;

return
