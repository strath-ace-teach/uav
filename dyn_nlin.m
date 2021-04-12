function dyn_nlin(block)
%MSFUNTMPL_BASIC A Template for a Level-2 MATLAB S-Function
%   The MATLAB S-function is written as a MATLAB function with the
%   same name as the S-function. Replace 'msfuntmpl_basic' with the
%   name of your S-function.

%   Copyright 2003-2018 The MathWorks, Inc.

%
% The setup method is used to set up the basic attributes of the
% S-function such as ports, parameters, etc. Do not add any other
% calls to the main body of the function.
setup(block);

%endfunction

%% Function: setup ===================================================
% Abstract:
%   Set up the basic characteristics of the S-function block such as:
%   - Input ports
%   - Output ports
%   - Dialog parameters
%   - Options
%
%   Required         : Yes
%   C MEX counterpart: mdlInitializeSizes
%
function setup(block)

% Register number of ports
block.NumInputPorts  = 4; %controls
block.NumOutputPorts = 5; % pos(n,e), alt, head, vel

% Setup port properties to be inherited or dynamic
block.SetPreCompInpPortInfoToDynamic;
block.SetPreCompOutPortInfoToDynamic;

% Override input port properties
for i=1:block.NumInputPorts
block.InputPort(i).Dimensions        = 1;
block.InputPort(i).DatatypeID  = 0;  % double
block.InputPort(i).Complexity  = 'Real';
block.InputPort(i).DirectFeedthrough = true;
block.InputPort(i).SamplingMode = 'Sample';
% SetInputPortSamplingMode
end

% Override output port properties
for i=1:block.NumOutputPorts
block.OutputPort(i).Dimensions       = 1;
block.OutputPort(i).DatatypeID  = 0; % double
block.OutputPort(i).Complexity  = 'Real';
block.OutputPort(i).SamplingMode = 'Sample';
end

% Register parameters
block.NumDialogPrms     = 0;

% Set up the continuous states.
block.NumContStates = 12;

% Register sample times
%  [0 offset]            : Continuous sample time
%  [positive_num offset] : Discrete sample time
%
%  [-1, 0]               : Inherited sample time
%  [-2, 0]               : Variable sample time
block.SampleTimes = [0 0];


% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'CustomSimState',  < Has GetSimState and SetSimState methods
%    'DisallowSimState' < Error out when saving or restoring the model sim state
block.SimStateCompliance = 'DefaultSimState';

%% -----------------------------------------------------------------
% The MATLAB S-function uses an internal registry for all
% block methods. You should register all relevant methods
% (optional and required) as illustrated below. You may choose
% any suitable name for the methods and implement these methods
% as local functions within the same file. See comments
% provided for each function for more information.
%% -----------------------------------------------------------------

% block.RegBlockMethod('PostPropagationSetup',    @DoPostPropSetup);
block.RegBlockMethod('InitializeConditions', @InitializeConditions);
% block.RegBlockMethod('Start', @Start);
block.RegBlockMethod('Outputs', @Outputs);     % Required
% block.RegBlockMethod('Update', @Update);
block.RegBlockMethod('Derivatives', @Derivatives);
block.RegBlockMethod('Terminate', @Terminate); % Required

%end setup


%% PostPropagationSetup:
%   Functionality    : Setup work areas and state variables. Can
%                      also register run-time methods here
%   Required         : No
%   C MEX counterpart: mdlSetWorkWidths
%
function DoPostPropSetup(block)
% block.NumDworks = 1;

%   block.Dwork(1).Name            = 'x1';
%   block.Dwork(1).Dimensions      = 1;
%   block.Dwork(1).DatatypeID      = 0;      % double
%   block.Dwork(1).Complexity      = 'Real'; % real
%   block.Dwork(1).UsedAsDiscState = false; %true;


%% InitializeConditions:
%   Functionality    : Called at the start of simulation and if it is
%                      present in an enabled subsystem configured to reset
%                      states, it will be called when the enabled subsystem
%                      restarts execution to reset the states.
%   Required         : No
%   C MEX counterpart: mdlInitializeConditions
%
function InitializeConditions(block)
% block.ContStates.Data = block.DialogPrm(1).Data;

block.ContStates.Data(1) = 0; %n;
block.ContStates.Data(2) = 0; %e;
block.ContStates.Data(3) = 100; %h;

block.ContStates.Data(4) = 30; % u;
block.ContStates.Data(5) = 20; %v;
block.ContStates.Data(6) = 0; %w;

block.ContStates.Data(7) = 0; %phi;
block.ContStates.Data(8) = 0; %theta;
block.ContStates.Data(9) = 0; %psi;

block.ContStates.Data(10) = 0; %p;
block.ContStates.Data(11) = 0; %q;
block.ContStates.Data(12) = 0; %r;

%end InitializeConditions


%% Start:
%   Functionality    : Called once at start of model execution. If you
%                      have states that should be initialized once, this
%                      is the place to do it.
%   Required         : No
%   C MEX counterpart: mdlStart
%
 function Start(block)
% 
% block.Dwork(1).Data = 0;

%end Start

%% Outputs:
%   Functionality    : Called to generate block outputs in
%                      simulation step
%   Required         : Yes
%   C MEX counterpart: mdlOutputs
%
function Outputs(block)

block.OutputPort(1).Data = block.ContStates.Data(1); %pos-n
block.OutputPort(2).Data = block.ContStates.Data(2); %pos-e
block.OutputPort(3).Data = block.ContStates.Data(3); %pos-h
block.OutputPort(4).Data = norm([block.ContStates.Data(4); block.ContStates.Data(5); block.ContStates.Data(6)]); % vel

u = block.ContStates.Data(4); % u;
v= block.ContStates.Data(5); %v;
w= block.ContStates.Data(6); %w;
phi=block.ContStates.Data(7); %phi;
theta=block.ContStates.Data(8); %theta;
psi=block.ContStates.Data(9); %psi;
vn = u*cos(theta)*cos(psi) + v*(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(phi)) + w*(cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi));
ve = u*cos(theta)*sin(psi) + v*(sin(phi)*sin(theta)*sin(psi)+cos(phi)*sin(phi)) + w*(cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi));
block.OutputPort(5).Data = atan2(ve,vn); % heading

%end Outputs

%% Derivatives:
%   Functionality    : Called to update derivatives of
%                      continuous states during simulation step
%   Required         : No
%   C MEX counterpart: mdlDerivatives
function Derivatives(block)

delta_t = block.InputPort(1).Data;
delta_e = block.InputPort(2).Data;
delta_r = block.InputPort(3).Data;
delta_a = block.InputPort(4).Data;

n = block.ContStates.Data(1); %n;
e = block.ContStates.Data(2); %e;
h = block.ContStates.Data(3); %h;

u = block.ContStates.Data(4); % u;
v= block.ContStates.Data(5); %v;
w= block.ContStates.Data(6); %w;

phi=block.ContStates.Data(7); %phi;
theta=block.ContStates.Data(8); %theta;
psi=block.ContStates.Data(9); %psi;

p=block.ContStates.Data(10); %p;
q=block.ContStates.Data(11); %q;
r=block.ContStates.Data(12); %r;


%%
Va = sqrt(u^2+v^2+w^2);
alpha = atan2(w,u);
beta  = asin(v/Va);
g = 9.80665; % gravitational acceleration, m/s^2
rho =  1.2133;  %kg/P.mass^3, density of air at 100 m height (assume to be constant)

%% Physical parameters of airframe
P.mass = 13.5; %mass, kg
P.S_wing        = 0.55;  %wing area, m^2
P.b             = 2.8956;  %wing span, m
P.c             = 0.18994;  %, mean aerodynamic chord, m
P.e             = 0.9;  %Oswald efficiency
P.F_prop = 1.2*P.mass*g;  %max thrust from propeller (** specific to this assignment)
P.AR = (P.b^2)/P.S_wing;  % aspect ratio

% moments of inertia (matrix J)
P.Jx   = 0.8244;  %kg-m^2
P.Jy   = 1.135; %kg-m^2
P.Jz   = 1.759; %kg-m^2
P.Jxz  = 0.1204; %kg-m^2

Gamma = P.Jx*P.Jz - P.Jxz^2;
P.Gamma_1 = (P.Jxz*(P.Jx-P.Jy+P.Jz))/Gamma;
P.Gamma_2 = (P.Jz*(P.Jz-P.Jy)+P.Jxz^2)/Gamma;
P.Gamma_3 = P.Jz/Gamma;
P.Gamma_4 = P.Jxz/Gamma;
P.Gamma_5 = (P.Jz - P.Jx)/P.Jy;
P.Gamma_6 = P.Jxz/P.Jy;
P.Gamma_7 = ((P.Jx-P.Jy)*P.Jx + P.Jxz^2)/Gamma;
P.Gamma_8 = P.Jx/Gamma;

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
P.C_L = P.C_L_0 + P.C_L_alpha*alpha;
P.C_D = P.C_D_0 + P.C_D_alpha*alpha;

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
P.C_p_0 = P.Gamma_3*P.C_ell_0 + P.Gamma_4*P.C_n_0;
P.C_p_beta = P.Gamma_3*P.C_ell_beta + P.Gamma_4*P.C_n_beta;
P.C_p_p = P.Gamma_3*P.C_ell_p + P.Gamma_4*P.C_n_p;
P.C_p_r = P.Gamma_3*P.C_ell_r + P.Gamma_4*P.C_n_r;
P.C_p_delta_a = P.Gamma_3*P.C_ell_delta_a + P.Gamma_4*P.C_n_delta_a;
P.C_p_delta_r = P.Gamma_3*P.C_ell_delta_r + P.Gamma_4*P.C_n_delta_r;

P.C_r_0 = P.Gamma_4*P.C_ell_0 + P.Gamma_8*P.C_n_0;
P.C_r_beta = P.Gamma_4*P.C_ell_beta + P.Gamma_8*P.C_n_beta;
P.C_r_p = P.Gamma_4*P.C_ell_p + P.Gamma_8*P.C_n_p;
P.C_r_r = P.Gamma_4*P.C_ell_r + P.Gamma_8*P.C_n_r;
P.C_r_delta_a = P.Gamma_4*P.C_ell_delta_a + P.Gamma_8*P.C_n_delta_a;
P.C_r_delta_r = P.Gamma_4*P.C_ell_delta_r + P.Gamma_8*P.C_n_delta_r;

% [X, Y, Z]
P.C_X = -P.C_D*cos(alpha) + P.C_L*sin(alpha);
% P.C_X_0 = -P.C_D_0*cos(alpha) + P.C_L_0*sin(alpha);
% P.C_X_alpha = -P.C_D_alpha*cos(alpha) + P.C_L_alpha*sin(alpha);
P.C_X_q =  -P.C_D_q*cos(alpha) + P.C_L_q*sin(alpha);
P.C_X_delta_e =  -P.C_D_delta_e*cos(alpha) + P.C_L_delta_e*sin(alpha);

P.C_Z = -P.C_D*sin(alpha) - P.C_L*cos(alpha);
% P.C_Z_0= -P.C_D_0*sin(alpha) - P.C_L_0*cos(alpha);
% P.C_Z_alpha = -P.C_D_alpha*sin(alpha) + P.C_L_alpha*cos(alpha);
P.C_Z_q =  -P.C_D_q*sin(alpha) - P.C_L_q*cos(alpha);
P.C_Z_delta_e =  -P.C_D_delta_e*sin(alpha) - P.C_L_delta_e*cos(alpha);

P.C_Y_0         = 0.0;
P.C_Y_beta      = -0.98;
P.C_Y_p         = 0.0;
P.C_Y_r         = 0.0;
P.C_Y_delta_a   = 0.0;
P.C_Y_delta_r   = -0.17;

% Coupled non-linear equations of motion: (d/dt) of position (x3) + velocity (x3), attitude (x3) + angular rotation (x3)
ndot = u*cos(theta)*cos(psi) + v*(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(phi)) + w*(cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi));
edot = u*cos(theta)*sin(psi) + v*(sin(phi)*sin(theta)*sin(psi)+cos(phi)*sin(phi)) + w*(cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi));
hdot = u*sin(theta) - v*sin(phi)*cos(theta) - w*cos(phi)*cos(theta);  % h = -pz

pdot = P.Gamma_1*p*q - P.Gamma_2*q*r + (rho*P.b^2*Va*P.S_wing/4)*(P.C_p_p*p + P.C_p_r*r) ...
    + (P.S_wing*P.b*rho*Va^2/2)*(P.C_p_0 + P.C_p_beta*beta + P.C_p_delta_a*delta_a + P.C_p_delta_r*delta_r);
qdot = P.Gamma_5*p*r - P.Gamma_6*(p^2-r^2) + (rho*Va^2*P.S_wing*P.c/(2*P.Jy))*(P.C_m_0 + P.C_m_alpha*alpha + P.C_m_q*P.c*q/(2*Va) + P.C_m_delta_e*delta_e);
rdot = P.Gamma_7*p*q - P.Gamma_1*q*r + (rho*P.b^2*Va*P.S_wing/4)*(P.C_r_p*p + P.C_r_r*r) ...
    + (P.S_wing*P.b*rho*Va^2/2)*( P.C_r_0 + P.C_r_beta*beta + P.C_r_delta_a*delta_a + P.C_r_delta_r*delta_r);

thetadot = q*cos(phi) - r*sin(phi);
phidot = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
psidot = q*sin(phi)*sec(theta) + r*cos(phi)*sec(theta);

udot = r*v - q*w - g*sin(theta) + (rho*Va^2*P.S_wing/(2*P.mass))*(P.C_X + P.C_X_q*P.c*q/(2*Va) + P.C_X_delta_e*delta_e) ...
    + (P.F_prop/P.mass)*delta_t;   
vdot = p*w - r*u + g*cos(theta)*sin(phi) + (rho*P.S_wing*Va*P.b/(4*P.mass))*(P.C_Y_p*p + P.C_Y_r*r) ...
    + (rho*P.S_wing*Va^2/(2*P.mass))*(P.C_Y_0 + P.C_Y_beta*beta + P.C_Y_delta_a*delta_a + P.C_Y_delta_r*delta_r);
wdot = q*u - p*v + g*cos(theta)*cos(phi) + (rho*Va^2*P.S_wing/(2*P.mass))*(P.C_Z + P.C_Z_q*P.c*q/(2*Va) + P.C_Z_delta_e*delta_e);

block.Derivatives.Data(1) = ndot;
block.Derivatives.Data(2) = edot;
block.Derivatives.Data(3) = hdot;

block.Derivatives.Data(4) = udot;
block.Derivatives.Data(5) = vdot;
block.Derivatives.Data(6) = wdot;

block.Derivatives.Data(7) = phidot;
block.Derivatives.Data(8) = thetadot;
block.Derivatives.Data(9) = psidot;

block.Derivatives.Data(10) = pdot;
block.Derivatives.Data(11) = qdot;
block.Derivatives.Data(12) = rdot;

%end Derivatives

%%
%% Terminate:
%   Functionality    : Called at the end of simulation for cleanup
%   Required         : Yes
%   C MEX counterpart: mdlTerminate
%
 function Terminate(block)

%end Terminate

