%% ME528 Control Systems Design
%{
Code fragment listing the equations of motion for the Group Assignment
You can change variables names as appropriate for your code. 
%}

alpha=0.15;
% Coupled non-linear equations of motion: (d/dt) of position (x3) + velocity (x3), attitude (x3) + angular rotation (x3)
ndot = u*cos(theta)*cos(psi) + v*(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(phi)) + w*(cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi));
edot = u*cos(theta)*sin(psi) + v*(sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(phi)) + w*(cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi));
hdot = u*sin(theta) - v*sin(phi)*cos(theta) - w*cos(phi)*cos(theta);  % h = -pz

pdot = P.Gamma_1*p*q - P.Gamma_2*q*r + (rho*P.b^2*Va*P.S_wing/4)*(P.C_p_p*p + P.C_p_r*r) ...
    + (P.S_wing*P.b*rho*Va^2/2)*(P.C_p_0 + P.C_p_beta*beta + P.C_p_delta_a*delta_a + P.C_p_delta_r*delta_r);

qdot = P.Gamma_5*p*r - P.Gamma_6*(p^2-r^2) + (rho*Va^2*P.S_wing*P.c/(2*P.Jy))*(P.C_m_0 + P.C_m_alpha*alpha + P.C_m_q*P.c*q/(2*Va) + P.C_m_delta_e*delta_e);

rdot = P.Gamma_7*p*q - P.Gamma_1*q*r + (rho*P.b^2*Va*P.S_wing/4)*(P.C_r_p*p + P.C_r_r*r) ...
    + (P.S_wing*P.b*rho*Va^2/2)*( P.C_r_0 + P.C_r_beta*beta + P.C_r_delta_a*delta_a + P.C_r_delta_r*delta_r);

thetadot = q*cos(phi) - r*sin(phi);
phidot = p + q*sin(phi)*tan(theta) + r*cos(phi)*tan(theta);
psidot = q*sin(phi)*sec(theta) + r*cos(phi)*sec(theta); %when wind=0, psidot=chidot

udot = r*v - q*w - g*sin(theta) + (rho*Va^2*P.S_wing/(2*P.mass))*(P.C_X + P.C_X_q*P.c*q/(2*Va) + P.C_X_delta_e*delta_e) ...
    + (P.F_prop/P.mass)*delta_t;   
vdot = p*w - r*u + g*cos(theta)*sin(phi) + (rho*P.S_wing*Va*P.b/(4*P.mass))*(P.C_Y_p*p + P.C_Y_r*r) ...
    + (rho*P.S_wing*Va^2/(2*P.mass))*(P.C_Y_0 + P.C_Y_beta*beta + P.C_Y_delta_a*delta_a + P.C_Y_delta_r*delta_r);
wdot = q*u - p*v + g*cos(theta)*cos(phi) + (rho*Va^2*P.S_wing/(2*P.mass))*(P.C_Z + P.C_Z_q*P.c*q/(2*Va) + P.C_Z_delta_e*delta_e);