% constants
g = 9.80665; % gravitational acceleration, P.mass/s^2
rho =  1.2133;  %kg/P.mass^3, density of air at 100 P.mass height (assume to be constant)
alpha=0;
% beta=0;
psi=0;

% vehicle specific parameters
P = param_uav_updated2(alpha);

%% LATERAL state-space, x = [v p r phi, psi]
Yv = (rho*P.S_wing*P.b*v/(4*P.mass*Va))*(P.C_Y_p*p + P.C_Y_r*r) ...
    + (rho*P.S_wing*v/P.mass)*(P.C_Y_0+P.C_Y_beta*beta*+P.C_Y_delta_a*delta_a+P.C_Y_delta_r*delta_r)...
    + (rho*P.S_wing*P.C_Y_beta/(2*P.mass))*norm([u, w]);
Yp = w + (rho*Va*P.S_wing*P.b/(4*P.mass))*P.C_Y_p;
Yr = -u + (rho*Va*P.S_wing*P.b/(4*P.mass))*P.C_Y_r;
Yda = (rho*Va^2*P.S_wing/(2*P.mass))*P.C_Y_delta_a;
Ydr = (rho*Va^2*P.S_wing/(2*P.mass))*P.C_Y_delta_r;

Lv = (rho*P.S_wing*P.b^2*v/(4*Va))*(P.C_p_p*p + P.C_p_r*r) ...
    + (rho*P.S_wing*P.b*v)*(P.C_p_0+P.C_p_beta*beta+P.C_p_delta_a*delta_a+P.C_p_delta_r*delta_r)...
    + (rho*P.S_wing*P.b*P.C_p_beta/2)*norm([u, w]);
Lp = P.Gamma_1*q + (rho*Va*P.S_wing*P.b^2/4)*P.C_p_p;
Lr = -P.Gamma_2*q + (rho*Va*P.S_wing*P.b^2/4)*P.C_p_r;
Lda = (rho*Va^2*P.S_wing*P.b/2)*P.C_p_delta_a;
Ldr = (rho*Va^2*P.S_wing*P.b/2)*P.C_p_delta_r;

Nv = (rho*P.S_wing*P.b^2*v/(4*Va))*(P.C_r_p*p + P.C_r_r*r) ...
    + (rho*P.S_wing*P.b*v)*(P.C_r_0+P.C_r_beta*beta+P.C_r_delta_a*delta_a+P.C_r_delta_r*delta_r)...
    + (rho*P.S_wing*P.b*P.C_r_beta)*norm([u, w]);
Np = P.Gamma_7*q + (rho*Va*P.S_wing*P.b^2/4)*P.C_r_p;
Nr = -P.Gamma_1*q + (rho*Va*P.S_wing*P.b^2/4)*P.C_r_r;
Nda = (rho*Va^2*P.S_wing*P.b/2)*P.C_r_delta_a;
Ndr = (rho*Va^2*P.S_wing*P.b/2)*P.C_r_delta_r;

%% LONGITUDINAL state space, x = [u w q theta h]
Xu = (u*rho*P.S_wing/P.mass)*(P.C_X_0+P.C_X_alpha*alpha+P.C_X_delta_e*delta_e) - (rho*P.S_wing*w*P.C_X_alpha/(2*P.mass)) ...
    + (rho*P.S_wing*P.c*P.C_X_q*u*q/(4*P.mass*Va)); %- (rho*Sprop*Cprop*u/P.mass); 
Xw = -q +  (w*rho*P.S_wing/P.mass)*(P.C_X_0+P.C_X_alpha*alpha+P.C_X_delta_e*delta_e) + (rho*P.S_wing*u*P.C_X_alpha/(2*P.mass)) ...
    + (rho*P.S_wing*P.c*P.C_X_q*w*q/(4*P.mass*Va)); %(rho*Sprop*Cprop*w/P.mass); 
Xq = -w + (rho*Va*P.S_wing*P.C_X_q*P.c)/(4*P.mass);
Xde = (rho*Va^2*P.S_wing*P.C_X_delta_e)/(2*P.mass);
% Xdt = rho*Sprop*Cprop*k^2*delta_t/P.mass;
Xdt = P.F_prop/P.mass;

Zu = q +  (u*rho*P.S_wing/P.mass)*(P.C_Z_0+P.C_Z_alpha*alpha+P.C_Z_delta_e*delta_e) - (rho*P.S_wing*w*P.C_Z_alpha/(2*P.mass)) ...
    + (u*rho*P.S_wing*P.c*P.C_Z_q*q/(4*P.mass*Va));  
Zw = (w*rho*P.S_wing/P.mass)*(P.C_Z_0+P.C_Z_alpha*alpha+P.C_Z_delta_e*delta_e) + (rho*P.S_wing*u*P.C_Z_alpha/(2*P.mass)) ...
    + (w*rho*P.S_wing*P.c*P.C_Z_q*q/(4*P.mass*Va));  
Zq = u + (rho*Va*P.S_wing*P.C_X_q*P.c)/(4*P.mass);
Zde = (rho*Va^2*P.S_wing*P.C_X_delta_e)/(2*P.mass);

Mu = (u*rho*P.S_wing*P.c/P.Jy)*(P.C_m_0+P.C_m_alpha*alpha+P.C_m_delta_e*delta_e) - (rho*P.S_wing*P.c*P.C_m_alpha*w)/(2*P.Jy) ...
    + (rho*P.S_wing*P.c^2*P.C_m_q*q*u)/(4*P.Jy*Va);
Mw = (w*rho*P.S_wing*P.c/P.Jy)*(P.C_m_0+P.C_m_alpha*alpha+P.C_m_delta_e*delta_e) + (rho*P.S_wing*P.c*P.C_m_alpha*u)/(2*P.Jy) ...
    + (rho*P.S_wing*P.c^2*P.C_m_q*q*w)/(4*P.Jy*Va);
Mq = (rho*Va*P.S_wing*P.c^2*P.C_m_q)/(4*P.Jy);
Mde = (rho*Va^2*P.S_wing*P.c*P.C_m_delta_e)/(4*P.Jy);
