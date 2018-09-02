x_St = [ 22.22; 0; 0]; %Anfangszustand x_k0
x_dot_St = zeors(3,1); % Anfangszustand der Ableitung von x_k0
x_vec_ist_k = [22; 2; 3]; % x_k Zustand
x_dot_ist_k = [-3; 0; 1]; % Ableitung von x_k
a_x0 = 0; % Anfangsbeschleunigung;
u = -ones(4,1); % aktuelle Bremskraft vom MPC opt
d_St = zeros(4,1); % Anfangsseitenkräfte
d = [20; 30; -13; -34]; % aktuelle Seitenkräfte
c_w = 0.5; % Luftwiderstandsbeiwert;
A_L = 4; % LKW Querschnittsfläche von vorne
rho_L = 1.24; % Luftdichte;
m_Fzg = 6800; % Fahrzeugmasse
delta_rel = zeros(4,1); %Lenkwinkel pro Rad
n_L = ones(4,1); % Reifennachlauf in Längsrichtung
n_L0 = zeros(4,1);
n_S = ones(4,1); % Reifennachlauf in Querrichtung
n_S0 = zeros(4,1);
l_v = 1.5; % Radstand vorne
l_h = 1.52; % Radstand hinten
B_v = 2.3; % Fahrzeugbreite vorne
B_h = 2.2; %Fahrzeugbreite hinten
t_sim = tic; %Simulationslaufzeit
T_ab = 0.001; % Abtastzeit;
c_yv = 45000; % Reifensteifigkeit vorne(links und rechts sind gleich)
c_yh = 48000; % Reifensteifigkeit hinten (links und rechts sind gleich)

J_z = 0.5*m_Fzg*((l_v+l_h)/2)^2; % [kgm^2] Massenträgheit

% Lenkwinkel am Rad
delta_rel_VL = delta_rel(1);
delta_rel_VR = delta_rel(2);
delta_rel_HL = delta_rel(3);
delta_rel_HR = delta_rel(4);

 % Anfangs- Bremskräfte
F_B_VL_0 = a_x0*(m_Fzg/4);
F_B_VR_0 = a_x0*(m_Fzg/4);
F_B_HL_0 = a_x0*(m_Fzg/4);
F_B_HR_0 = a_x0*(m_Fzg/4);

u_St = [F_B_VL_0; F_B_VR_0; F_B_HL_0; F_B_HR_0];

v_0 = x_St(1);
beta_0 = x_St(2);
psi_dot_0 = x_St(3);

% Anfangs-Seitenkräfte am Rad
F_S_VL_0 = d_St(1);
F_S_VR_0 = d_St(2);
F_S_HL_0 = d_St(3);
F_S_HR_0 = d_St(4);
    
% Kräfte in x-Richtung pro Rad
F_x_VL_0 = (F_B_VL_0*cosd(delta_rel_VL))-(F_S_VL_0*sind(delta_rel_VL));
F_x_VR_0 = (F_B_VR_0*cosd(delta_rel_VR))-(F_S_VR_0*sind(delta_rel_VR));
F_x_HL_0 = (F_B_HL_0*cosd(delta_rel_HL))-(F_S_HL_0*sind(delta_rel_HL));
F_x_HR_0 = (F_B_HR_0*cosd(delta_rel_HR))-(F_S_HR_0*sind(delta_rel_HR));
    
% Kräfte in y-Richtung pro Rad
F_y_VL_0 = (F_S_VL_0*cosd(delta_rel_VL))+(F_B_VL_0*sind(delta_rel_VL));
F_y_VR_0 = (F_S_VR_0*cosd(delta_rel_VR))+(F_B_VR_0*sind(delta_rel_VR));
F_y_HL_0 = (F_S_HL_0*cosd(delta_rel_HL))+(F_B_HL_0*sind(delta_rel_HL));
F_y_HR_0 = (F_S_HR_0*cosd(delta_rel_HR))+(F_B_HR_0*sind(delta_rel_HR));
    
F_x0 = [F_x_VL_0; F_x_VR_0; F_x_HL_0; F_x_HR_0];
F_y0 = [F_y_VL_0; F_y_VR_0; F_y_HL_0; F_y_HR_0];

n_LV0 = (n_L0(1)+n_L0(2))/2;
n_LH0 = (n_L0(3)+n_L0(4))/2;
n_S_VL0 = n_S0(1);
n_S_VR0 = n_S0(2);

A110 = -(cosd(beta_0)*c_w*A_L*rho_L*v_0)/(2*m_Fzg);
A120 = (-(sind(beta_0)/m_Fzg)*((F_B_VL_0*cosd(delta_rel_VL))+(F_B_VR_0*cosd(delta_rel_VR))-F_S_VL_0*(sind(delta_rel_VL))+(F_S_VR_0*sind(delta_rel_VR))+F_B_HL_0+F_B_HR_0-((c_w*A_L*rho_L*v_0^2)/2)))+...
    ((cosd(beta_0)/m_Fzg)*((F_S_VL_0*(cosd(delta_rel_VL))+(F_S_VR_0*cosd(delta_rel_VR))+(F_B_VL_0*sind(delta_rel_VL))+(F_B_VR_0*sind(delta_rel_VR))+F_S_HL_0+F_S_HR_0)));
A130 = 0;
A210 = (-(cosd(beta_0)/(m_Fzg*v_0^2))*((F_S_VL_0*cosd(delta_rel_VL))+(F_B_VL_0*sind(delta_rel_VL))+(F_S_VR_0*cosd(delta_rel_VR))+(F_B_VR_0*sind(delta_rel_VR))+F_S_HL_0+F_S_HR_0))+...
    ((sind(beta_0)/(m_Fzg*v_0^2))*((F_B_VL_0*cosd(delta_rel_VL))-(F_S_VL_0*sind(delta_rel_VL))+(F_B_VR_0*cosd(delta_rel_VR))-(F_S_VR_0*sind(delta_rel_VR))+F_B_HL_0+F_B_HR_0))+...
    ((sind(beta_0)*c_w*A_L*rho_L)/(2*m_Fzg));
A220 = ((-sind(beta_0)/(m_Fzg*v_0^2))*((F_S_VL_0*cosd(delta_rel_VL))+(F_B_VL_0*sind(delta_rel_VL))+(F_S_VR_0*cosd(delta_rel_VR))+(F_B_VR_0*sind(delta_rel_VR))+F_S_HL_0+F_S_HR_0))-...
    ((cosd(beta_0)/(m_Fzg*v_0^2))*((F_B_VL_0*cosd(delta_rel_VL))-(F_S_VL_0*sind(delta_rel_VL))+(F_B_VR_0*cosd(delta_rel_VR))-(F_S_VR_0*sind(delta_rel_VR))+F_B_HL_0+F_B_HR_0-((c_w*A_L*rho_L*v_0^2)/(2))));
A230 = -1;
A310 = 0;
A320 = 0;
A330 = 0;

A0 = [ A110 A120 A130; 
    A210 A220 A230; 
    A310 A320 A330];

B110 = ((cosd(beta_0)*cosd(delta_rel_VL))+(sind(beta_0)*sind(delta_rel_VL)))/(m_Fzg);
B120 = ((cosd(beta_0)*cosd(delta_rel_VR))+(sind(beta_0)*sind(delta_rel_VR)))/(m_Fzg);
B130 = (cosd(beta_0))/m_Fzg;
B140 = (cosd(beta_0))/m_Fzg;
B210 = ((cosd(beta_0)*sind(delta_rel_VL))-(sind(beta_0)*cosd(delta_rel_VL)))/(m_Fzg*v_0);
B220 = ((cosd(beta_0)*sind(delta_rel_VR))-(sind(beta_0)*cosd(delta_rel_VR)))/(m_Fzg*v_0);
B230 = (-sind(beta_0))/(m_Fzg*v_0);
B240 = (-sind(beta_0))/(m_Fzg*v_0);
B310 = (1/J_z)*((sind(delta_rel_VL)*(l_v-n_LV0))-((B_v*cosd(delta_rel_VL))/2)-(n_S_VL0*sind(delta_rel_VL)*cosd(delta_rel_VL)));
B320 = (1/J_z)*((sind(delta_rel_VR)*(l_v*n_LV0))+((B_v*cosd(delta_rel_VR))/2)-(n_S_VR0*sind(delta_rel_VR)*cosd(delta_rel_VR)));
B330 = -((B_h)/(2*J_z));
B340 = (B_h)/(2*J_z);

B0 = [ B110 B120 B130 B140;
    B210 B220 B230 B240;
    B310 B320 B330 B340];

E110 = ((sind(beta_0)*cosd(delta_rel_VL))-(cosd(beta_0)*sind(delta_rel_VL)))/(m_Fzg);
E120 = ((sind(beta_0)*cosd(delta_rel_VR))-(cosd(beta_0)*sind(delta_rel_VR)))/(m_Fzg);
E130 = (sind(beta_0))/(m_Fzg*v_0);
E140 = (sind(beta_0))/(m_Fzg*v_0);
E210 = ((cosd(beta_0)*cosd(delta_rel_VL))+(sind(beta_0)*sind(delta_rel_VR)))/(m_Fzg*v_0);
E220 = ((cosd(beta_0)*cosd(delta_rel_VR))+(sind(beta_0)*sind(delta_rel_VR)))/(m_Fzg*v_0);
E230 = (cosd(beta_0))/(m_Fzg*v_0);
E240 = (cosd(beta_0))/(m_Fzg*v_0);
E310 = (1/J_z)*((cosd(delta_rel_VL)*(l_v*n_LV0))+((B_v*sind(delta_rel_VL))/2)+(n_S_VL0*(sind(delta_rel_VL))^2));
E320 = (1/J_z)*((cosd(delta_rel_VR)*(l_v*n_LV0))+((B_v*sind(delta_rel_VR))/2)+(n_S_VR0*(sind(delta_rel_VR))^2));
E330 = -(l_h+n_LH0);
E340 = -(l_h+n_LH0);

E0 = [ E110 E120 E130 E140;
    E210 E220 E230 E240;
    E310 E320 E330 E340];

C0 = [1 0 0;
    0 0 1];

D0 = [0 0 0 0;
    0 0 0 0;
    0 0 0 0];


if t_sim == 0
    A = A0;
    B = B0;
    E = E0;
    C = C0;
    u = u_St;
    d = d_St; 
    x_vec_ist_k = x_St;
    x_dot_ist_k = x_dot_St;
else
    % Berechnung der werte bei k>0   
    
    v_Fzg = x_vec_ist_k(1);
    beta = x_vec_ist_k(2);
    psi_dot = x_vec_ist_k(3);

    F_B_VL = u(1);
    F_B_VR = u(2);
    F_B_HL = u(3);
    F_B_HR = u(4);

    F_S_VL = d(1);
    F_S_VR = d(2);
    F_S_HL = d(3);
    F_S_HR = d(4);

    n_LV = (n_L(1)+n_L(2))/2;
    n_LH = (n_L(3)+n_L(4))/2;
    n_S_VL = n_S(1);
    n_S_VR = n_S(2);

    A11 = -((c_w*A_L*rho_L*v_Fzg*cosd(beta))/(m_Fzg))...
        +((c_yv*l_v*psi_dot*sind(beta-delta_rel_VL))/(m_Fzg*v_Fzg^2))...
        -((c_yv*l_v*psi_dot*sind(delta_rel_VR-beta))/(m_Fzg*v_Fzg^2))...
        -((sind(beta)*c_yh*l_h*psi_dot)/(v_Fzg^2))...
        -((sind(beta)*c_yh*l_h*psi_dot)/(v_Fzg^2));
    
    A12 = -((F_B_VL*sind(beta-delta_rel_VL)*F_B_VR*cosd(delta_rel_VR-beta))/(m_Fzg))...
        -((F_B_VL*cosd(beta-delta_rel_VL)*F_B_VR*sind(delta_rel_VR-beta))/(m_Fzg))...
        -(((F_B_HL+F_B_HR-(c_w*A_L*(rho_L/2)*v_Fzg^2))*sind(beta))/(m_Fzg))...
        +((c_yv*delta_rel_VL*sind(beta-delta_rel_VL))/(m_Fzg))...
        -((c_yv*sind(beta-delta_rel_VR))/(m_Fzg))...
        -((c_yv*beta*cosd(beta-delta_rel_VL))/(m_Fzg))...
        -((c_yv*l_v*psi_dot*sind(beta_delta_rel_VL))/(m_Fzg*v_Fzg))...
        +((c_yv*sind(delta_rel_VR-beta))/(m_Fzg))...
        +((c_yv*beta*cosd(delta_rel_VR-beta))/(m_Fzg))...
        -((c_yv*delta_rel_VR*cosd(delta_rel_VR-beta))/(m_Fzg))...
        +((c_yv*l_v*psi_dot*cosd(delta_rel_VR-beta))/(m_Fzg*v_Fzg))...
        +((c_yh*l_h*psi_dot*cosd(beta))/(m_Fzg*v_Fzg))...
        -((c_yh*sind(beta))/(m_Fzg))...
        -((c_yh*beta*cosd(beta))/(m_Fzg))...
        -((c_yh*sind(beta))/(m_Fzg))...
        -((c_yh*beta*cosd(beta))/(m_Fzg))...
        +((c_yh*l_h*psi_dot*cosd(beta))/(v_Fzg*m_Fzg));
    
    A13 = -((c_yv*l_v*sind(beta-delta_rel_VL))/(v_Fzg*m_Fzg))...
        +((c_yv*l_v*sind(delta_rel_VR-beta))/(v_Fzg*m_Fzg))...
        +((sind(beta)*c_yh*l_h)/(m_Fzg*v_Fzg))...
        +((sind(beta)*c_yh*l_h)/(m_Fzg*v_Fzg));
    
    A21 = ((F_B_VL*sind(beta-delta_rel_VL))/(m_Fzg*v_Fzg^2))...
        -((F_B_VR*sind(delta_rel_VR-beta))/(m_Fzg*v_Fzg^2))...
        -((c_yv*delta_rel_VK*cosd(beta-delta_rel_VL))/(m_Fzg*v_Fzg^2))...
        +((c_yv*beta*cosd(beta-delta_rel_VL))/(m_Fzg*v_Fzg^2))...
        +((2*c_yv*l_v*psi_dot*cosd(beta-delta_rel_VL))/(m_Fzg*v_Fzg^3))...
        -((c_yv*delta_rel_VR*cosd(delta_rel_VR-beta))/(m_Fzg*v_Fzg^2))...
        +((c_yv*beta*cosd(delta_rel_VR-beta))/(m_Fzg*v_Fzg^2))...
        +((2*c_yv*l_v*psi_dot*cosd(delta_rel_VR-beta))/(m_Fzg*v_Fzg^3))...
        -((2*cosd(beta)*c_yh*l_h*psi_dot)/(m_Fzg*v_Fzg^3))...
        +((c_yh*beta*cosd(beta))/(m_Fzg*v_Fzg^2))...
        -((2*c_yh*l_h*psi_dot*cosd(beta))/(m_Fzg*v_Fzg^3))...
        +((c_yh*beta*cosd(beta))/(m_Fzg*v_Fzg^2))...
        +((sind(beta)*c_w*A_L*rho_L)/(2*m_Fzg))...
        +((sind(beta)*F_B_HL)/(m_Fzg*v_Fzg^2))...
        +((sind(beta)*F_B_HR)/(m_Fzg*v_Fzg^2));
    
    A22 = -((F_B_VL*cosd(beta_delta_rel_VL))/(m_Fzg*v_Fzg))...
        +((F_B_VR*cosd(delta_rel_VR-beta))/(m_Fzg*v_Fzg))...
        -((c_yv*delta_rel_VL*sind(beta-delta_rel_VL))/(m_Fzg*v_Fzg))...
        +((c_yv*beta*sind(beta-delta_rel_VL))/(m_Fzg*v_Fzg))...
        -((c_yv*cosd(beta-delta_rel_VL))/(m_Fzg*v_Fzg))...
        +((c_yv*l_v*psi_dot*sind(beta-delta_rel_VL))/(m_Fzg*v_Fzg^2))...
        +((c_yv*delta_rel_VR*sind(delta_rel_VR-beta))/(m_Fzg*v_Fzg))...
        -((c_yv*cosd(delta_rel_VR-beta))/(m_Fzg*v_Fzg))...
        -((c_yh*l_h*psi_dot*sind(beta))/(m_Fzg*v_Fzg^2))...
        -((c_yh*beta*sind(beta))/(m_Fzg*v_Fzg))...
        -((c_yh*cosd(beta))/(m_Fzg*v_Fzg))...
        +((c_yh*l_h*psi_dot*sind(beta))/(m_Fzg*v_Fzg^2))...
        -((c_yh*cosd(beta))/(m_Fzg*v_Fzg))...
        -((c_yh*beta*sind(beta))/(m_Fzg*v_Fzg))...
        -((c_yh*cosd(beta))/(m_Fzg*v_Fzg))...
        +((c_yh*beta*sind(beta))/(m_Fzg*v_Fzg))...
        +((c_w*A_L*rho_L*v_Fzg*cosd(beta))/(2*m_Fzg))...
        -((F_B_HL*cosd(beta))/(m_Fzg*v_Fzg))...
        -((F_B_HR*cosd(beta))/(m_Fzg*v_Fzg));
    
    A23 = -((c_yv*l_v*cosd(beta-delta_rel_VL))/(m_Fzg*v_Fzg^2))...
        -((c_yv*l_v*cosd(delta_rel_VR-beta))/(m_Fzg*v_Fzg^2))...
        +((c_yh*l_h*cosd(beta))/(m_Fzg*v_Fzg^2))...
        +((c_yh*l_h*cosd(beta))/(m_Fzg*v_Fzg^2))-1;
    
    A31 = ((c_yv*l_v*psi_dot*B_v*sind(delta_rel_VL))/(2*J_z*v_Fzg^2))...
        +((2*c_yv*(l_v^2)*psi_dot*cosd(delta_rel_VL))/(2*J_z*v_Fzg^2))...
        -((c_yv*l_v*psi_dot*B_v*sind(delta_rel_VR))/(2*J_z*v_Fzg^2))...
        +((2*c_yv*(l_v^2)*psi_dot*cosd(delta_rel_VR))/(2*J_z*v_Fzg^2))...
        -((((2*c_yv)+(2*c_yv))*l_v*psi_dot*n_LV)/(2*J_z*v_Fzg^2))...
        -(((-(2*cyh*l_h)-(2*c_yh*n_LH)-(2*c_yh*n_LR)-(2*c_yh*l_h))*l_h*psi_dot)/(2*J_z*v_Fzg^2));

   A32 = ((2*c_yv*n_LV)+(2*c_yv*n_LV)+(2*c_yh*n_LH)+(2*c_yh*l_h)+(2*c_yh*n_LH)+(2*c_yh*l_h))/(2*J_z);
   
   A33 = -((c_yv*l_v*B_v*sind(delta_rel_VL))/(2*J_z*v_Fzg))...
       -((c_yv*(l_v^2)*cosd(delta_rel_VL))/(J_z*v_Fzg))...
       +((c_yv*l_v*B_v*sind(delta_rel_VR))/(2*J_z*v_Fzg))...
       -((2*c_yv*(l_v^2)*cosd(delta-rel_VR))/(2*J_z*v_Fzg))...
       +((((2*c_yv)+(2*c_yv))*l_v*n_LV)/(2*J_z*v_Fzg))...
       +((-(2*c_yh*l_h)-(2*c_yh*n_LH)-(2*c_yh*n_LR)-(2*c_yh*l_h))*l_h/(2*J_z*v_Fzg));
   
    A = [ A11 A12 A13; 
        A21 A22 A23; 
        A31 A32 A33];
    
% 
%     B11 = ((cosd(beta)*cosd(delta_rel_VL))+(sind(beta)*sind(delta_rel_VL)))/(m_Fzg);
%     B12 = ((cosd(beta)*cosd(delta_rel_VR))+(sind(beta)*sind(delta_rel_VR)))/(m_Fzg);
%     B13 = (cosd(beta))/m_Fzg;
%     B14 = (cosd(beta))/m_Fzg;
%     B21 = ((cosd(beta)*sind(delta_rel_VL))-(sind(beta)*cosd(delta_rel_VL)))/(m_Fzg*v_Fzg);
%     B22 = ((cosd(beta)*sind(delta_rel_VR))-(sind(beta)*cosd(delta_rel_VR)))/(m_Fzg*v_Fzg);
%     B23 = (-sind(beta))/(m_Fzg*v_Fzg);
%     B24 = (-sind(beta))/(m_Fzg*v_Fzg);
%     B31 = (psi_dot/M_Gier_ist)*((sind(delta_rel_VL)*(l_v-n_LV))-((B_v*cosd(delta_rel_VL))/2)-(n_S_VL*sind(delta_rel_VL)*cosd(delta_rel_VL)));
%     B32 = (psi_dot/M_Gier_ist)*((sind(delta_rel_VR)*(l_v*n_LV))+((B_v*cosd(delta_rel_VR))/2)-(n_S_VR*sind(delta_rel_VR)*cosd(delta_rel_VR)));
%     B33 = -((B_h*psi_dot)/(2*M_Gier_ist));
%     B34 = (B_h*psi_dot)/(2*M_Gier_ist);
% 
%     B = [ B11 B12 B13 B14;
%         B21 B22 B23 B24;
%         B31 B32 B33 B34];
%     
%     
%     C = [1 0 0; 0 0 1];
% 
%     E11 = ((sind(beta)*cosd(delta_rel_VL))-(cosd(beta)*sind(delta_rel_VL)))/(m_Fzg);
%     E12 = ((sind(beta)*cosd(delta_rel_VR))-(cosd(beta)*sind(delta_rel_VR)))/(m_Fzg);
%     E13 = (sind(beta))/(m_Fzg*v_Fzg);
%     E14 = (sind(beta))/(m_Fzg*v_Fzg);
%     E21 = ((cosd(beta)*cosd(delta_rel_VL))+(sind(beta)*sind(delta_rel_VR)))/(m_Fzg*v_Fzg);
%     E22 = ((cosd(beta)*cosd(delta_rel_VR))+(sind(beta)*sind(delta_rel_VR)))/(m_Fzg*v_Fzg);
%     E23 = (cosd(beta))/(m_Fzg*v_Fzg);
%     E24 = (cosd(beta))/(m_Fzg*v_Fzg);
%     E31 = (psi_dot/M_Gier_ist)*((cosd(delta_rel_VL)*(l_v*n_LV))+((B_v*sind(delta_rel_VL))/2)+(n_S_VL*(sind(delta_rel_VL))^2));
%     E32 = (psi_dot/M_Gier_ist)*((cosd(delta_rel_VR)*(l_v*n_LV))+((B_v*sind(delta_rel_VR))/2)+(n_S_VR*(sind(delta_rel_VR))^2));
%     E33 = -(l_h+n_LH);
%     E34 = -(l_h+n_LH);
% 
%     E = [ E11 E12 E13 E14;
%         E21 E22 E23 E24;
%         E31 E32 E33 E34];
    
end

%     x_vec_ist = x_vec_ist_k+T_ab *x_dot_ist_k;
%     x_dot_ist = A*x_vec_ist + B*u+E*d;



