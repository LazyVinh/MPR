N2 = 5;
Nu = 3;
T_ab = 0.01;
y_soll = 3;
Y_E = 3;
I_k = eye(4);
C = [1 0 0; 0 0 1];
A = [2 3 4; 3 4 5; 4 5 6];
B = [ 12 3 4 5; 2 3 4 5 ; 2 3 4 5];
g_u = [ 1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
Q = [1 0; 0 1];
O = ones(2,3);
x_vec_soll = [16.67; 0.5; 4];
x_p = [22.33; 3; 4];
v_0 = x_p(1);
v_Fzg = dot(x_p, [1; 0; 0]);
u_k = [ -0.3; -0.3; -0.3; -0.3];
Qk = [2 0 0; 0 3 0; 0 0 4];
S = [2 0 0; 0 3 0; 0 0 4];
x_dot_p = [-2; 1; 1];
t = 1;
f1 = [1; 0; 0];
f2 = [0; 1; 0];
f3 = [0; 0; 1];
v_des = 16.67;
psi = zeros(1,1);
R = eye(4);
y_soll = [ 5 ; 1];
y_p = [0; 0];
A_kn = zeros(3,3);
dM_Gier = -100;
u = ones(4,1);
H = ones(3, 3);
t_sim = tic;
beta_lim = 10;

% Berechnung der lower und upper bound der Bremskräfte
m_Fzg = 6800;
a_x0 = 0;
a_x_TB2 = -6;
a_x_max = -2.3;
F_B_max = (a_x_TB2+a_x_max)*(m_Fzg/4);
u_max = [F_B_max;F_B_max; F_B_max; F_B_max];
F_B0 = a_x0*m_Fzg;
u_0 = [F_B0; F_B0; F_B0; F_B0];
dF_B_max = -a_x_max*m_Fzg/4;
du_max = [dF_B_max; dF_B_max; dF_B_max; dF_B_max];

lambda = 0.6;

r_k = berechneReferenztrajektorie(T_ab, N2, y_soll, y_p);

addpath MPR

f_k = berechneFreieRegelgroesse(A, B, C, N2, x_p, u_k);
phi = berechnePhiBlockMatrix(A, B, C, N2, Nu);
J = @(du) berechneKostenfunktion(f_k, phi, du, r_k, lambda);
f_x = berechneFreieRegelgroesseOhneC(A, B, N2, x_p, u_k);
B_matrix = berechneBBlockMatrix(A, B, N2, Nu);
x_ub_vec = berechneUntereZustandsgrenze(v_0, beta_lim, N2);
x_lb_vec = berechneObereZustandsgrenze(v_des, beta_lim, N2);
u_k_spalte = berechneU_kSpalte(u_k, Nu);
u_max_spalte = berechneU_max_spalte(u_max, Nu);
I_Dreieck = berechneI_Dreieck(Nu);

A_ineq1 = B_matrix;
A_ineq2 = -B_matrix;
A_ineq3 = I_Dreieck;
A_ineq4 = I_Dreieck;

A_ineq = [A_ineq1; A_ineq2; A_ineq3; A_ineq4];

b_ineq1 = x_ub_vec-f_x;
b_ineq2 = x_lb_vec+f_x;
b_ineq3 = -u_k_spalte;
b_ineq4 = u_k_spalte-u_max_spalte;

b_ineq = [b_ineq1; b_ineq2; b_ineq3; b_ineq4];

du_0 = zeros(4*Nu,1);
 
du = fmincon(J, du_0, A_ineq, b_ineq);

du_k = du(1:4);

u = u_k+du_k;
 
 
         