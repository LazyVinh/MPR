N2 = 4;
Nu = 3;
T_ab = 0.01;
y_soll = 3;
Y_E = 3;
I_k = eye(4);
C = [1 0 0; 0 0 1];
A = [0 0 0; 0 -1.489 -0.9887; 0 5.921 -2.439];
B = [ 0.0001471 0.0001471 0.0001471 0.0001471; 0 0 0 0 ; -0.0001483 0.0001483 -0.000114 0.000114];
x_vec_soll = [16.67; 0.5; 4];
x_p = [22.33; 0; 0];
v_0 = x_p(1);
v_Fzg = dot(x_p, [1; 0; 0]);
u_k = [ -0.3; -0.3; -0.3; -0.3];
x_dot_p = [-2; 1; 1];
t = 1;
f1 = [1; 0; 0];
f2 = [0; 1; 0];
f3 = [0; 0; 1];
v_des = 16.67;
psi = zeros(1,1);
y_soll = [ 13 ; 4];
y_p = [0; 0];
dM_Gier = -100;
u = ones(4,1);
t_sim = tic;
beta_lim = 10;
S = eig(A);

% Berechnung der lower und upper bound der Bremskräfte
m_Fzg = 6800;
a_x0 = 0;
a_x_TB2 = -6;
a_x_max = -2.3;
F_B_max = (a_x_TB2+a_x_max)*(m_Fzg/4);
u_max = [F_B_max;F_B_max; F_B_max; F_B_max];
F_B0 = a_x0*m_Fzg;
u_0 = [F_B0; F_B0; F_B0; F_B0];
dF_B_max = a_x_max*m_Fzg/4;
du_max = [dF_B_max; dF_B_max; dF_B_max; dF_B_max];

lambda = 0.3;
v_max = 100/3.6;
v_min = 0;

r_k = berechneReferenztrajektorie(T_ab, N2, y_soll, y_p);

addpath MPR

f_k = berechneFreieRegelgroesse(A, B, C, N2, x_p, u_k);
phi = berechnePhiBlockMatrix(A, B, C, N2, Nu);
J = @(du) berechneKostenfunktion(f_k, phi, du, r_k, lambda);
f_x = berechneFreieRegelgroesseOhneC(A, B, N2, x_p, u_k);
B_matrix = berechnePhiBlockMatrix(A, B, eye(3), N2, Nu);
x_ub_vec = berechneObereZustandsgrenze(v_max, beta_lim, N2);
x_lb_vec = berechneUntereZustandsgrenze(v_min, beta_lim, N2);
u_k_spalte = berechneU_kSpalte(u_k, Nu);
u_max_spalte = berechneU_max_spalte(u_max, Nu);
I_Dreieck = berechneI_Dreieck(Nu);

A_ineq1 = B_matrix;
A_ineq2 = -B_matrix;
A_ineq3 = I_Dreieck;
A_ineq4 = -I_Dreieck;

A_ineq = [A_ineq1; A_ineq2; A_ineq3; A_ineq4];

b_ineq1 = x_ub_vec-f_x;
b_ineq2 = -x_lb_vec+f_x;
b_ineq3 = zeros(4*Nu,1)-u_k_spalte;
b_ineq4 = u_k_spalte-u_max_spalte;

b_ineq = [b_ineq1; b_ineq2; b_ineq3; b_ineq4];

du_0 = zeros(4*Nu,1);

f = zeros(size(du_0));

du_new = linprog(f, A_ineq, b_ineq);
 
du = fmincon(J, du_new, A_ineq, b_ineq);
 
du_k = du(1:4);
 
u = u_k+du_k;    