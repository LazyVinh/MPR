function r_k = berechneReferenztrajektorie(T_ab, N2, y_soll, y_p)

AA = zeros(1,1);
I = zeros(1,1);
for n = 1:2
    AA(n,n) = exp(-T_ab/n);
    I(n,n) = 1;
end

r_k = zeros(2*N2,1);
for n = 1:N2
    r_k((2*(n-1)+1):2*n,1) = (I-(AA)^n)*y_soll+(AA^n)*y_p;
end
