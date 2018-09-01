function x_ub_vec = berechneUntereZustandsgrenze(v_0, beta_lim, N2)

x_ub_vec = zeros(3*N2,1);
x_ub = [v_0; beta_lim; 10];

    for i = 1:N2
         zeilen = (3*(i-1)+1);
        x_ub_vec(zeilen:(zeilen+2),1) = x_ub;
    end

end