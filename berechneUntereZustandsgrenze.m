function x_lb_vec = berechneUntereZustandsgrenze(v_des, beta_lim, N2)

x_lb_vec = zeros(3*N2,1);
x_lb = [v_des; -beta_lim; -10];

    for i = 1:N2
         zeilen = (3*(i-1)+1);
        x_lb_vec(zeilen:(zeilen+2),1) = x_lb;
    end

end