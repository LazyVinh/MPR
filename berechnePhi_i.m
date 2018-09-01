function phi_i = berechnePhi_i(A, B, C, i)

    if i < 1
        phi_i = zeros(2, 4);
    else
        phi_i = C * berechnePotenzSumme(A, 0, i - 1) * B;
    end

end