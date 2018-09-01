function CA_spalte = berechneCABlockSpalte(A, C, N2)

    CA_spalte = zeros(2 * N2, 3);
 
    for m = 1:N2
        zeilen = (2 * (m - 1) + 1):(2 * m);
        CA_spalte(zeilen, :) = C * A ^ m;
    end

end