function A_spalte = berechneABlockSpalte(A, N2)

    A_spalte = zeros(3 * N2, 3);
 
    for m = 1:N2
        zeilen = (3 * (m - 1) + 1):(3 * m);
        A_spalte(zeilen, :) = A ^ m;
    end

end