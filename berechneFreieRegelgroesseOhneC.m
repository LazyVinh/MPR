function f_x = berechneFreieRegelgroesseOhneC(A, B, N2, x, u)

    A_spalte = berechneABlockSpalte(A, N2);
    AB_spalte = berechneABBlockSpalte(A, B, N2);
 
    f_x = A_spalte * x + AB_spalte * u;

end