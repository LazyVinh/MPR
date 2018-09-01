%
%   Siehe Formel (10.14) im Skript, Seite 341 H
%
function f = berechneFreieRegelgroesse(A, B, C, N2, x, u)

    CA_spalte = berechneCABlockSpalte(A, C, N2);
    phi_spalte = berechnePhiBlockSpalte(A, B, C, N2, 1);
 
    f = CA_spalte * x + phi_spalte * u;

end