%
%    Berechnet die grosse Block-Matrix in Formel (10.13),
%    siehe Skript, Seite 341 H
%
function phi_matrix = berechnePhiBlockMatrix(A, B, C, N2, Nu)
  
phi_matrix = zeros(2*N2, 4*Nu);

  for blockSpaltenIndex = 1:Nu
    spalten = (4*(blockSpaltenIndex - 1) + 1):(4*blockSpaltenIndex);
    phi_matrix(:, spalten) = berechnePhiBlockSpalte(A, B, C, N2, blockSpaltenIndex);
  end

end