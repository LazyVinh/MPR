%
%  Berechnet die grosse Block-Matrix in Formel (10.13),
%  siehe Skript, Seite 341 H
%
function phi_matrix = berechnePhiBlockMatrix(A, B, C, N2, Nu)
 
  m = size(C, 1);
  n = size(B, 2);

  phi_matrix = zeros(m*N2, n*Nu);

  for blockSpaltenIndex = 1:Nu
    spalten = (n*(blockSpaltenIndex - 1) + 1):(n*blockSpaltenIndex);
    phi_matrix(:, spalten) = berechnePhiBlockSpalte(A, B, C, N2, blockSpaltenIndex);
  end

end