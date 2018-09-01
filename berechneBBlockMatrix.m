function B_matrix = berechneBBlockMatrix(A, B, N2, Nu)
  
B_matrix = zeros(3*N2, 4*Nu);

  for blockSpaltenIndex = 1:Nu
    spalten = (4*(blockSpaltenIndex - 1) + 1):(4*blockSpaltenIndex);
    B_matrix(:, spalten) = berechneABBlockSpalte(A, B, N2);
  end

end