function phi_spalte = berechnePhiBlockSpalte(A, B, C, N2, spaltenIndex)
  
  [m, w] = size(C);
  [v, n] = size(B);
  
  phi_spalte = zeros(m*N2, n);
  
  for j = 1:N2
      zeilen = (m*(j - 1) + 1):(m*j);
      i = j - (spaltenIndex - 1);
      phi_i = berechnePhi_i(A, B, C, i);
      phi_spalte(zeilen, :) = phi_i;
  end
  
end