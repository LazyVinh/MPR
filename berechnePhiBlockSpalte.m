function phi_spalte = berechnePhiBlockSpalte(A, B, C, N2, spaltenIndex)
  
  phi_spalte = zeros(2*N2, 4);
  
  for m = 1:N2
    zeilen = (2*(m - 1) + 1):(2*m);
      i = m - (spaltenIndex - 1);
      phi_spalte(zeilen, 1:4) = berechnePhi_i(A, B, C, i);
  end
  
end