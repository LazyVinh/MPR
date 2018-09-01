function AB_spalte = berechneABBlockSpalte(A, B, N2)
  
  AB_spalte = zeros(3*N2, 4);
  
  for m = 1:N2
      zeilen = (3*(m - 1) + 1):(3*m);
      AB_spalte(zeilen, 1:4) = berechneAB(A, B, m);
  end
  
end