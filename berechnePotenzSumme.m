function sum_A_j = berechnePotenzSumme(A, start, ende)

  [m,n] = size(A);

  sum_A_j = zeros(m, n);

  for j = start:ende
      sum_A_j = sum_A_j + A ^ j;
  end

end