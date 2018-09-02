function phi_i = berechnePhi_i(A, B, C, i)

  m = size(C, 1);
  n = size(B, 2);

  if i < 1
      phi_i = zeros(m, n);
  else
      phi_i = C * berechnePotenzSumme(A, 0, i - 1) * B;
  end

end