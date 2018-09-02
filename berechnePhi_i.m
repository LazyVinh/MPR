function phi_i = berechnePhi_i(A, B, C, i)

  [m, _] = size(C);
  [_, n] = size(B);

  if i < 1
      phi_i = zeros(m, n);
  else
      phi_i = C * berechnePotenzSumme(A, 0, i - 1) * B;
  end

end