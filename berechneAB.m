function AB = berechneAB(A, B, i)

    if i < 1
        AB = zeros(3,4);
    else
        AB = berechnePotenzSumme(A, 0, i - 1) * B;
    end

end