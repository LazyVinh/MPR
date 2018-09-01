function sum_A_j = berechnePotenzSumme(A, start, ende)

    sum_A_j = zeros(3, 3);
    
    for j = start:ende
        sum_A_j = sum_A_j + A ^ j;
    end

end