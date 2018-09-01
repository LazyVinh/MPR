function I_Dreieck = berechneI_Dreieck(Nu)

I_Dreieck = zeros(4*Nu,4*Nu);
I = eye(4);
O = zeros(4,4);
    for n = 1:Nu
        for m = 1:Nu 
            if n<m
                I_Dreieck((4*(n-1)+1):((4*(n-1)+1)+3),(4*(m-1)+1):((4*(m-1)+1)+3)) = O;
            else
                I_Dreieck((4*(n-1)+1):((4*(n-1)+1)+3),(4*(m-1)+1):((4*(m-1)+1)+3)) = I;
            end
        end
    end

end