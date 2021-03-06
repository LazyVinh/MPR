%
%   siehe Skript, Seite 341 H, Formel (10.15)
%
function J = berechneKostenfunktion(f, Phi, du, r, lambda, mue, gamma)
 
    % ueberpruefe Dimensionen
    [anzahlZeilen, anzahlSpalten] = size(Phi);
 
    if any(size(f) ~= [anzahlZeilen, 1])
        error("f hat falsche Dimensionen!")
    end
 
    if any(size(r) ~= [anzahlZeilen, 1])
        error("r hat falsche Dimensionen!")
    end
 
    if any(size(du) ~= [anzahlSpalten, 1])
        error("du hat falsche Dimensionen!")
    end
 
    % Q kann man in den Berechnungen eigntl weglassen,
    % weil Einheitsmatrix
    % Q = eye(anzahlZeilen, anzahlZeilen);
    Q = zeros(anzahlZeilen, anzahlZeilen);
    q = zeros(1,1);
    
    for i = 1:anzahlZeilen
        if rem(i,2)
            q(i) = gamma;
        else
            q(i) = mue;
        end
    end
    Q = diag(q);
 
    % R wird eigntl auch nicht benoetigt, einfach lambda
    % multiplizieren ist auch ok (sogar effizienter!)
    R = lambda * eye(anzahlSpalten, anzahlSpalten);
 
    % berechne freie Regelabweichung
    e = f - r;
 
    summand_1 = e' * Q * e;
    summand_2 = 2 * (e' * Q) * (Phi * du);
    summand_3 = du' * (Phi' * Q * Phi + R) * du;
 
    Q = isnan(summand_1 + summand_2 + summand_3);
    
    if Q ==true
        J = 0;
    else
        J = summand_1 + summand_2 + summand_3;
    end
 
end
