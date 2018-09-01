function u_k_spalte = berechneU_kSpalte(u_k, Nu)

u_k_spalte = zeros(4*Nu,1);

for i = 1:Nu
    zeilen = (4*(i-1)+1);
    u_k_spalte(zeilen:(zeilen+3),1) = u_k;
end