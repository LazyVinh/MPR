function u_max_spalte = berechneU_max_spalte(u_max, Nu)

u_max_spalte = zeros(4*Nu,1);

for i = 1:Nu
    zeilen = 4*(i-1)+1;
    u_max_spalte(zeilen:(zeilen + 3), 1) = u_max;
end