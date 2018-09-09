Y_E = 3; % Endquerversatz
y_ist = 0.2; % aktuelle Querposition
a_y_g = 5; % maximale Querbeschleunigung bis zur Kippgrenze
t_sim = tic; % Simulationslaufzeit;
dx_min_A_n = 60; % minimal benötigter Abstand für ein Ausweichmanöver
y_0 = 0; % Anfangsquerposition
v_y0 = 0; % Anfangsquergeschwindigkeit
a_y0 = 0; % Anfangsquerbeschleunigung
dx1 = 3; %  Abstand von Hindernis und Egofahrzeug bis zum Crash
%  zählt ab erreichen des gebrauchten Mindestabstand für
%  Ausweichen von 0 bis minimal benötigten Abstand hoch
Y_soll_k = 0.4; % Vergangene Sollquerposition der Trajektorie


X_E = dx_min_A_n; % Endzustand X-Position
x = X_E;
v_yE = 0;       % Endquergeschwindigkeit der Trajektorie
a_yE = 0;       % Endquerbeschleunigung der Trajektorie


if t_sim == 0       % Simulationsstart
    Y_soll = 0;     % Definitere Anfangs-Querposition, -geschwindigkeit und -beschleunigung
    Y_soll_dot = 0;
    Y_soll_2dot = 0;
    
else
    y = Y_E ;
    v_y = v_yE ;
    a_y = a_yE ;
    
    if abs(a_y)>=abs(a_y_g/2) % Begrenzung der maximalen Querbeschleunigung auf hälfte der Kippgrenzbeschleunigung
        a_y = a_y_g/2;
    end
    
    % Koeffizienten für die S-Funktion Ausweichbahn
    a0 = y_0;
    a1 = (v_y0/3.6); % m/s
    a2 = (a_y0/2); % m/s^2
    a3 = ((x^2)*(-a_y+(3*a_y0))-x*((8*v_y/3.6)+(12*v_y0/3.6))+(20*y))/(2*x^3);
    a4 = ((x^2)*((-2*a_y)+(3*a_y0))+x*((14*v_y/3.6)+(16*v_y0/3.6))-(30*y))/(2*x^4);
    a5 = ((x^2)*(a_y-a_y0)-(6*x)*((v_y+v_y0)/3.6)+(12*y))/(2*x^5);
    
if dx1 <=((dx_min_A_n)/2) && dx1>=0   % Abschnitt der S-kurvenanfang und vor dem Wendepunkt
    x = ((dx_min_A_n)/2)-dx1;
    if Y_E>0 || Y_E<0       % wenn links oder rechts ausgewichen wird soll s-kurve ausgegeben werden
        % Sollgrößen in y-Richtung
        Y_soll = (a5*x^5)+(a4*x^4)+(a3*x^3)+(a2*x^2)+(a1*x)+a0;
        Y_soll_dot = ((5*a5*x^4)+(4*a4*x^3)+(3*a3*x^2)+(2*a2*x)+a1)*3.6;
        Y_soll_2dot = (20*a5.*x^3)+(12*a4.*x^2)+(6*a3.*x)+(a2.*2);
    else
        Y_soll = 0;         % andernfalls keine s-kurve als sollgröße
        Y_soll_dot = 0;
        Y_soll_2dot = 0;
    end
elseif dx1>(dx_min_A_n/2) % Abschnitt vor beginn der S-Kurve 
    x = ((dx_min_A_n)/2)-dx1;
    Y_soll = 0;
    Y_soll_dot = 0;
    Y_soll_2dot = 0;
else            % Abschnitt der S-kurve und nach dem Wendepunkt
        x = ((dx_min_A_n)/2)-dx1;
        if (Y_E >0) || (Y_E<0)      % Sollbahn unabhängig ob links oder rechts ausgewichen wird
            if (Y_soll_k < Y_E)||(Y_soll_k > Y_E) % wenn gerade angegebene sollgröße den endquerversatz noch nicht erreicht hat
                % Sollgrößen in y-Richtung
                Y_soll = (a5*x^5)+(a4*x^4)+(a3*x^3)+(a2*x^2)+(a1*x)+a0;
                Y_soll_dot = ((5*a5*x^4)+(4*a4*x^3)+(3*a3*x^2)+(2*a2*x)+a1)*3.6;
                Y_soll_2dot = (20*a5.*x^3)+(12*a4.*x^2)+(6*a3.*x)+(a2.*2);
            else
                Y_soll = Y_E;       % andernfalls endzustand als sollgröße
                Y_soll_dot = 0;
                Y_soll_2dot = 0;
            end
        else
            Y_soll = 0;
            Y_soll_dot = 0;
            Y_soll_2dot = 0;
        end
end
X_soll = x;
end






