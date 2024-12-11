%% Parameterek

close all;
clear all;
clc;

L = input('Adja meg az induktivitas (L) [H]: ');
while L <= 0
    disp('Az induktivitasnak pozitiv valos szamnak kell lennie.');
    L = input('Adja meg az induktivitas (L) [H]: ');
end

C = input('Adja meg a kapacitas (C) [F]: ');
while C <= 0
    disp('A kapacitasnak pozitiv valos szamnak kell lennie.');
    C = input('Adja meg a kapacitas (C) [F]: ');
end

R0 = input('Adja meg a forras belso ellenallasat (R0) [Ohm]: ');
while R0 <= 0
    disp('A forras belso ellenallasanak pozitiv valos szamnak kell lennie.');
    R0 = input('Adja meg a forras belso ellenallasat (R0) [Ohm]: ');
end
Rt = 2 * sqrt(L / C); % Lezaro ellenallas [Ohm]
U0 = input('Adja meg a feszultsegforras amplitudojat (U0) [V]: ');
n = input('Adja meg a letrafokok szamat (n): ');
while n <= 0
    disp('A letrafokok szamanak pozitiv egesz szamnak kell lennie.');
    n = input('Adja meg a letrafokok szamat (n): ');
end

tmax = input('Adja meg az idotartomany veget (tmax) [s] (0 ha automatikusan legyen beallitva): ');
while tmax < 0
    disp('Az idotartomany vegenek nem lehet negativ.');
    tmax = input('Adja meg az idotartomany veget (tmax) [s] (0 ha automatikusan legyen beallitva): ');
end
if tmax == 0
    tmax = 3 * n * sqrt(L * C);
end

% Szimulacio parametereinek kiirasa
disp('---------------------------------------------------------------------------------');
disp('Szimulacio parameterei:');
fprintf('Induktivitas (L): %f H\n', L);
fprintf('Kapacitas (C): %f F\n', C);
fprintf('Forras belso ellenallasa (R0): %f Ohm\n', R0);
fprintf('Lezaro ellenallas (Rt): %f Ohm\n', Rt);
fprintf('Feszultsegforras amplitudoja (U0): %2f V\n', U0);
fprintf('Letrafokok szama (n): %d\n', n);
fprintf('Idotartomany vege (tmax): %f s\n', tmax);

% Idotartomany
tspan = [0 tmax]; % 0 - 3*letrafokszam*tau

% Kezdoallapot: minden aram es feszultseg nulla
y0 = zeros(2 * n, 1);

%% Numerikus megoldas ode45 solverrel
[t, y] = ode45(@(t, y) odefun(t, y, L, C, R0, Rt, U0, n), tspan, y0);

%% Abrazolas

% kondenzatorok feszultsegei

% Elso, kozepso es utolso kondenzator kivalasztasa
U_1 = y(:, n+1); % Elso kondenzator feszultsege
mid_idx = ceil(n/2); % Kozepso kondenzator indexe
U_mid = y(:, n+mid_idx); % Kozepso kondenzator feszultsege
U_last = y(:, end); % Utolso kondenzator feszultsege

figure;
plot(t, U_1, 'b', 'DisplayName', 'Elso kondenzator (U_{C1}');
hold on;
plot(t, U_mid, 'g', 'DisplayName', sprintf('Kozepso kondenzator (U_{C%d})', mid_idx));
plot(t, U_last, 'r', 'DisplayName', sprintf('Utolso kondenzator (V_C_{%d})', n));
xlabel('Ido [s]');
ylabel('Feszultseg [V]');
legend;
grid on;
title('Kondenzatorok feszultsegei a tavvezeteken');

% Tekercsek aramai
I_1 = y(:, 1); % Elso tekercs arama
I_mid = y(:, mid_idx); % Kozepso tekercs arama
I_last = y(:, n); % Utolso tekercs arama

figure;
plot(t, I_1, 'b', 'DisplayName', 'Elso tekercs (I_{L1})');
hold on;
plot(t, I_mid, 'g', 'DisplayName', sprintf('Kozepso tekercs (I_{L%d})', mid_idx));
plot(t, I_last, 'r', 'DisplayName', sprintf('Utolso tekercs (I_{L%d})', n));
xlabel('Ido [s]');
ylabel('Aram [A]');
legend;
grid on;
title('Tekercsek aramai a tavvezeteken');

%% Egyeb user requestek

while true
    choice = input('Kivan egyeb dolgokat plotolni? (y/n): ', 's');
    if choice == 'n'
        break;
    elseif choice == 'y'
        plot_choice = input('Mit szeretne plotolni? (L ha tekercs aramot, C ha kondenzator feszultseget): ', 's');
        if plot_choice == 'L'
            idx = input('Adja meg a tekercs indexet (1-tol n-ig): ');
            if idx >= 1 && idx <= n
                figure;
                plot(t, y(:, idx), 'DisplayName', sprintf('Tekercs %d arama', idx));
                xlabel('Ido [s]');
                ylabel('Aram [A]');
                legend;
                grid on;
                title(sprintf('Tekercs %d arama a tavvezeteken', idx));
            else
                disp('Ervenytelen index.');
            end
        elseif plot_choice == 'C'
            idx = input('Adja meg a kondenzator indexet (1-tol n-ig): ');
            if idx >= 1 && idx <= n
                figure;
                plot(t, y(:, n+idx), 'DisplayName', sprintf('kondenzator %d feszultsege', idx));
                xlabel('Ido [s]');
                ylabel('Feszultseg [V]');
                legend;
                grid on;
                title(sprintf('kondenzator %d feszultsege a tavvezeteken', idx));
            else
                disp('Ervenytelen index.');
            end
        else
            disp('Ervenytelen valasztas.');
        end
    else
        disp('Ervenytelen valasztas.');
    end
end

%% odefun segedfuggveny
% Differencialegyenletek definicioja
% Visszaadja az allapotvaltozok derivaltjara rendezett egyenletrendszer jobboldalat
% (AVLNA)
function dydt = odefun(t, y, L, C, R0, Rt, U0, n)
    % Allapotvaltozok
    iL = y(1:n); % Tekercsek aramai
    uC = y(n+1:end); % kondenzatorok feszultsegei
    
    % Forras feszultsege (0V t<0, 5V t>=0)
    if t >= 0
        Us = U0;
    else
        Us = 0;
    end
    
    % Egyenletek
    diL_dt = zeros(n, 1);
    duC_dt = zeros(n, 1);
    
    for i = 1:n
        if i == 1
            U_be = Us - R0 * iL(i); % Elso letrafok bemeneti feszultsege
        else
            U_be = uC(i-1); % Elozo kondenzator feszultsege
        end
        
        if i == n
            I_ki = uC(i) / Rt; % Lezaro ellenallas a vegen
        else
            I_ki = iL(i+1); % Kovetkezo tekercs arama
        end
        
        % Tekercs karakterisztikaja
        diL_dt(i) = (U_be - uC(i)) / L;
        
        % Kondenzator karakterisztikaja
        duC_dt(i) = (iL(i) - I_ki) / C;
    end
    
    % Kombinalt eredmeny
    dydt = [diL_dt; duC_dt];
end
