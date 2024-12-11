%% Paraméterek

close all;
clear all;
clc;

L = input('Adja meg az induktivitást (L) [H]: ');
while L <= 0
    disp('Az induktivitásnak pozitív valós számnak kell lennie.');
    L = input('Adja meg az induktivitást (L) [H]: ');
end

C = input('Adja meg a kapacitást (C) [F]: ');
while C <= 0
    disp('A kapacitásnak pozitív valós számnak kell lennie.');
    C = input('Adja meg a kapacitást (C) [F]: ');
end

R0 = input('Adja meg a forrás belső ellenállását (R0) [Ohm]: ');
while R0 <= 0
    disp('A forrás belső ellenállásának pozitív valós számnak kell lennie.');
    R0 = input('Adja meg a forrás belső ellenállását (R0) [Ohm]: ');
end
Rt = 2 * sqrt(L / C); % Lezáró ellenállás [Ohm]
U0 = input('Adja meg a feszültségforrás amplitúdóját (U0) [V]: ');
n = input('Adja meg a létrafokok számát (n): ');
while n <= 0
    disp('A létrafokok számának pozitív egész számnak kell lennie.');
    n = input('Adja meg a létrafokok számát (n): ');
end

tmax = input('Adja meg az időtartomány végét (tmax) [s] (0 ha automatikusan legyen beállítva): ');
while tmax < 0
    disp('Az időtartomány végének nem lehet negatív.');
    tmax = input('Adja meg az időtartomány végét (tmax) [s] (0 ha automatikusan legyen beállítva): ');
end
if tmax == 0
    tmax = 3 * n * sqrt(L * C);
end

% Szimuláció paramétereinek kiírása
disp('---------------------------------------------------------------------------------');
disp('Szimuláció paraméterei:');
fprintf('Induktivitás (L): %f H\n', L);
fprintf('Kapacitás (C): %f F\n', C);
fprintf('Forrás belső ellenállása (R0): %f Ohm\n', R0);
fprintf('Lezáró ellenállás (Rt): %f Ohm\n', Rt);
fprintf('Feszültségforrás amplitúdója (U0): %2f V\n', U0);
fprintf('Létrafokok száma (n): %d\n', n);
fprintf('Időtartomány vége (tmax): %f s\n', tmax);

% Időtartomány
tspan = [0 tmax]; % 0 - 3*létrafokszám*tau

% Kezdőállapot: minden áram és feszültség nulla
y0 = zeros(2 * n, 1);

% Numerikus megoldás
[t, y] = ode45(@(t, y) odefun(t, y, L, C, R0, Rt, U0, n), tspan, y0);

%% Ábrázolás

% Kondenátorok feszültségei

% Első, középső és utolsó kondenzátor kiválasztása
U_1 = y(:, n+1); % Első kondenzátor feszültsége
mid_idx = ceil(n/2); % Középső kondenzátor indexe
U_mid = y(:, n+mid_idx); % Középső kondenzátor feszültsége
U_last = y(:, end); % Utolsó kondenzátor feszültsége

figure;
plot(t, U_1, 'b', 'DisplayName', 'Első kondenzátor (U_{C1}');
hold on;
plot(t, U_mid, 'g', 'DisplayName', sprintf('Középső kondenzátor (U_{C%d})', mid_idx));
plot(t, U_last, 'r', 'DisplayName', sprintf('Utolsó kondenzátor (V_C_{%d})', n));
xlabel('Idő [s]');
ylabel('Feszültség [V]');
legend;
grid on;
title('Kondenzátorok feszültségei a távvezetéken');

% Tekercsek áramai
I_1 = y(:, 1); % Első tekercs árama
I_mid = y(:, mid_idx); % Középső tekercs árama
I_last = y(:, n); % Utolsó tekercs árama

figure;
plot(t, I_1, 'b', 'DisplayName', 'Első tekercs (I_{L1})');
hold on;
plot(t, I_mid, 'g', 'DisplayName', sprintf('Középső tekercs (I_{L%d})', mid_idx));
plot(t, I_last, 'r', 'DisplayName', sprintf('Utolsó tekercs (I_{L%d})', n));
xlabel('Idő [s]');
ylabel('Áram [A]');
legend;
grid on;
title('Tekercsek áramai a távvezetéken');

%% Egyéb user requestek

while true
    choice = input('Kíván egyéb dolgokat plotolni? (y/n): ', 's');
    if choice == 'n'
        break;
    elseif choice == 'y'
        plot_choice = input('Mit szeretne plotolni? (L ha tekercs áramot, C ha kondenzátor feszültséget): ', 's');
        if plot_choice == 'L'
            idx = input('Adja meg a tekercs indexét (1-től n-ig): ');
            if idx >= 1 && idx <= n
                figure;
                plot(t, y(:, idx), 'DisplayName', sprintf('Tekercs %d árama', idx));
                xlabel('Idő [s]');
                ylabel('Áram [A]');
                legend;
                grid on;
                title(sprintf('Tekercs %d árama a távvezetéken', idx));
            else
                disp('Érvénytelen index.');
            end
        elseif plot_choice == 'C'
            idx = input('Adja meg a kondenzátor indexét (1-től n-ig): ');
            if idx >= 1 && idx <= n
                figure;
                plot(t, y(:, n+idx), 'DisplayName', sprintf('Kondenzátor %d feszültsége', idx));
                xlabel('Idő [s]');
                ylabel('Feszültség [V]');
                legend;
                grid on;
                title(sprintf('Kondenzátor %d feszültsége a távvezetéken', idx));
            else
                disp('Érvénytelen index.');
            end
        else
            disp('Érvénytelen választás.');
        end
    else
        disp('Érvénytelen választás.');
    end
end

%% odefun segédfüggvény
% Differenciálegyenletek definíciója
% Visszaadja az állapotváltozók deriváltjára rendezett egyenletrendszer jobboldalát
% (ÁVLNA)
function dydt = odefun(t, y, L, C, R0, Rt, U0, n)
    % Állapotváltozók
    iL = y(1:n); % Tekercsek áramai
    uC = y(n+1:end); % Kondenzátorok feszültségei
    
    % Forrás feszültsége (0V t<0, 5V t>=0)
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
            U_be = Us - R0 * iL(i); % Első létrafok bemeneti feszültsége
        else
            U_be = uC(i-1); % Előző kondenzátor feszültsége
        end
        
        if i == n
            I_ki = uC(i) / Rt; % Lezáró ellenállás a végén
        else
            I_ki = iL(i+1); % Következő tekercs árama
        end
        
        % Tekercs karakterisztikája
        diL_dt(i) = (U_be - uC(i)) / L;
        
        % Kondenzátor karakterisztikája
        duC_dt(i) = (iL(i) - I_ki) / C;
    end
    
    % Kombinált eredmény
    dydt = [diL_dt; duC_dt];
end
