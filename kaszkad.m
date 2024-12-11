close all;
clear all;
clc;

% Paraméterek
L = 1e-3; % Tekercs induktivitása (H)
C = 1e-6; % Kondenzátor kapacitása (F)
R0 = 10;  % Bemeneti ellenállás (Ohm)
Rt = 2 * sqrt(L / C); % Kimeneti ellenállás (Ohm)
U0 = 5;  % Bemeneti feszültség (V)
n = 20;

% Időtartomány
tspan = [0 0.001]; % 0-tól 5*tau-ig

% Kezdeti feltételek: iL(0) és uC(0)
initial_conditions = zeros(2*n, 1);

% y(paratlan) = iL, y(páros) = uC

% Differenciálegyenletek definíciója
odefun = @(t, y) [
    -R0/L * y(1) - 1/L * y(2) + 1/L * U0; % 1. egyenlet: diL/dt = -R0/L * iL1 - uC1/L + U0/L

    % 2. - (n-1). egyenletek
    transpose(arrayfun(@(k) 1/L * y(2*k-2) - 1/L * y(2*k), 2:n-1))

    % További egyenletek
    transpose(arrayfun(@(k) 1/C * y(2*k-1) - 1/C * y(2*k+1), 2:n-1))

    1/C * y(2*n-1) - y(2*n)/(Rt*C); % 2*n. egyenlet: duC/dt = iLn/C - uCn/(Rt*C)
];

%%
% Numerikus megoldás
[T, Y] = ode45(@(t,y) odefun(t,y), tspan, initial_conditions);

% Eredmények
iL = Y(:, 1); % Tekercs áramai
uC = Y(:, 2); % Kondenzátor feszültségei

% Ábrázolás
figure;
subplot(2, 1, 1);
plot(T, iL, 'b-', 'LineWidth', 1.5);
xlabel('Idő (s)');
ylabel('Tekercs áram (A)');
title('Tekercs áram időfüggése');
grid on;

subplot(2, 1, 2);
plot(T, uC, 'r-', 'LineWidth', 1.5);
xlabel('Idő (s)');
ylabel('Kondenzátor feszültség (V)');
title('Kondenzátor feszültség időfüggése');
grid on;
