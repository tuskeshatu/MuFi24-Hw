close all;
clear all;
clc;

% Paraméterek
L = 1e-5; % Tekercs induktivitása (H)
C = 1e-8; % Kondenzátor kapacitása (F)
R0 = 10;  % Bemeneti ellenállás (Ohm)
Rt = 2 * sqrt(L / C); % Kimeneti ellenállás (Ohm)
U0 = 5;  % Bemeneti feszültség (V)
n = 500;   % Szakaszok száma

% Időtartomány
tau = 2 * sqrt(L * C); % Időállandó
tspan = [0 5*tau]; % 0-tól 5*tau-ig

% Kezdeti feltételek: iL(0) és uC(0) minden szakaszra
initial_conditions = zeros(2*n, 1);

% Differenciálegyenletek definíciója
odefun = @(t, y) [
    -R0/L * y(1) - y(2)/L + U0/L  ; % 1. szakasz diL/dt
    y(1)/C - y(2)/(Rt*C)           ; % 1. szakasz duC/dt
    % További szakaszok
    arrayfun(@(k) [
        -R0/L * y(2*k-1) - y(2*k)/L + U0/L  ; % k. szakasz diL/dt
        y(2*k-1)/C - y(2*k)/(Rt*C)           % k. szakasz duC/dt
    ], 2:n, 'UniformOutput', false)
];

% Numerikus megoldás
[t, y] = ode45(@(t, y) odefun_combined(t, y, n, R0, L, C, Rt, U0), tspan, initial_conditions);

% Eredmények
iL = y(:, 1:2:end); % Tekercs áramai minden szakaszra
uC = y(:, 2:2:end); % Kondenzátor feszültségei minden szakaszra

% Ábrázolás
figure;
subplot(3, 1, 1);
plot(t, iL(:, 1), 'b', t, uC(:, 1), 'r');
title('1st Section');
xlabel('Time (s)');
ylabel('Current (A) / Voltage (V)');
legend('iL', 'uC');

subplot(3, 1, 2);
plot(t, iL(:, ceil(n/2)), 'b', t, uC(:, ceil(n/2)), 'r');
title('Middle Section');
xlabel('Time (s)');
ylabel('Current (A) / Voltage (V)');
legend('iL', 'uC');

subplot(3, 1, 3);
plot(t, iL(:, end), 'b', t, uC(:, end), 'r');
title('Last Section');
xlabel('Time (s)');
ylabel('Current (A) / Voltage (V)');
legend('iL', 'uC');

function dydt = odefun_combined(t, y, n, R0, L, C, Rt, U0)
    dydt = zeros(2*n, 1);
    for k = 1:n
        dydt(2*k-1) = -R0/L * y(2*k-1) - y(2*k)/L + U0/L; % diL/dt
        dydt(2*k) = y(2*k-1)/C - y(2*k)/(Rt*C);           % duC/dt
    end
end