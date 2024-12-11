close all;
clear all;
clc;

% Paraméterek
L = 1e-3; % Tekercs induktivitása (H)
C = 1e-6; % Kondenzátor kapacitása (F)
R0 = 10;  % Bemeneti ellenállás (Ohm)
Rt = 2 * sqrt(L / C); % Kimeneti ellenállás (Ohm)
U0 = 5;  % Bemeneti feszültség (V)

fprintf('Az áramkör paraméterei: L = %f H, C = %f F, R0 = %f Ohm, Rt = %f Ohm, U0 = %f V\n', L, C, R0, Rt, U0);

A = [-R0 / L, -1 / L; 1 / C, -1/(Rt*C)];

lambda = eig(A);

tau = -1 ./ real(lambda);

fprintf('A rendszer sajátértékei: %s\n', mat2str(lambda));

fprintf('A rendszer időállandója: tau = %.4f s\n', tau);

% Időtartomány
tspan = [0 5*max(tau)]; % 0-tól 5*tau-ig

% Állandó bemenő feszültség
Usf = U0;

% Kezdeti feltételek: iL(0) és uC(0)
initial_conditions = [0; 0];

% y(1) = iL, y(2) = uC
% Differenciálegyenletek definíciója
odefun = @(t, y, us) [
    -R0/L * y(1) - y(2)/L + us/L  ; % diL/dt = -R0/L * iL - uC/L + U0/L
    y(1)/C - y(2)/(Rt*C)           ; % duC/dt = iL/C - uC/(Rt*C)
];

% Numerikus megoldás
[T, Y] = ode45(@(t,x) odefun(t,x,Usf), tspan, initial_conditions);

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
