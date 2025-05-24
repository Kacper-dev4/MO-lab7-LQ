clear all; close all; clc;

% Parametry
N = 20;                     % liczba iteracji
Q = [25 -5; -5 1];          % z (-5*x1 + x2)^2 = x'*Q*x
R = 19;
A = [1 1; 0 2];
B = [1; 2];

% --- 1. Sprawdzenie sterowalności ---
Co = ctrb(A, B);
if rank(Co) == size(A,1)
    disp('Układ jest sterowalny.');
else
    disp('Układ NIE jest sterowalny.');
end

% --- Rekurencja Riccatiego dla zadanej wartości R ---
P = zeros(2,2,N+1);
P(:,:,N+1) = zeros(2); % terminal cost = 0

for k = N:-1:1
    Pk = P(:,:,k+1);
    K = (R + B'*Pk*B) \ (B'*Pk*A);
    P(:,:,k) = Q + A'*Pk*A - A'*Pk*B*K;
end

% --- 3. Wpływ różnych x0 ---
x0_vals = {[10; 15], [5; 5], [-5; 10], [20; -10]};
colors = ['r', 'g', 'b', 'k'];
figure;

for idx = 1:length(x0_vals)
    x = zeros(2, N+1);
    u = zeros(1, N);
    x(:,1) = x0_vals{idx};
    for k = 1:N
        K = (R + B' * P(:,:,k+1) * B) \ (B' * P(:,:,k+1) * A);
        u(k) = -K * x(:,k);
        x(:,k+1) = A * x(:,k) + B * u(k);
    end
    subplot(3,1,1); hold on;
    plot(0:N, x(1,:), 'Color', colors(idx)); title('x_1 dla różnych x_0');
    subplot(3,1,2); hold on;
    plot(0:N, x(2,:), 'Color', colors(idx)); title('x_2 dla różnych x_0');
    subplot(3,1,3); hold on;
    plot(0:N-1, u, 'Color', colors(idx)); title('u dla różnych x_0');
end

% --- 4. Wpływ wartości R ---
R_vals = [1, 10, 19, 100];
figure;

for idx = 1:length(R_vals)
    R = R_vals(idx);
    % Rekurencja Riccatiego
    P = zeros(2,2,N+1);
    P(:,:,N+1) = zeros(2);
    for k = N:-1:1
        Pk = P(:,:,k+1);
        K = (R + B'*Pk*B) \ (B'*Pk*A);
        P(:,:,k) = Q + A'*Pk*A - A'*Pk*B*K;
    end
    x = zeros(2, N+1);
    u = zeros(1, N);
    x(:,1) = [10; 15];
    for k = 1:N
        K = (R + B' * P(:,:,k+1) * B) \ (B' * P(:,:,k+1) * A);
        u(k) = -K * x(:,k);
        x(:,k+1) = A * x(:,k) + B * u(k);
    end
    subplot(3,1,1); hold on;
    plot(0:N, x(1,:)); title('x_1 dla różnych R');
    subplot(3,1,2); hold on;
    plot(0:N, x(2,:)); title('x_2 dla różnych R');
    subplot(3,1,3); hold on;
    plot(0:N-1, u); title('u dla różnych R');
end
legend('R=1','R=10','R=19','R=100');

% --- 5. Elementy macierzy K w czasie ---
R = 19;
P = zeros(2,2,N+1);
P(:,:,N+1) = zeros(2);
for k = N:-1:1
    Pk = P(:,:,k+1);
    K = (R + B'*Pk*B) \ (B'*Pk*A);
    P(:,:,k) = Q + A'*Pk*A - A'*Pk*B*K;
end

K_hist = zeros(N, 2);  % zapisujemy obie składowe K
for k = 1:N
    K = (R + B' * P(:,:,k+1) * B) \ (B' * P(:,:,k+1) * A);
    K_hist(k,:) = K;
end

figure;
plot(1:N, K_hist(:,1), '-o'); hold on;
plot(1:N, K_hist(:,2), '-x');
legend('K(1)','K(2)');
title('Elementy macierzy K w czasie');
xlabel('i'); ylabel('K_{i}');
