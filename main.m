clear all

J = @(x,u)sum((-5*x(1) + x(2))^2 + 19*u^2);

A = [1,1;
     0,2];
B = [1;
     2];

C = [B A*B];
sterowalnosc = 2==rank(C);
R = 19;
F = [0,0;
     0,0];
Q = [25,-5;  
     -5,1];
N = 20;
x0 = [10;
      15];
x = zeros(2,N+1);
u = zeros(1,N);
x(:,1) = x0;


K(:,:,N+1) = F;


for i=N:-1:1 
    K(:,:,i) = A'*(K(:,:,i+1) - K(:,:,i+1)*B*((R+B'*K(:,:,i+1)*B)^(-1))*B'*K(:,:,i+1))*A + Q;
end

for i=1:N
    S = -(R+B'*K(:,:,i+1)*B)^(-1)*B'*K(:,:,i+1)*A;

    u(i) = S*x(:,i);
    x(:,i+1) = A*x(:,i) + B*u(i);
end

J0 = (1/2)*x0'*K(:,:,1)*x0;

figure
plot(x(1,:),'x','LineStyle','--')
xlabel('Iteracja')
ylabel('Wartość x1 w danej iteracji')
title('x1')

figure
plot(x(2,:),'x','LineStyle','--')
xlabel('Iteracja')
ylabel('Wartość x2 w danej iteracji')
title('x2')

figure
plot(u,'o','LineStyle','--')
xlabel('Iteracja')
ylabel('Wartość u w danej iteracji')
title('u')

figure
hold on
for i=1:2
    for j=1:2
        Kpom = K(i,j,:);
        Kpom = Kpom(:);
        plot(Kpom)
    end
end
xlabel('Iteracja')
ylabel('Wartości elementów K w danej iteracji')
title('Zmiana elementów macierzy K')
legend('w:1 k:1', 'w:1 k:2', 'w:2 k:1', 'w:2 k:2');
hold off