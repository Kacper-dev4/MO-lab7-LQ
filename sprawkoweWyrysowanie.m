clear all
close all

A = [1,0;
     1,2];
B = [1;
     2];

C = [B A*B];
sterowalnosc = 2==rank(C);

F = [0,0;
     0,0];
Q = [16,-12;  
     -12,9];
N = 20;


R = 12;

x0 = [2,10,30;
      3,15,45];


for ii=1:3
x = zeros(2,N+1);
u = zeros(1,N+1);
x(:,1) = x0(:,ii);
K(:,:,N+1) = F;
for i=N:-1:1 
    K(:,:,i) = A'*(K(:,:,i+1) - K(:,:,i+1)*B*((R+B'*K(:,:,i+1)*B)^(-1))*B'*K(:,:,i+1))*A + Q;
end

for i=1:N
    S = -(R+B'*K(:,:,i+1)*B)^(-1)*B'*K(:,:,i+1)*A;

    u(i) = S*x(:,i);
    x(:,i+1) = A*x(:,i) + B*u(i);
    if i==20
     S = -(R+B'*K(:,:,i+1)*B)^(-1)*B'*K(:,:,i+1)*A;
     u(i+1) = S*x(:,i+1);
    end
end

J0 = (1/2)*x0'*K(:,:,1)*x0;

figure(1)
plot(0:20,x(1,:),'x','LineStyle','--')
xlabel('Iteracja')
ylabel('Wartość x1 w danej iteracji')
title('x1, R=12')
hold on

figure(2)
plot(0:20,x(2,:),'x','LineStyle','--')
xlabel('Iteracja')
ylabel('Wartość x2 w danej iteracji')
title('x2, R=12')
hold on

figure(3)
plot(0:20,u,'o','LineStyle','--')
xlabel('Iteracja')
ylabel('Wartość u w danej iteracji')
title('u, R=12')
hold on

% figure
% hold on
% for i=1:2
%     for j=1:2
%         Kpom = K(i,j,:);
%         Kpom = Kpom(:);
%         if ~(i==1 && j==2)
%             stairs(0:20,Kpom,'o--')
%         else
%             stairs(0:20,Kpom,'*:')    
%         end
%     end
% end
% xlabel('Iteracja')
% ylabel('Wartości elementów K w danej iteracji')
% title('Zmiana elementów macierzy K')
% legend('w:1 k:1', 'w:1 k:2', 'w:2 k:1', 'w:2 k:2');
% hold off

end

labels_x0 = cell(1,3);
for i = 1:3
    labels_x0{i} = sprintf('x0=[%d;%d]', x0(1,i), x0(2,i));
end

figure(1)
legend(labels_x0);

figure(2)
legend(labels_x0);

figure(3)
legend(labels_x0);









%%%%%%% 4 

Rumba = [1,12,144];

x0 = [10;
      15];


for ii=1:3
R = Rumba(ii);
x = zeros(2,N+1);
u = zeros(1,N+1);
x(:,1) = x0;
K(:,:,N+1) = F;
for i=N:-1:1 
    K(:,:,i) = A'*(K(:,:,i+1) - K(:,:,i+1)*B*((R+B'*K(:,:,i+1)*B)^(-1))*B'*K(:,:,i+1))*A + Q;
end

for i=1:N
    S = -(R+B'*K(:,:,i+1)*B)^(-1)*B'*K(:,:,i+1)*A;

    u(i) = S*x(:,i);
    x(:,i+1) = A*x(:,i) + B*u(i);
    if i==20
     S = -(R+B'*K(:,:,i+1)*B)^(-1)*B'*K(:,:,i+1)*A;
     u(i+1) = S*x(:,i+1);
    end
end

J0 = (1/2)*x0'*K(:,:,1)*x0;

figure(4)
plot(0:20,x(1,:),'x','LineStyle','--')
xlabel('Iteracja')
ylabel('Wartość x1 w danej iteracji')
title('x1, x0=[10;15]')
hold on

figure(5)
plot(0:20,x(2,:),'x','LineStyle','--')
xlabel('Iteracja')
ylabel('Wartość x2 w danej iteracji')
title('x2, x0=[10;15]')
hold on

figure(6)
plot(0:20,u,'o','LineStyle','--')
xlabel('Iteracja')
ylabel('Wartość u w danej iteracji')
title('u, x0=[10;15]')
hold on

% figure
% hold on
% for i=1:2
%     for j=1:2
%         Kpom = K(i,j,:);
%         Kpom = Kpom(:);
%         if ~(i==1 && j==2)
%             stairs(0:20,Kpom,'o--')
%         else
%             stairs(0:20,Kpom,'*:')    
%         end
%     end
% end
% xlabel('Iteracja')
% ylabel('Wartości elementów K w danej iteracji')
% title('Zmiana elementów macierzy K')
% legend('w:1 k:1', 'w:1 k:2', 'w:2 k:1', 'w:2 k:2');
% hold off

end

labels_R = cell(1,3);
for i = 1:3
    labels_R{i} = sprintf('R=%d', Rumba(i));
end

figure(4)
legend(labels_R);

figure(5)
legend(labels_R);

figure(6)
legend(labels_R);