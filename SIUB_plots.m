subplot(2,2,1) % S
hold on
grid on
box on 
plot(t/T,u(:,1),color,'linewidth',2);
xlabel('$\bf Time\;(years)$','Interpreter','latex')
xticks([1:kmax]);
title('$\bf Healthy \; leaves \;S $ ','Interpreter','latex')
    
subplot(2,2,2) % I
hold on
grid on
box on 
p=plot(t/T,u(:,2),color,'linewidth',2);
xlabel('$\bf Time\;(years)$','Interpreter','latex')
xticks([1:kmax]);
title('$\bf Infected \; leaves \;I$','Interpreter','latex'); 

subplot(2,2,3) % U
hold on
grid on
box on 
plot(t/T,u(:,3),color,'linewidth',2);
xlabel('$\bf Time\;(years)$','Interpreter','latex')
xticks([1:kmax]);
title('$\bf Uredospores \; U$','Interpreter','latex')

subplot(2,2,4) % B
grid on
box on 
hold on
plot(t/T,u(:,6),color,'linewidth',2);
xlabel('$\bf Time\;(years)$','Interpreter','latex')
xticks([1:kmax]);
title('$\bf Berries \; B$','Interpreter','latex' )

