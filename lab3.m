% Lab 3 - Intitial value problem
%% TASK 1
clear all, clc, clf
epsSTUD = 54 / 1000;
p0= 5;
k = 0.1 + epsSTUD;
tExact = 0:0.05:20;
pExact = p0 * exp(k*tExact);
plot(tExact,pExact)
hold on


[tEuler1,pEuler1] = ForwardEuler(0.5,p0,k);
[tEuler2,pEuler2] = ForwardEuler(0.25,p0,k);
[tEuler3,pEuler3] = ForwardEuler(0.1,p0,k);
plot(tEuler1,pEuler1)
plot(tEuler2,pEuler2)
plot(tEuler3,pEuler3)


xlabel('Time','FontSize',15)
ylabel('Nr Parasites','FontSize',15)
legend('Exact solution', 'h=0.5', 'h=0.25', 'h=0.1' )
saveas(gcf,'lab3.task1.png')
hold off

hplot = [0.5 0.25 0.1];
pExact_end = pExact(end);
pEuler_end = [pEuler1(end), pEuler2(end), pEuler3(end)];
Errplot = [abs(pEuler_end-pExact_end)];
loglog(hplot,Errplot,'b-o',hplot,hplot,'--r')
xlabel('Time step sizes', 'FontSize', 15)
ylabel('Errors Euler', 'FontSize', 15)
legend('Euler', 'Reference line slope 1','Location','NorthWest')
saveas(gcf,'lab3.task1.error.png')


function  [tEuler,pEuler] = ForwardEuler(h,p0,k)
    N = 20 / h + 1;
    tEuler =  0:h:20;
    pEuler = zeros(1,length(tEuler));
    pEuler(1) = p0;
    for n=2:N
        pEuler(n) = pEuler(n-1) + h*k*pEuler(n-1);
    end
end

%% TASK 2
clf, clc, clear
epsSTUD = 54 / 1000;
alpha = 0.005*epsSTUD;
beta = 0.01;
zeta = 0.2;
%h = 0.65;
h = 0.1;
T_end = 10;


H_0 = 500;
Z_0 = 10;
R_0 = 0;


N = T_end/h+1;
t = 0:h:T_end;
H = zeros(1,length(t));
Z = zeros(1,length(t));
R = zeros(1,length(t));

H(1) = H_0;
Z(1) = Z_0;
R(1) = R_0;

for n=1:(N-1)
    H(n+1) = H(n) - h * beta * H(n) * Z(n);
    Z(n+1) = Z(n) + h * ( (beta - alpha) * H(n) * Z(n) + zeta * R(n));
    R(n+1) = R(n) + h * (alpha * H(n) * Z(n) - zeta * R(n));
end

plot(t,H,t,Z,t,R)
xlabel('Time','FontSize',15)
ylabel('Nr people','FontSize',15)
legend('Human','Zombies', 'Dead Zombies')
saveas(gcf,'lab3.task2.png')

T = H+Z+R;
max(T)-min(T)
plot(t,T)
saveas(gcf,'lab3.task2.totalpopulation.png')










