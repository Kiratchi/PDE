%% Task 1
clf, clear, clc

subplot(2,2,1)
ploting(10,0.01)

subplot(2,2,2)
ploting(100,0.01)

subplot(2,2,3)
ploting(10,1)

subplot(2,2,4)
ploting(20,1)
saveas(gcf,'lab4.task1.png')

function [u_numeric, x_numeric] = nummeric(m,D)
    h= pi / (m+1);

    P = diag(ones(1,m-1),1) - diag(ones(1,m-1),-1);
    P = P/4;

    S = 2*diag(ones(1,m)) - diag(ones(1,m-1),1) - diag(ones(1,m-1),-1);
    S = S * D/h;

    b = h * ones(m,1);
    zeta = (P+S)\b;

    u_numeric = [0; zeta; 0];
    x_numeric = 0:h:pi;
end

function [u_exact, x_exact] = exact(D)
    x_exact = linspace(0, pi);
    u_exact = 2 .* pi ./ (exp(pi ./ (2.*D) ) - 1) .* (1 - exp(x_exact ./ (2.*D) )) +2*x_exact;
end

function ploting(m,D)
    [u_numeric, x_numeric] = nummeric(m,D);
    [u_exact, x_exact] = exact(D);

    plot(x_numeric, u_numeric, x_exact,u_exact)
    legend("Gc(1)", "Exact")
    title(sprintf("D = %.3f, m = %d", D, m));
    hold off
end

%% Task 2 - cG(1)
clc, clear
f = @(x) pi^2 * sin(pi*x);
u =@( x ) sin ( pi * x );

% Showing one solution with m=100
m = 10;
[u_cG1, x_cG1] = cG1(m,f);
plot(x_cG1, u_cG1)
hold on
x = linspace(0,1);
plot(x,sin(pi*x))
hold off
legend("cG(1)", "Exact")
saveas(gcf,'lab4.task2.cG1.png')

max_m = 5;
err = zeros(1,max_m);
for l=1:max_m
    m=2^l;
    h = (1-0) / (m+1);   
    [zeta_cG1, x_cG1] = cG1(m,f);

    x = linspace(0,1,10*m);
    u_cG1 = interp1(x_cG1,zeta_cG1,x ,'linear');
    u_exact = u(x);

    err(l) = norm(u_cG1-u_exact,Inf);
end

hl = 1./ (2.^(1:max_m)+1);
loglog(hl,err,'b*-')
hold on
loglog ( hl , hl .^1 , 'b--' ,hl , hl .^2 , 'r--') 
hold off
xlabel('m')
ylabel('Error')
legend('error','h', 'h^2')
saveas(gcf,'lab4.task2.cG1.error.png')

function [u_cG1, x_cG1] = cG1(m,f)
    h = (1-0) / (m+1);
    x_cG1 = h:h:(1-h);


    % CONSTRUCTING STIFFNESS MATRIX A
    A = 2*diag(ones(1,m)) - diag(ones(1,m-1),1) - diag(ones(1,m-1),-1);
    A = 1/h*A;

    % CONSTRUCTING LOAD VECTOR
    J = 1:m;
    phi = @(x,j) ((x./h  - j + 1) .* ( (j-1 <= x./h) & (x./h <=  j) )) + ...
                 ((-x./h + j + 1) .* ( (j < x./h)    & (x./h < j+1) ));
    b = integral( @(x) f(x) .* phi(x,J) , 0,1,'ArrayValued',true);
    


    % CALCULATING SOLUTION
    u_cG1 = A \b';
    u_cG1 = [0;u_cG1;0];
    x_cG1 = [0,x_cG1,1];
end

%% Task 2 -cG(2)
% Create a program that computes cG(1) and cG(2) FEM approximation
clc, clear
f = @(x) pi^2 * sin(pi*x);
m = 4;

h = (1-0) / m;
x = h/2:h/2:(1-h/2);


% CONSTRUCTING STIFFNESS MATRIX A
% diag
d = repmat([-2, 79], 1, m-1);
d = [d , -2];
A = diag(d);

% diag +1
A = A + diag(8*ones(1,2*m-2), 1) + diag(8*ones(1,2*m-2), -1);

% diag +2
d = repmat([0, 3], 1, m-2);
d = [d , 0];
A = A + diag(d, 2) + diag(d,-2);


A = A / (3*h);

% CONSTRUCTING LOAD VECTOR

% The integer approximations are *2
H = repmat([1 2],1, m-1);
H = [H, 1];

b = h.*H .* f(x);


% CALCULATING SOLUTION
zeta = A \ b';
zeta = [0;zeta;0];
x = [0,x,1];
plot(x,zeta)

plot(x,zeta)
hold on
plot(x,sin(pi*x))
hold off
legend("cG(2)", "Exact")
saveas(gcf,'lab4.task2.cG2.png')








%%
clc; clear; close all;

% Define f(x) (Example: f(x) = sin(pi*x) )
f = @(x) sin(pi*x);

% Define domain parameters
N = 10;         % Number of subintervals
h = 1/N;        % Step size
x = (0:N) * h;  % Node points

% Compute integral for half-integer basis functions
I_half = zeros(1, N-1);
for i = 2:N  % Half-integer indices (i-1/2)
    x_i = i * h;
    x_im1 = (i - 1) * h;
    
    phi_half = @(x) (4/h^2) * (x - x_im1) .* (x_i - x) .* (x >= x_im1 & x <= x_i);
    integrand = @(x) f(x) .* phi_half(x);
    
    I_half(i-1) = integral(integrand, x_im1, x_i);
end

% Compute integral for integer basis functions
I_int = zeros(1, N-1);
for i = 2:N-1
    x_i = i * h;
    x_ip = (i + 1) * h;
    x_im = (i - 1) * h;
    
    x_ihalf_p = (i + 1/2) * h;
    x_ihalf_m = (i - 1/2) * h;
    
    phi_int = @(x) (2/h^2) * ((x - x_ihalf_p) .* (x - x_ip) .* (x >= x_im & x <= x_i)) + ...
                    (2/h^2) * ((x - x_ihalf_m) .* (x - x_im) .* (x >= x_i & x <= x_ip));
    
    integrand = @(x) f(x) .* phi_int(x);
    
    I_int(i) = integral(integrand, x_im, x_ip);
end

% Display results
disp('Integral values for half-integer basis functions:');
disp(I_half);

disp('Integral values for integer basis functions:');
disp(I_int);
