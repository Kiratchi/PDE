clc, clear, clf
epsSTUD = 54 / 1000;

m = 3;
n = 5;
T = 1;

u_0 = @(x) 0.*x;
f = @(x,t) - epsSTUD ./ 2 .* ( 2.*t + x.^2 -1 );

% Defining time partition
k = T / n;
t_part = 0:k:T;

% Defining spacial partition
h = 1 / (m + 1);
x_part = (0+h):h:1;



% Defining Mass matrix
M = 4*diag(ones(1,m)) + diag(ones(1,m-1),1) + diag(ones(1,m-1),-1);
M = M * h/6;

% Defining Stiffness matrix
S = 2*diag(ones(1,m)) - diag(ones(1,m-1),1) - diag(ones(1,m-1),-1);
S = S/h;


zeta = zeros(n,m+1);

% Calculating zeta_0
zeta_0 = u_0(x_part);
zeta(1,:) = zeta_0;


phi = @(x,j) ( (x./h - j + 1) .* ( (j-1 <= x./h) & (x./h <= j) ) ) + ...
               ( (-x./h + j + 1) .* ( (j < x./h) & (x./h < j+1) ) );

% Kan vara bugg här
%hold on
%for j = 0:m
%    plot(x_part, phi(x_part, j))
%end
%hold off



% Begin loop
% At a specified "l"
l = 2

% Calculate F

J = 0:m
F = integral( @(x) f(x,k*l) .* phi(x,J) , 0,1,'ArrayValued',true)

% calculate b
b = M.*zeta(l-1, 1:m) + k .* F

