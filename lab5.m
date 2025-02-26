clc, clear, clf
epsSTUD = 54 / 1000;

m = 4;
n = 100;
T = 1;


u_0 = @(x) 0.*x;
f = @(x,t) - epsSTUD ./ 2 .* ( 2.*t + x.^2 -1 );

u_exact = @(x,t) -epsSTUD ./2 .*(t.*x.^2 -t);


% Defining time partition
k = T / n;
t_part = 0:k:T;

% Defining spacial partition
h = 1 / (m + 1);
x_part = (0+h):h:1;



% Defining Mass matrix
M = 4*diag(ones(1,m+1)) + diag(ones(1,m),1) + diag(ones(1,m),-1);
M = M * h/6;
M(1,1) = M(1,1) / 2; %Neumann BC

% Defining Stiffness matrix
S = 2*diag(ones(1,m+1)) - diag(ones(1,m),1) - diag(ones(1,m),-1);
S = S/h;
S(1,1) = S(1,1) / 2; %Neumann BC

% Defining zeta
zeta = zeros(m+1, n);

% Calculating zeta_0
zeta_0 = u_0(x_part);
zeta(:,1) = zeta_0;


% Defining hat functions   HALF HAT FUNCTION!!!!
phi = @(x,j) ( (x./h - j + 1) .* ( (j-1 <= x./h) & (x./h <= j) ) ) + ...
               ( (-x./h + j + 1) .* ( (j < x./h) & (x./h < j+1) & (j < m+1)  ) );

% Plot the hat functions
% hold on
% x_line = linspace(0-0.1, 1+0.1, 1000);
% for j = 1:(m+1)
%     plot(x_line, phi(x_line, j));
% end
% hold off


% Itterating over time
for l = 2:(n+1)

    % Calculate F
    %F = f(x_part', k*l) * 2*h;
    J = (0:m)';
    F = integral( @(x) f(x,k*l) .* phi(x,J) , 0,1,'ArrayValued',true)


    % calculate b
    b = M*zeta(:, l-1) + k .* F;

    % Calculate A
    A = M + k.*S;
    
    % Solve for zeta
    zeta(:,l) = A\b;

end

zeta = [zeta ; zeros(1,n+1)];
x_part = [0, x_part];
plot(x_part, zeta(:,n))
hold on
x = linspace(0,1);
plot(x,u_exact(x,T));

disp("Done")