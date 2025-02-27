clc, clear, clf
epsSTUD = 54 / 1000;

% Defining data for FE
m = 100;
n = 100;
T = 1;
u_0 = @(x) 0.*x;
f = @(x,t) - epsSTUD ./ 2 .* ( -2.*t + x.^2 -1 );

% Solving FE
[zeta, x_part, t_part] = FE_IE_solver(u_0, f, m, n, T);

% Surf of FEM solution
[X_mesh, T_mesh] = meshgrid(x_part, t_part);
surf(X_mesh,T_mesh, zeta');
xlabel('Space (x)');
ylabel('Time (t)');
zlabel('u_{approx}');
shading interp;
colorbar;
saveas(gcf,'lab5.surf.png')

% Exact solution
u_exact = @(x,t) -epsSTUD ./2 .*(t.*x.^2 -t);

% Surf of error
surf(X_mesh, T_mesh, zeta' - u_exact(X_mesh, T_mesh))
xlabel('Space (x)');
ylabel('Time (t)');
zlabel('Error (u_{approx} - u_{exact})');
shading interp;
colorbar;
saveas(gcf,'lab5.surf.error.png')

% 2d plot at 5 diffrent t
for l = 1:n/5:(n+2)
    t_current = (l-1)*T / n;

    plot(x_part, zeta(:,l))
    hold on

    x = linspace(0,1);
    plot(x,u_exact(x,t_current));
    hold off

    legend('Approx','Exact')
    xlabel('Space (x)')
    ylabel('u')

    filename = "lab5_t" + t_current + ".png";
    saveas(gcf, filename);
end


function [zeta, x_part, t_part] = FE_IE_solver(u_0, f, m, n, T)
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
    
    % Defining hat functions
    phi = @(x,j) ( (x./h - j + 1) .* ( (j-1 <= x./h) & (x./h <= j) ) ) + ...
        ( (-x./h + j + 1) .* ( (j < x./h) & (x./h < j+1) & (j < m+1)  ) );
    
    % Itterating over time
    for l = 2:(n+1)
    
        % Calculate F
        %F = f(x_part', k*l) *h; % Midpoint rule
        J = (0:m)';
        F = integral( @(x) f(x,k*l) .* phi(x,J) , 0,1,'ArrayValued',true);
    
        % calculate b
        b = M*zeta(:, l-1) + k .* F;
    
        % Calculate A
        A = M + k.*S;
    
        % Solve for zeta
        zeta(:,l) = A\b;
    
    end 

    zeta = [zeta ; zeros(1,n+1)];
    x_part = [0, x_part];

end
