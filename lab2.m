% Numerical Integration in 1d

%% TASK 1
% Create a function for the trapezoidal rule
clc, clear

fprintf("Results for int_0^1 xe^x/(x+1)^2 \n")
fprintf("%-6s : %6.6f\n", "N=10", mytrapezoidalrule(@(x) x .* exp(x) ./ (x+1).^2, 0, 1, 10));
fprintf("%-6s : %6.6f\n", "N=100", mytrapezoidalrule(@(x) x .* exp(x) ./ (x+1).^2, 0, 1, 100));
fprintf("%-6s : %6.6f\n\n", "Exact", (exp(1)-2)/2);

fprintf("Results for int_0^pi sin(x) \n")
fprintf("%-6s : %6.6f\n", "N=10", mytrapezoidalrule(@(x) sin(x), 0, pi, 10));
fprintf("%-6s : %6.6f\n", "N=100", mytrapezoidalrule(@(x) sin(x), 0, pi, 100));
fprintf("%-6s : %6.6f\n\n", "Exact", 2);


function t = mytrapezoidalrule(f,a,b,N)
    h = (b-a)/N;
    x_interior = (a+h):h:(b-h);
    t = h * ( f(a)/2 +f(b)/2 + sum(f(x_interior)) );
end

%% TASK 2
% Create a function for the midpoint rule
clc, clear

fprintf("Results for int_0^1 xe^x/(x+1)^2 \n")
fprintf("%-6s : %6.6f\n", "N=10", mymidpointrule(@(x) x .* exp(x) ./ (x+1).^2, 0, 1, 10));
fprintf("%-6s : %6.6f\n", "N=100", mymidpointrule(@(x) x .* exp(x) ./ (x+1).^2, 0, 1, 100));
fprintf("%-6s : %6.6f\n\n", "Exact", (exp(1)-2)/2);

fprintf("Results for int_0^pi sin(x) \n")
fprintf("%-6s : %6.6f\n", "N=10", mymidpointrule(@(x) sin(x), 0, pi, 10));
fprintf("%-6s : %6.6f\n", "N=100", mymidpointrule(@(x) sin(x), 0, pi, 100));
fprintf("%-6s : %6.6f\n\n", "Exact", 2);


function t = mymidpointrule(f,a,b,N)
    h = (b-a)/N;
    x = (a+h):h:b;
    t = h * sum(f(x-h/2));
end
%% TASK 3
clc, clear
exact_value = 2;
for n=0:7
    N = 2^n;
    T = mytrapezoidalrule(@(x) sin(x), 0, pi, N);
    err = abs(T- exact_value);
    if n ~= 0
        fprintf("N=%-4d err=%6.6f Reduction factor=%6.6f\n", N, err, err / err_prev)

    end
    err_prev = err;
end

