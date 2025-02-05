% POLYNOMIAL INTERPOLATION IN 1D


%% TASK 1 
clf

x = [2,3,4,5,6];
y = [2,6,5,5,6];
%figure(),
plot(x,y,'r*', 'MarkerSize',16)


c = polyfit(x,y,4)

xx = linspace(min(x),max(x),50);
yy = polyval(c,xx)

hold on
plot(xx,yy, 'k')
legend("Data points", "4 deg interpolant");
saveas(gcf,'lab1.task1.png')


%% Task 2
clf
x = [1975,1980,1985,1990];
y_west = [72.8, 74.2, 75.2, 76.4];
y_east = [70.2, 70.2, 70.3, 71.2];

plot(x,y_west,'b*', x,y_east,'r*')
hold on
c_west = polyfit(x,y_west,length(x)-1);
c_east = polyfit(x,y_east,length(x)-1);

% Interpolant
xx = linspace(1965,1995,50);
yy_west = polyval(c_west,xx);
yy_east = polyval(c_east,xx);
plot(xx,yy_west, 'b', xx, yy_east, 'r')


% Adding estimation for 1970
y_west_70 = polyval(c_west,1970)
y_east_70 = polyval(c_east,1970)
plot(1970, y_west_70, 'bo', 1970, y_east_70, 'ro')

% Adding true values for 1970
plot(1970, 71.8, "bpentagram")
plot(1970, 69.6, "rpentagram")

legend("Data west", "Data east", "Interp. west", "Interp. east")
xlim([1965 1995])
saveas(gcf,'lab1.task2.png')

%% Task 3

x = [1975,1980,1985,1990];
y_west = [72.8, 74.2, 75.2, 76.4];
y_east = [70.2, 70.2, 70.3, 71.2];


y_west_70 = MyLagrangeInterpol(x,y_west,1970)
y_easy_70 = MyLagrangeInterpol(x,y_east,1970)

function yy = MyLagrangeInterpol(x,y,xx)
    n = length(x);
    nn = length(xx);
    for i=1:nn
        yy(i)=0;
        for k=1:n
            x_nk = x([1:k-1 k+1:end]);
            yy(i)=yy(i)+y(k) * prod( (xx(i)-x_nk) ./ (x(k)-x_nk ) );
        end
    end
end
