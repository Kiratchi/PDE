clc, clear,clf
[N, T, P] = mygrid(2,1);
[N, T, P] = mygridrefinment(N,T,P);
[N, T, P] = mygridrefinment(N,T,P);
plotmygrid(N,T,P) 


t = [1 2;
     4 5; 
     7 1];

element_stiffness_matrix(t)
element_mass_matrix(t)

function [N, T, P] = mygrid(G,w)
epsAir=8.85e-12;
muAir=4*pi*1e-7;
sigmaAir=0;
epsChicken=4.43e-11;
muChicken=4*pi*1e-7;
sigmaChicken=3e-11;

if G == 0 %Unit square
    N = [0 0;
        1 0 ;
        0 1 ;
        1 1 ];
    T = [ 1 2 3 1 0 1 ;
        2 4 3 1 1 0 ];
    P = [muAir*epsAir*ones(1,4)];
elseif G == 1 %Triangle
    N = [0 0;
        1 0;
        0.5 sqrt(3)/2];
    T = [1 2 3 1 1 1];
    P = [muAir*epsAir*ones(1,3)];
elseif G == 2 % Chicken
    N = [.08 .09; .14 .04; .2 .03; .3 .03; .35 .04 % inside
        .4 .1 ; .35 .16; .3 .2 ; .25 .22; .2 .24
        .14 .24; .1 .2 ; .05 .2 ; .03 .22; .02 .19
        .03 .16; .05 .18; .1 .18; .13 .14; .15 .18
        .2 .2 ; .15 .1 ; .2 .1 ; .3 .1 ; .35 .1
        0 0 ; .1 0 ; .2 0 ; .3 0 ; .4 0 % outside
        .5 0 ; .5 .3 ; .4 .3 ; .25 .3 ; .15 .3
        .08 .3 ; 0 .3 ; 0 .19; 0 .1 ; .5 .15];
    T=[13 14 15 0 0 0; 15 16 17 0 0 0; 15 17 13 0 0 0; 17 12 13 0 0 0
        17 18 12 0 0 0; 18 20 12 0 0 0; 18 19 20 0 0 0; 12 20 11 0 0 0
        20 21 11 0 0 0; 21 10 11 0 0 0; 21 9 10 0 0 0; 1 22 19 0 0 0
        22 20 19 0 0 0; 22 23 20 0 0 0; 2 22 1 0 0 0; 2 23 22 0 0 0
        20 23 21 0 0 0; 23 9 21 0 0 0; 2 3 23 0 0 0; 3 24 23 0 0 0
        23 24 9 0 0 0; 3 4 24 0 0 0; 24 8 9 0 0 0; 24 7 8 0 0 0
        24 25 7 0 0 0; 24 5 25 0 0 0; 4 5 24 0 0 0; 5 6 25 0 0 0
        6 7 25 0 0 0; 26 27 1 1 0 0; 27 2 1 0 0 0; 27 28 2 1 0 0
        28 3 2 0 0 0; 28 29 3 1 0 0; 29 4 3 0 0 0; 29 30 4 1 0 0
        30 5 4 0 0 0; 30 6 5 0 0 0; 30 31 6 1 0 0; 31 40 6 1 0 0
        32 33 40 1 0 1; 33 7 6 0 0 0; 33 8 7 0 0 0; 33 34 8 1 0 0
        34 9 8 0 0 0; 34 10 9 0 0 0; 34 35 10 1 0 0; 35 11 10 0 0 0
        35 36 11 1 0 0; 36 12 11 0 0 0; 36 13 12 0 0 0; 36 14 13 0 0 0
        36 37 14 1 0 0; 37 38 14 1 0 0; 38 15 14 0 0 0; 38 16 15 0 0 0
        38 39 16 1 0 0; 39 1 16 0 0 0; 1 17 16 0 0 0; 1 18 17 0 0 0
        1 19 18 0 0 0; 39 26 1 1 0 0; 40 33 6 0 0 0];
    P=[muChicken*(epsChicken+sigmaChicken/i/w)*ones(1,29) ...
        muAir*epsAir*ones(1,34)];
end
end


function plotmygrid(N,T,P)

% Looping over triangles
for i = 1:size(T,1)
    T_i = T(i,:);
    N_i = N(T_i(1:3),:);

    % Draw fill
    patch(N_i([1,2,3,1],1),N_i([1,2,3,1],2), abs(P(i)))

    % Drawing outer border
    edges = {[1,2], [2,3], [3,1]};
    for k = 1:3
        if T_i(k+3) == 1
            line(N_i(edges{k}, 1), N_i(edges{k}, 2), 'LineWidth', 2, 'Color', 'r');
        end
    end


end

x_shift_lable = 0.01;
y_shift_lable = 0.03;

x_range = max(N(:,1)) - min(N(:,1));
y_range = max(N(:,2)) - min(N(:,2));
% Draw numbers
if size(T,1) < 100
    text(N(:,1)+x_range*x_shift_lable, N(:,2)+y_range*y_shift_lable, string(1:size(N,1)))
end

% Adding padding
padding = 0.1;
axis([min(N(:,1))-x_range*padding  max(N(:,1))+x_range*padding ...
    min(N(:,2)) - y_range*padding max(N(:,2)) + y_range*padding])
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';


end

function [Nr,Tr,Pr] = mygridrefinment(N,T,P)
Nr = N;
nn = size(N,1); %nr of nodes
Tr = [];
nt = 0;
Pr = [];

% Each triangle becomes 4
for j=1:size(T,1)
    i = T(j,1:3); %Nodes of triangle j
    n = N(i,:);

    % New nodes
    n(4,:) = ( n(1,:) + n(2,:) ) / 2;
    n(5,:) = ( n(1,:) + n(3,:) ) / 2;
    n(6,:) = ( n(2,:) + n(3,:) ) / 2;
    
    % Checking if new nodes actually are new
    for k = 4:6
        l = find(Nr(:,1)==n(k,1));
        m = find(Nr(l,2)==n(k,2));
        
        if isempty(m)
            nn = nn +1;
            Nr(nn,:) = n(k,:);
            i(k) = nn;
        
        else
            i(k)=l(m);
        end
    end
    Tr(nt+1,:) = [i(1) i(4) i(5) T(j,4) 0 T(j,6)];
    Tr(nt+2,:) = [i(4) i(5) i(6) 0 0 0];
    Tr(nt+3,:) = [i(6) i(4) i(2) 0 T(j,4) T(j,5)];
    Tr(nt+4,:) = [i(6) i(3) i(5) T(j,5) T(j,6) 0];
    
    Pr(nt+1:nt+4) = P(j);
    nt=nt+4;
end

end

function S = element_stiffness_matrix(t)
    x_1 = t(1,1); x_2 = t(2,1); x_3 = t(3,1);
    y_1 = t(1,2); y_2 = t(2,2); y_3 = t(3,2);

    J = [x_2-x_1 x_3-x_1;
        y_2-y_1 y_3-y_1];
    Jd = (-x_1*y_3 -x_2*y_1 + x_2*y_3 + x_1*y_2 + x_3*y_1 - x_3*y_2);
    %Jd = abs(x_1 * (y_2 - y_3) + x_2 * (y_3 - y_1) + x_3 * (y_1 - y_2))
    gradphi_local = [-1 -1;
                      1  0 ;
                      0  1];
    gradphi_global = gradphi_local/J;

    S = (gradphi_global*gradphi_global') * Jd /2;

    
end

function M = element_mass_matrix(t)
    x_1 = t(1,1); x_2 = t(2,1); x_3 = t(3,1);
    y_1 = t(1,2); y_2 = t(2,2); y_3 = t(3,2);
    
    Jd = abs(x_1 * (y_2 - y_3) + x_2 * (y_3 - y_1) + x_3 * (y_1 - y_2));
    
    M = Jd / 24 * [2 1 1;
                   1 2 1;
                   1 1 2];
end

function [u,K,M] = FEHelmhlotz2d(g, N, T, w, P)

 

end