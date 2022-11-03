%**************************************************************************
%                          Finite Element Method
%                                Part A.2
%                         Evripidis Sidiropoulos
%**************************************************************************

clc; clear all; close all;

% Declaration of Initial Values

x0 = 0;                               % x coordinate of center
y0 = 0;                               % y coordinate of center
w = 4e-2;                             % Capacitor Plate Width
tt = 2e-3;                            % Capacitor Plate Height
d = 1e-2;                             % Capacitor Plates' distance
delta = 1e-5;                         % Margin for Geometry purposes
V = 100;                              % Voltage
er = 2.2;                             % Dielectric Constant

% Construction and plotting of Geometry

gd = [3 3 3 3;                        % Geometry Description Matrix
      4 4 4 4;                         
      -5*w/2 -w/2 -w/2     -w/2;      % First column is the whole area R4
      5*w/2  w/2  w/2      w/2;       % Second column is the Dielectric
      5*w/2  w/2  w/2      w/2;       % Third column is the upper
      -5*w/2 -w/2 -w/2     -w/2;      % Capacitor Plate R2
      -5*w/2 d/2  (d/2)+tt -d/2;      % Fourth column is the lower 
      -5*w/2 d/2  (d/2)+tt -d/2;      % Capacitor Plate R1  
      5*w/2  -d/2 d/2      (-d/2)-tt;    
      5*w/2  -d/2 d/2      (-d/2)-tt]; 
sf = '(R4+R3)-R2-R1';                 % Set Formula for capacitor plates
ns = [82 82 82 82; 52 51 50 49];      % Name-space Matrix, used with sf
                                      % as inputs to decsg alongside gd

d1 = decsg(gd, sf, ns);               % Decomposed Solid Geometry Implementation
[p, e, t] = initmesh(d1);             % Creation of 2-D Triangular Mesh
[p, e, t] = refinemesh(d1, p, e, t);  % Refinement of Triangular Mesh x2
[p, e, t] = refinemesh(d1, p, e, t);
figure;                               % Plot of Mesh
pdeplot(p, e, t);
axis equal;
axis tight;
hold on;

% Definition of known numbering using the node_id vector

Nn = size(p, 2);                      % Number of Nodes
Ne = size(t, 2);                      % Number of elements (triangles)
Nd = size(e, 2);                      % Number of (boundary) edges

node_id = ones(Nn, 1);                % Initialization of node flag
X0 = zeros(Nn, 1);                    % Voltage values at every node

for id = 1:Nd                         % Construction of node_id vector
    n(1:2) = e(1:2, id);
    x(1:2) = p(1, n(1:2));
    y(1:2) = p(2, n(1:2));
    r1 = e(6, id);
    r2 = e(7, id);
    if (r1 == 0 || r2 == 0)           % Assignment of known voltage values
        for i = 1:2
            if ((y(i) < tt+(d/2)+delta) && (y(i) > (d/2)-delta) && (abs(x(i)) < (w/2)+delta))
                X0(n(i)) = V/2;
                node_id(n(i)) = 0;
            end
            if ((y(i) > -tt-(d/2)-delta) && (y(i) < -(d/2)+delta) && (abs(x(i)) < (w/2)+delta))
                X0(n(i)) = -V/2;
                node_id(n(i)) = 0;
            end
        end
    end
end

for in = 1:Nn
    text(p(1, in), p(2, in), num2str(node_id(in)));
end

% Definition of unknown numbering using the index vector

ic = 0;                               % Define counter to count unknowns
index = zeros(Nn, 1);                 % Define index vector with unknown's 
                                      % numbering for each node

for in = 1:Nn
    if (node_id(in) == 1)
        ic = ic + 1;
        index(in) = ic;
    end
end

Nf = ic;                              % Total number of unknowns
for in = 1:Nn
    text(p(1, in), p(2, in), num2str(index(in)));
end

% Initialization of Sparse Array

S = spalloc(Nf, Nf, 7*Nf);            % Sparse Array
B = zeros(Nf, 1);

for ie = 1:Ne                         % Scan all elements
    n(1:3) = t(1:3, ie);              % n(1), n(2), n(3) are the three nodes of element ie  
    rg = t(4, ie);                    % Element region
    x(1:3) = p(1, n(1:3));            % x Simplex Coordinate
    y(1:3) = p(2, n(1:3));            % y Simplex Coordinate
    D = det([1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3)]); % Simplex to Cartesian
    b(1) = ((y(2) - y(3))/D);
    b(2) = ((y(3) - y(1))/D);
    b(3) = ((y(1) - y(2))/D);
    c(1) = ((x(3) - x(2))/D);
    c(2) = ((x(1) - x(3))/D);
    c(3) = ((x(2) - x(1))/D);
    Ae = abs(D)/2;                    % Element area
    for i = 1:3
        for j = 1:3
            Se(i, j) = (b(i)*b(j) + c(i)*c(j))*Ae;
            if (node_id(n(i)) == 1)
                if (node_id(n(j)) == 1)
                    S(index(n(i)), index(n(j))) = S(index(n(i)), index(n(j))) + Se(i, j);
                else
                    B(index(n(i))) = B(index(n(i))) - er*Se(i, j)*X0(n(j));
                end
            end
        end
    end   
end

X = S\B;                              % System Solution via Gauss-Jordan Elimination

for in = 1:Nn
    if (node_id(in) == 1)
        X0(in) = X(index(in));
    end
end

% Plot of Solution and Sparse Array S

figure;                              
pdeplot(p, e, t, 'xydata', X0, 'contour', 'off', 'mesh', 'on');
axis equal;
axis tight;
hold on;
colormap jet;

[ux, uy] = pdegrad(p, t, X0);          % Computation of Electric Field
figure;                                % Plot of Electric Field
pdeplot(p, e, t, 'xydata', X0, 'FlowData', [-ux; -uy]);
axis equal;
axis tight;
hold on;
colormap white;

figure;                              
spy(S);
hold on;

% Calculation of Energy and Capacitance

E = -[ux, uy];
W = er*norm(E*E')*Ae/2;
C = 2*W/(V^2);
Cexact = 2*w*tt*er/d;
Wexact = Cexact*(V^2)/2;