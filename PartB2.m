%**************************************************************************
%                          Finite Element Method
%                                Part B.2
%                         Evripidis Sidiropoulos
%**************************************************************************

clc; clear all; close all;

% Declaration of Initial Values

x0 = 0;                              % x coordinate of center
y0 = 0;                              % y coordinate of center
E0 = 1;                              % Amplitude of Electric Field
f0 = 3e8;                            % Frequency
wl = 1;                              % Wavelength
k0 = 2*pi/wl;                        % Wave number
aa = 5*wl/2;                         % Scatterer's Radius
bb = 9*wl/2;                         % ABC Radius
Et = @(x)E0*exp(-1i*k0*x);           % Total Field

% Construction and plotting of Geometry

gd = [1 1;                           % Geometry Description Matrix
      x0 x0;                         
      y0 y0;
      bb aa];
sf = 'R2-R1';                        % Set Formula for subtraction of inner circle
ns = [82 82; 50 49];                 % Name-space Matrix, used with sf
                                     % as inputs to decsg alongside gd

d1 = decsg(gd, sf, ns);              % Decomposed Solid Geometry Implementation
[p, e, t] = initmesh(d1);            % Creation of 2-D Triangular Mesh
[p, e, t] = refinemesh(d1, p, e, t); % Refinement of Triangular Mesh x2
[p, e, t] = refinemesh(d1, p, e, t);
figure;                              % Plot of Mesh
pdeplot(p, e, t);
axis equal;
axis tight;
hold on;

% Definition of known numbering using the node_id vector

Nn = size(p, 2);                     % Number of Nodes
Ne = size(t, 2);                     % Number of elements (triangles)
Nd = size(e, 2);                     % Number of (boundary) edges

node_id = ones(Nn, 1);               % Initialization of node flag
X0 = zeros(Nn, 1);                   % Voltage values at every node
E = zeros(Nn, 1);

for id = 1:Nd                        % Construction of node_id vector
    n(1:2) = e(1:2, id);
    x(1:2) = p(1, n(1:2));
    y(1:2) = p(2, n(1:2));
    r1 = e(6, id);
    r2 = e(7, id);
    if (r1 == 0 || r2==0)            % Assignment of known voltage values
        node_id(n(1)) = 0;
        node_id(n(2)) = 0;
        radius(1:2) = sqrt(x(1:2).^2 + y(1:2).^2);
        if (radius(1) == bb && radius(2) == bb)
            node_id((n(1))) = 2;
            node_id((n(2))) = 2;
            E(n(1)) = -E0;
            E(n(2)) = -E0;
        end
        if (radius(1) == aa && radius(2) == aa)
            node_id((n(1))) = 0;
            node_id((n(2))) = 0;
            E(n(1)) = E0;
            E(n(2)) = E0;
        end    
    end
end

% Definition of unknown numbering using the index vector

ic = 0;                              % Define counter to count unknowns
index = zeros(Nn, 1);                % Define index vector with unknown's 
                                     % numbering for each node

for in = 1:Nn
    if (node_id(in) == 1)
        ic = ic + 1;
        index(in) = ic;
    end
end

Nf = ic;                             % Total number of unknowns

% Initialization of Sparse Array

S = spalloc(Nf, Nf, 7*Nf);           % Stiffness Sparse Array
T = spalloc(Nf, Nf, 7*Nf);           % Mass Sparse Array
A = spalloc(Nf, Nf, 7*Nf);           % Completed Array of the problem
B = zeros(Nf, 1);

for ie = 1:Ne                        % Scan all elements
    n(1:3) = t(1:3, ie);             % n(1), n(2), n(3) are the three nodes of element ie  
    rg = t(4, ie);                   % Element region
    x(1:3) = p(1, n(1:3));           % x Simplex Coordinate
    y(1:3) = p(2, n(1:3));           % y Simplex Coordinate
    D = det([1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3)]); % Simplex to Cartesian
    b(1) = ((y(2) - y(3))/D);
    b(2) = ((y(3) - y(1))/D);
    b(3) = ((y(1) - y(2))/D);
    c(1) = ((x(3) - x(2))/D);
    c(2) = ((x(1) - x(3))/D);
    c(3) = ((x(2) - x(1))/D);
    Aeel = abs(D)/2;                 % Element area
    for i = 1:3
        for j = 1:3
            Se(i, j) = (b(i)*b(j) + c(i)*c(j))*Aeel;
            if (i == j)
                Te(i, j) = Aeel/6;
            else
                Te(i, j) = Aeel/12;
            end
            Ae(i, j) = Se(i, j) - ((2*pi*f0)^2)*Te(i, j);
            if (node_id(n(i)) == 1)
                if (node_id(n(j)) == 1)
                    S(index(n(i)), index(n(j))) = S(index(n(i)), index(n(j))) + Se(i, j);
                    T(index(n(i)), index(n(j))) = T(index(n(i)), index(n(j))) + Te(i, j);
                    A(index(n(i)), index(n(j))) = A(index(n(i)), index(n(j))) + Ae(i, j);
                else
                    if (i == j)
                        Te(i, j) = 1/3;
                        Se(i, j) = (bb^2)/(2*pi*bb);
                    else
                        Te(i, j) = 1/6;
                        Se(i, j) = -(bb^2)/(2*pi*bb);
                    end
                    B(index(n(i))) = B(index(n(i))) - Ae(i, j)*E(n(j)) - Se(i, j)*E(n(j))...
                        + Te(i, j)*E(n(j));
                end
            end
        end
    end   
end

X = S\B;                             % System Solution via Gauss-Jordan Elimination

for in = 1:Nn
    if (node_id(in) == 1)
        E(in) = X(index(in));
    end
end

figure;                              
pdeplot(p, e, t, 'xydata', E, 'contour', 'on', 'mesh', 'off');
axis equal;
axis tight;
hold on;
colormap jet;