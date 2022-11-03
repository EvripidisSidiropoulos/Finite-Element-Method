%**************************************************************************
%                          Finite Element Method
%                                Part B.1
%                                TM Modes
%                         Evripidis Sidiropoulos
%**************************************************************************

clc; clear all; close all;

% Declaration of Initial Values

aa = 1e-2;                           % Radius of Cylindrical Waveguide
x0 = 0;                              % x coordinate of center
y0 = 0;                              % y coordinate of center
k = 12;                              % Number of modes

% Construction and plotting of Geometry

gd = [1;                             % Geometry Description Matrix
      x0;                         
      y0;
      aa];
  
d1 = decsg(gd);                      % Decomposed Solid Geometry Implementation
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
        if (radius(1) == aa)
            if (radius(2) == aa)
                X0(n(1)) = 0;
                X0(n(2)) = 0;
            end
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
for in = 1:Nn
    text(p(1, in), p(2, in), num2str(index(in)));
end

% Initialization of Sparse Array

S = spalloc(Nf, Nf, 7*Nf);           % Stiffness Sparse Array
T = spalloc(Nf, Nf, 7*Nf);           % Mass Sparse Array

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
    Ae = abs(D)/2;                   % Element area
    for i = 1:3
        for j = 1:3
            Se(i, j) = (b(i)*b(j) + c(i)*c(j))*Ae;
            if (i == j)
                Te(i, j) = Ae/6;
            else
                Te(i, j) = Ae/12;
            end    
            if (node_id(n(i)) == 1)
                if (node_id(n(j)) == 1)
                    S(index(n(i)), index(n(j))) = S(index(n(i)), index(n(j))) + Se(i, j);
                    T(index(n(i)), index(n(j))) = T(index(n(i)), index(n(j))) + Te(i, j);
                end
            end
        end
    end   
end

[V, D] = eigs(S, T, k, 0);           % Solution of Generalised Eigenvalue
                                     % Problem of TM Modes

for i = 1:k                          % TM Modes Solution and Plot
    X1 = V(:,i);
    for in = 1:Nn
        if (node_id(in) == 1)
            X0(in) = X1(index(in));
        end
    end
    fc(i) = sqrt(D(i,i))*(3e8/(2*pi));
    figure;                              
    pdeplot(p, e, t, 'xydata', X0, 'contour', 'on', 'mesh', 'off');
    axis equal;
    axis tight;
    hold on;
    colormap jet;
    disp(fc(i));
end

figure;                              % Plot of Sparse Array S
spy(S);
hold on;