%**************************************************************************
%                          Finite Element Method
%                                Part A.1
%                         Evripidis Sidiropoulos
%**************************************************************************

clc; clear all; close all;

% Declaration of Initial Values

x0 = 0;                              % x coordinate of center
y0 = 0;                              % y coordinate of center
bb = 1.75e-3;                        % External Radius b
aa = exp(log(bb) - 5/6);             % Internal Radius a
V = 1;                               % Voltage

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
%[p, e, t] = refinemesh(d1, p, e, t);
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
        if (radius(1) > (aa+bb)/2)
            X0(n(1)) = 0;
        else
            X0(n(1)) = V;
        end    
        if (radius(2) > (aa+bb)/2)
            X0(n(2)) = 0;
        else
            X0(n(2)) = V;
        end
    end
end

for in = 1:Nn
    text(p(1, in), p(2, in), num2str(node_id(in)));
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

S = spalloc(Nf, Nf, 7*Nf);           % Sparse Array
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
    Ae = abs(D)/2;                   % Element area
    for i = 1:3
        for j = 1:3
            Se(i, j) = (b(i)*b(j) + c(i)*c(j))*Ae;
            if (node_id(n(i)) == 1)
                if (node_id(n(j)) == 1)
                    S(index(n(i)), index(n(j))) = S(index(n(i)), index(n(j))) + Se(i, j);
                else
                    B(index(n(i))) = B(index(n(i))) - Se(i, j)*X0(n(j));
                end
            end
        end
    end   
end

tic;
X = S\B;                             % System Solution via Gauss-Jordan Elimination
toc;

for in = 1:Nn
    if (node_id(in) == 1)
        X0(in) = X(index(in));
    end
end

% Plot of Solution and Electric Field E

figure;                              % Plot of Solution
pdeplot(p, e, t, 'xydata', X0, 'contour', 'off', 'mesh', 'on');
axis equal;
axis tight;
hold on;
colormap jet;

[ux, uy] = pdegrad(p, t, X0);        % Computation of Electric Field
figure;                              % Plot of Electric Field
pdeplot(p, e, t, 'xydata', X0, 'FlowData', [-ux; -uy]);
axis equal;
axis tight;
hold on;
colormap jet;

figure;                              % Plot of Sparse Array S
spy(S);
hold on;

% System Solution via Generalized Minimum Residual Method (GMRES)

tic;
X = gmres(S, B);
toc;

% System Solution via Biconjugate Gradient Method

tic;
X = bicg(S, B);
toc;

% Calculation of Energy and Capacitance

E = -[ux, uy];
W = norm(E*E')*Ae/2;
C = 2*W/(V^2);
Cexact = 2*pi/log(bb/aa);
Wexact = Cexact*(V^2)/2;