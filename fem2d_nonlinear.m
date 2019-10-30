function [] = fem2d_nonlinear(m, n)

% create a point grid of mn =: N points

[x,y]=ndgrid((0:m-1)/(m-1),(0:n-1)/(n-1)); 
p=[x(:),y(:)]; % N by 2 matrix listing x,y coordinates of all N=mn nodes
T = 2*(m-1)*(n-1);
% generate mesh of T=2(m-1)(n-1) right triangles in unit square
% includes boundary nodes, mesh spacing 1/(m-1) and 1/(n-1)
%t=[1,2,m+2; 1,m+2,m+1]; % 3 node numbers for two triangles in first square
t = zeros(T, 3);
count = 1;
for ii = 1:n-1
    for jj = 1:m-1
        t([count, count + 1], :) = [(m*(ii-1) + jj) , (m*(ii-1) + jj + 1) , (m*(ii) + jj +1) ; (m*(ii-1) + jj), (m*(ii) + jj + 1), (m*(ii) + jj) ]; 
        count = count + 2;
    end
end

% final t lists 3 node numbers of all triangles in T by 3 matrix 

%Node numbers for all boundary nodes
b = [1:m,m+1:m:m*n,2*m:m:m*n,m*n-m+2:m*n-1]; % bottom, left, right, top  

N = size(p,1); T = size(t,1); % number of nodes, number of triangles
% p lists x,y coordinates of N nodes, t lists triangles by 3 node numbers
K = sparse(N,N); % Declare stiffness matrix as sparse zero matrix 
F = zeros(N,1); % Load vector initialized 

U_init = zeros(N, 1); %initial coefficient vector intialized as zero

for ii = 1:20
for e = 1:T  % integration over one triangular element at a time
  nodes = t(e,:); % row of t = node numbers of the 3 corners of triangle e
  Pe = [ones(3,1), p(nodes,:)]; % 3 by 3 matrix with rows=[1 xcorner ycorner] 
  Area = abs(det(Pe))/2; % area of triangle e = half of parallelogram area
  C = inv(Pe); % columns of C are coeffs in a+bx+cy to give phi=1,0,0 at nodes
  % now compute 3 by 3 Ke and 3 by 1 Fe for element e
  grad = C(2:3,:); % element matrix from slopes b,c in grad
  Ke = 10*Area*grad'*grad + grad'*grad * (U_init(1) * quad_f(2, p(nodes(1), :), p(nodes(2),: ), p(nodes(3), :), C(:, 1)) + U_init(2) * quad_f(2, p(nodes(1), :), p(nodes(2),: ), p(nodes(3), :), C(:, 2)) + U_init(3) * quad_f(2, p(nodes(1), :), p(nodes(2),: ), p(nodes(3), :), C(:, 3)) )^2;
  Fe = [quad_f(2, p(nodes(1), :), p(nodes(2),: ), p(nodes(3), :), C(:, 1)); quad_f(2, p(nodes(1), :), p(nodes(2),: ), p(nodes(3), :), C(:, 2)); quad_f(2, p(nodes(1), :), p(nodes(2),: ), p(nodes(3), :), C(:, 3))];
  %Fe=Area/3; % integral of phi over triangle is volume of pyramid: f(x,y)=1
  % multiply Fe by f at centroid for load f(x,y): one-point quadrature!
  % centroid would be mean(p(nodes,:)) = average of 3 node coordinates
  K(nodes,nodes) = K(nodes,nodes)+Ke; % add Ke to 9 entries of global K
  F(nodes) = F(nodes)+Fe; % add Fe to 3 components of load vector F
end   % all T element matrices and vectors now assembled into K and F


% Put boundary conditions U(b)=0 at nodes numbers in list b
K(b,:) = 0; K(:,b) = 0; F(b) = 0; % put zeros in boundary rows/columns of K and F
%To remove singularity of K, we put I as a submatrix  
K(b,b) = speye(length(b),length(b));
Kb = K; Fb = F; % Assign Stiffness matrix Kb  and load vector Fb

% Solving for the vector U , putting identity at (b, b) will produce U(b)=0 at boundary nodes
U = Kb\Fb;
disp(norm(U - U_init));
U_init = U; 

end

% The approximation to solution is U_1 phi_1 + ... + U_N phi_N
%Plot all triangles
trisurf(t, p(:,1), p(:,2), U);

end