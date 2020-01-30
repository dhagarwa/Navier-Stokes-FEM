function [] = fem(M)
%Solving BVP -u''(x) = f(x), u(0) = u(1) = 0; 0 <= x <= 1
%by the finite element method
%Initial space is V. We search for solution in the Vh which is the space of
%continuous function composed of piecewise linear polynomials. 
%The finite element is  (K, S, P) where K = subinterval of length h
%S is the value of function at endpoints of each subinterval. P is the
%space of piecewise linear polynomials. S is clearly P-unisolvent.

%We for the linear system associated Ak = b
%A is associated stiffness matrix A(i, j) = (h(i), h(j)) where h(i) is hat function




%Stiffness matrix entries
A = zeros(M);
for ii = 1:M
    for jj = 1:M
        if ii == jj
            A(ii, jj) = 2;
        elseif ii - jj == 1 | jj - ii == 1
            A(ii, jj) = -1;
        end
    end
end


%b vector
b = zeros(M, 1);
for ii = 1:M
    b(ii) = 1/(M+1);
end

k = (M+1) * A\b;
disp(k);
%Approximate solution = k(1)h(1) + k(2)h(2) + ..+ k(M)h(M)
%Actual solution is 1/2(x - x^2)

%PLotting both actual and approximate solution
points_x = linspace(0, 1, 1000);
points_actual_y = 0.5*(points_x - points_x.^2);
len = length(points_x);
points_approx_y = zeros(len, 1);

for ii=1: len
    subinterval_index = ceil(points_x(ii)*(M+1));
    if subinterval_index == 0
        subinterval_index = 1;
    end
    
    if subinterval_index == 1
        points_approx_y(ii) = k(subinterval_index) * points_x(ii)*(M+1);
    elseif subinterval_index == M+1
        points_approx_y(ii) = k(subinterval_index-1)* (1 - points_x(ii))*(M+1);
        
    else
        points_approx_y(ii) = k(subinterval_index -1)*((subinterval_index)/(M+1) - points_x(ii))*(M+1) + k(subinterval_index)*(points_x(ii) - (subinterval_index-1)/(M+1))*(M+1);
    
    end   
end

plot(points_x, points_actual_y, 'b');
hold on;
plot(points_x, points_approx_y, 'r');


end