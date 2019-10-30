function [error, grad_error] = fem(M)
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
    b(ii) = ii/((M+1)^2);
end

k = (M+1) * A\b;
%disp(k);
%Approximate solution = k(1)h(1) + k(2)h(2) + ..+ k(M)h(M)
%Actual solution is -x^3/6 + x/6

%PLotting both actual and approximate solution
points_x = linspace(0, 1, 1000);
points_actual_y = (1/6)*(points_x - points_x.^3);
len = length(points_x);
points_approx_y = zeros(len, 1);
new_error = 0;
error = 0;


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
    
    new_error = abs(points_approx_y(ii) - points_actual_y(ii));
    
    if new_error > error
        error = new_error;
    end
end

sup_norm_u = 0.06415;
h = 1/(M+1);
sup_norm_gradu = -h^2/2 + 1/6;

plot(points_x, points_actual_y, 'b');
hold on;
plot(points_x, points_approx_y, 'r');
%Gradient of actual solution is -x^2/2 + 1/6
%Recovered gradient is n th interval is -k(n-1)* M+1 + k(n)*M+1 except end
%intervals where only one term is present. At node, gradient sampled is
%taken to be average of gradient at centroids of neighbouring two elements.
%hence gradient is obtained by interpolation these values. 

sampled_gradient_centroids = zeros(M+1, 1); %Store sampled gradient at centroid of each interval

for ii = 1: M+1

    if ii ==1
        sampled_gradient_centroids(ii) = k(ii)*(M+1);
        
    elseif ii == M+1
        sampled_gradient_centroids(ii) = -k(ii-1)*(M+1);
        
    else
        sampled_gradient_centroids(ii) = k(ii)*(M+1) - k(ii-1)*(M+1);
        
    end
    
end


%disp('Sampled Gradient at Centroids');
%disp(sampled_gradient_centroids);

approx_gradient = zeros(M, 1);  %To store the approximation to gradient at nodes , these are average of centroid values

for ii= 1: M
    approx_gradient(ii) = (0.5)*(sampled_gradient_centroids(ii) + sampled_gradient_centroids(ii+1)); 
    
end

disp('Approx_gradient:');
disp(approx_gradient);
points_x = linspace(1/(M+1), 1 - 1/(M+1), 1000);
points_actual_y = (1/6)*(1 - 3*points_x.^2);
len = length(points_x);
points_approx_y = zeros(len, 1);
grad_error = 0;
new_grad_error =0;
for ii=1: len
    subinterval_index = ceil(points_x(ii)*(M+1));
    if subinterval_index == 0
        subinterval_index = 1;
    end
    
    if subinterval_index == 1
        points_approx_y(ii) = approx_gradient(subinterval_index) * points_x(ii)*(M+1);
    elseif subinterval_index == M+1
        points_approx_y(ii) = approx_gradient(subinterval_index-1)* (1 - points_x(ii))*(M+1);
        
    else
        points_approx_y(ii) = approx_gradient(subinterval_index -1)*((subinterval_index)/(M+1) - points_x(ii))*(M+1) + approx_gradient(subinterval_index)*(points_x(ii) - (subinterval_index-1)/(M+1))*(M+1);
        new_grad_error = abs(points_approx_y(ii) - points_actual_y(ii));
    end 
    
     if new_grad_error > grad_error
         grad_error = new_grad_error;
     end
     
end

plot(points_x, points_actual_y, '--b');
hold on;
plot(points_x, points_approx_y,  '--r');

title('Apprximate solution and gradient vs true solution and gradient');
xlabel('x');
ylabel('values');
legend('True solution', 'Finite element solution', 'True gradient', 'Recovered gradient');
error = log((error)/ (sup_norm_u + sup_norm_gradu));
%disp('Error:');
%disp(error);
grad_error = log((grad_error)/(sup_norm_gradu + sup_norm_u));
%disp('Grad Error: ');
%disp(grad_error);

end