function [] = get_order()
h = zeros(20, 1);
error = zeros(20, 1);
grad_error = zeros(20, 1);
for M = 20:40

    h(M-19) = log(1/(M+1));
    [error_new, grad_error_new] = fem2(M);
    error(M-19) = error_new;
    grad_error(M-19) = grad_error_new;

end


%disp(size(h));
%disp(size(error));
plot(h, error, 'b');
hold on;
plot(h,  grad_error, 'r');

title('Log error vs Log h plot');
xlabel('log(h)');
ylabel('log(error)');
legend('Solution error plot', 'Gradient error plot');
xlim([-4, -2.8]);
ylim([-30 7]);
p1 = polyfit(h, error, 1);
p2 = polyfit(h, grad_error, 1);

disp(p1(1));
disp(p2(1));

end