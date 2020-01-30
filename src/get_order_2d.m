function [] = get_order_2d()
h = zeros(20, 1);
error = zeros(20, 1);

for M = 20:40

    h(M-19) = log(1/(M+1));
    error_new = fem2d(M, M);
    error(M-19) = error_new;
    

end


%disp(size(h));
%disp(size(error));
plot(h, error, 'b');


title('Log error vs Log h plot');
xlabel('log(h)');
ylabel('log(error)');
legend('Solution error plot');

p1 = polyfit(h, error, 1);


disp(p1(1));


end