function [error] = RBFErrorL2_nonpoly(x_array, y_array, target_func, phi)   % calcalate L^2 error
%c2 = 1.00;
%c = sqrt(c2);

%phi = @(x)(sqrt(c.^2+x.*x));

[row_num,col_num] = size(x_array);
error = 0.0;
s = 1;
if row_num == 1
    s = col_num;
elseif col_num == 1
    s = row_num;
end

h = x_array(2) - x_array(1);
mat = [phi(0),phi(h);phi(h),phi(0)];


for i = 1:(s-1)
    para = mat\([y_array(i), y_array(i+1)]');
    func = @(x)(para(1).*phi(x-x_array(i)) + para(2).*phi(x-x_array(i+1)));  %linear part
    minus = @(x)((func(x)-target_func(x)).^2);  % difference between linear function and accurate solution
    error = error + integral(minus, x_array(i), x_array(i+1));
end

error = sqrt(error);