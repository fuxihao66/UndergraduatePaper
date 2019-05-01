function [error] = ErrorCal(x_array, y_array, target_func)   % calcalate L^2 error for linear

[row_num,col_num] = size(x_array);
error = 0.0;
s = 1;
if row_num == 1
    s = col_num;
elseif col_num == 1
    s = row_num;
end

for i = 1:(s-1)
    func = @(x)((y_array(i+1)-y_array(i))/(x_array(i+1)-x_array(i)).*(x-x_array(i))+y_array(i));  %linear part
    minus = @(x)((func(x)-target_func(x)).^2);  % difference between linear function and accurate solution
    error = error + integral(minus, x_array(i), x_array(i+1));
end

error = sqrt(error);