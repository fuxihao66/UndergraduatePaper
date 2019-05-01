function [error] = RBF_ERROR_varyC(x_array, y_array, target_func, pphi, epsilon)   % calcalate L^2 error
%c2 = 0.01;
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

h = x_array(3) - x_array(1);


for i = 1:(s-1)/2
    phi = @(x)(pphi(x,epsilon(i)));
    mat = [phi(0),phi(h/2),phi(h),1;phi(h/2),phi(0),phi(h/2),1;phi(h), phi(h/2), phi(0),1;1,1,1,0];
    
    index = 2*i-1;
    
    para = mat\([y_array(index), y_array(index+1), y_array(index+2), 0]');
    
    func = @(x)(para(1).*phi(x-x_array(index)) + para(2).*phi(x-x_array(index+1)) + para(3).*phi(x-x_array(index+2)) + para(4) );  %linear part
    minus = @(x)((func(x)-target_func(x)).^2);  % difference between linear function and accurate solution
    error = error + integral(minus, x_array(index), x_array(index+2));
end

error = sqrt(error);