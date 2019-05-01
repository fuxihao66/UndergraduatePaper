function [c] = getOptimizePara(x_array, y_array, target_func, phi, para)   % calcalate L^2 error
%c2 = 1.00;
%c = sqrt(c2);

%phi = @(x)(sqrt(c.^2+x.*x));

[row_num,col_num] = size(x_array);

s = 1;
if row_num == 1
    s = col_num;
elseif col_num == 1
    s = row_num;
end

h = x_array(2) - x_array(1);
mat = [phi(0, para),phi(h, para);phi(h, para),phi(0, para)];


errorFunc = @(a)(0);

for i = 1:(s-1)
    para = mat\([y_array(i), y_array(i+1)]');
    func = @(a, x)(para(1).*phi(a, x-x_array(i)) + para(2).*phi(a, x-x_array(i+1)));  %linear part
    minus = @(a, x)((func(a, x)-target_func(x)).^2);  % difference between linear function and accurate solution
    
    g = @(a) integral(@(x) minus(a,x),x_array(i), x_array(i+1));
    errorFunc = @(a)(errorFunc(a)+ g(a));
    
end
errorFunc(0.5)
errorFunc(0.9999)
c = fminbnd(errorFunc,0.01,1.0);