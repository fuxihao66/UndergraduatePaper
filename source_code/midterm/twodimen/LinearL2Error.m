function [error] = LinearL2Error(x_array, y_array, numerical_solu, target_func)   % calcalate L^2 error
% x_array = [0, 0.2, 0.4, ..., 3]
% y_array = [0, -0.2, ..., -1]
% numerical_solu(1,2) = u( x_array(2),y_array(1) )

error = 0.0;
[row_num,col_num] = size(x_array);
N1 = 1;
if row_num == 1
    N1 = col_num;
elseif col_num == 1
    N1 = row_num;
end
[row_num,col_num] = size(y_array);
N2 = 1;
if row_num == 1
    N2 = col_num;
elseif col_num == 1
    N2 = row_num;
end
h1 = x_array(2) - x_array(1);
h2 = y_array(1) - y_array(2);

for i = 1:N1-1
    for j = 1:N2-1
        u1 = numerical_solu(j,i); %= u( x_array(i),y_array(j) )
        u2 = numerical_solu(j+1,i);
        u3 = numerical_solu(j,i+1);
        u4 = numerical_solu(j+1,i+1);
        
        func = @(x,y)( u1.*(1-(x-x_array(i))./(h1)).*(1-(y_array(j)-y)./(h2)) + u2.*(1-(x-x_array(i))./(h1)).*(1-(y-y_array(j+1))./(h2)) ...
                      + u3.*(1-(x_array(i+1)-x)./(h1)).*(1-(y_array(j)-y)./(h2)) + u4.*(1-(x_array(i+1)-x)./(h1)).*(1-(y-y_array(j+1))./(h2)) );
        minus = @(x,y)((func(x,y)-target_func(x,y)).^2);
        error = error + integral2(minus, x_array(i), x_array(i+1), y_array(j+1), y_array(j));
    end
end
error = sqrt(error);

