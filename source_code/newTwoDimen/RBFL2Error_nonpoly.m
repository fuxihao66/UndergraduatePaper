function [error] = RBFL2Error_nonpoly(x_array, y_array, numerical_solu, target_func, phi)   % calcalate L^2 error
% x_array = [0, 0.2, 0.4, ..., 3]
% y_array = [-1, -0.5, ..., 1]
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
h2 = y_array(2) - y_array(1);

for i = 1:N1-1
    for j = 1:N2-1
        u1 = numerical_solu(j,i); %= u( x_array(i),y_array(j) )
        u2 = numerical_solu(j+1,i);
        u3 = numerical_solu(j,i+1);
        u4 = numerical_solu(j+1,i+1);
        mat = [phi(0,0),phi(h2,0),phi(h1,0),phi(sqrt(h1.*h1+h2.*h2),0);phi(h2,0),phi(0,0),phi(sqrt(h1.*h1+h2.*h2),0),phi(h1,0); ...
               phi(h1,0),phi(sqrt(h1.*h1+h2.*h2),0),phi(0,0),phi(h2,0);phi(sqrt(h1.*h1+h2.*h2),0),phi(h1,0),phi(h2,0),phi(0,0)];
        lambda = mat\([u1,u2,u3,u4]');
        func = @(x,y)( lambda(1).*phi( x-x_array(i),y-y_array(j) )+lambda(2).*phi( x-x_array(i),y-y_array(j+1) ) ...
                      +lambda(3).*phi( x-x_array(i+1),y-y_array(j) )+lambda(4).*phi( x-x_array(i+1),y-y_array(j+1) ) );
        minus = @(x,y)((func(x,y)-target_func(x,y)).^2);
        error = error + integral2(minus, x_array(i), x_array(i+1), y_array(j), y_array(j+1));
    end
end
error = sqrt(error);

