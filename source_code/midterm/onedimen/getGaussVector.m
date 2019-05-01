function [result] = getGaussVector(i,x_0, h,f,x_num, c)
%c = 0.914;
%c = 1.2;    % best for L2
%c = 1.212;
%c = 0.1;
phi = @(x)(exp(-c*x*x));


mat = [phi(0),phi(h);phi(h),phi(0)];
mat_inv = inv(mat);

a = mat_inv(1,1);
b = mat_inv(1,2);
c = mat_inv(2,1);
d = mat_inv(2,2);


if i == x_num
    result = integral(@(x)( f(x)*(b*phi(x-(i-1)*h)+d*phi(i*h-x)) ),(i-1)*h,i*h,'ArrayValued',true);
else
    result =  integral(@(x)( f(x)*(  a*phi(x-i*h)+c*phi((i+1)*h-x)  ) ),i*h,(i+1)*h,'ArrayValued',true);

    result = result + integral(@(x)( f(x)*(b*phi(x-(i-1)*h)+d*phi(i*h-x)) ),(i-1)*h,i*h,'ArrayValued',true);
 
end