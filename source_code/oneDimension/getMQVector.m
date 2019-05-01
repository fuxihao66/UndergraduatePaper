function [result] = getMQVector(i,h,f)
c = 900;
phi = @(x)(sqrt(c.^2+x*x));

mat = [phi(0),phi(h);phi(h),phi(0)];
mat_inv = inv(mat);

a = mat_inv(1,1);
b = mat_inv(1,2);
c = mat_inv(2,1);
d = mat_inv(2,2);

result = integral(@(x)( f(x)*(b*phi(x-(i-1)*h)+d*phi(x-i*h)) ),(i-1)*h,i*h,'ArrayValued',true);
result = result + integral(@(x)( f(x)*(  a*phi(x-i*h)+c*phi(x-(i+1)*h)  ) ),i*h,(i+1)*h,'ArrayValued',true);