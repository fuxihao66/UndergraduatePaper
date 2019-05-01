function [result] = getGaussVectorPolyOneOrder(i,x_0, h,f,x_num)
c = 1.5;
phi = @(x)(exp(-c*x*x));


mat = [phi(0),phi(h),1,(i-1)*h;phi(h),phi(0),1,i*h;1,1,0,0;(i-1)*h,i*h,0,0];
mat_inv = inv(mat);

a1 = mat_inv(1,1);
b1 = mat_inv(1,2);
c1 = mat_inv(2,1);
d1 = mat_inv(2,2);
e1 = mat_inv(3,1);
g1 = mat_inv(3,2);
m1 = mat_inv(4,1);
n1 = mat_inv(4,2);

mat = [phi(0),phi(h),1, i*h;phi(h),phi(0),1,(i+1)*h;1,1,0,0;i*h,(i+1)*h,0,0];
mat_inv = inv(mat);
a2 = mat_inv(1,1);
b2 = mat_inv(1,2);
c2 = mat_inv(2,1);
d2 = mat_inv(2,2);
e2 = mat_inv(3,1);
g2 = mat_inv(3,2);
m2 = mat_inv(4,1);
n2 = mat_inv(4,2);


if i == x_num
    result = integral(@(x)( f(x)*(b1*phi(x-(i-1)*h)+d1*phi(i*h-x)+g1+n1*x) ),(i-1)*h,i*h,'ArrayValued',true);
else
    result =  integral(@(x)( f(x)*(  a2*phi(x-i*h)+c2*phi((i+1)*h-x)+e2+m2*x  ) ),i*h,(i+1)*h,'ArrayValued',true);

    result = result + integral(@(x)( f(x)*(b1*phi(x-(i-1)*h)+d1*phi(i*h-x)+g1+n1*x) ),(i-1)*h,i*h,'ArrayValued',true);
 
end