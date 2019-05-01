function [result] = getGaussMatrixEle(i, j, x_0, h, p, q, scale, c)
%c = 0.914;
%c = 1.2;   % best for L2
%c = 1.212;
%c = 0.1;
phi = @(x)(exp(-c*x*x));
phi_d = @(x)(exp(-c*x*x)*(-2*c*x));

mat = [phi(0),phi(h);phi(h),phi(0)];
mat_inv = inv(mat);

a = mat_inv(1,1);
b = mat_inv(1,2);
c = mat_inv(2,1);
d = mat_inv(2,2);


if abs(i-j)>=2
    result = 0;
elseif i == j && i == scale
    result = integral(@(x)( p(x)*(b*phi_d(x-(i-1)*h)+d*phi_d(x-i*h))*(b*phi_d(x-(i-1)*h)+d*phi_d(x-i*h)) ),(i-1)*h,i*h,'ArrayValued',true); 
    
    result = result + integral(@(x)(q(x)*( b*phi(x-(i-1)*h)+d*phi(i*h-x) )*( b*phi(x-(i-1)*h)+d*phi(i*h-x) )),(i-1)*h,i*h,'ArrayValued',true);
elseif abs(i-j) == 1
    index = min(i,j);
    result = integral(@(x)( p(x)*(a*phi_d(x-index*h)+c*phi_d(x-(index+1)*h))*(b*phi_d(x-index*h)+d*phi_d(x-(index+1)*h)) ),index*h,(index+1)*h,'ArrayValued',true);
    
    result = result + integral(@(x)(q(x)*(a*phi(x-index*h)+c*phi((index+1)*h-x))*( b*phi(x-index*h)+d*phi((index+1)*h-x) )),index*h,(index+1)*h,'ArrayValued',true);
else
    result = integral(@(x)( p(x)*(b*phi_d(x-(i-1)*h)+d*phi_d(x-i*h))*(b*phi_d(x-(i-1)*h)+d*phi_d(x-i*h)) ),(i-1)*h,i*h,'ArrayValued',true);
    result = result + integral(@(x)( p(x)*(a*phi_d(x-i*h)+c*phi_d(x-(i+1)*h))*(a*phi_d(x-i*h)+c*phi_d(x-(i+1)*h)) ),i*h,(i+1)*h,'ArrayValued',true);
    
    result = result + integral(@(x)(q(x)*( b*phi(x-(i-1)*h)+d*phi(i*h-x) )*( b*phi(x-(i-1)*h)+d*phi(i*h-x) )),(i-1)*h,i*h,'ArrayValued',true);
    result = result + integral(@(x)(q(x)*( a*phi(x-i*h)+c*phi((i+1)*h-x) )*( a*phi(x-i*h)+c*phi((i+1)*h-x) )),i*h,(i+1)*h,'ArrayValued',true);
end