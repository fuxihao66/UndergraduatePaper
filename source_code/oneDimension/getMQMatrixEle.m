function [result] = getMQMatrixEle(i,j,h,p,q)
c = 900;
phi = @(x)(sqrt(c.^2+x*x));
phi_d = @(x)(x/sqrt(c.^2+x*x));
mat = [phi(0),phi(h);phi(h),phi(0)];
mat_inv = inv(mat);

a = mat_inv(1,1);
b = mat_inv(1,2);
c = mat_inv(2,1);
d = mat_inv(2,2);
if abs(i-j)>=2
    result = 0;
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