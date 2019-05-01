function [result] = getGaussMatrixElePolyOneOrder(i, j, x_0, h, p, q, scale)
c = 1.5;
phi = @(x)(exp(-c*x*x));
phi_d = @(x)(exp(-c*x*x)*(-2*c*x));


mat = [phi(0),phi(h),1,(i-1)*h;phi(h),phi(0),1,i*h;1,1,0,0;(i-1)*h,i*h,0,0];
mat_inv = inv(mat);

b = mat_inv(1,2);
d = mat_inv(2,2);
f = mat_inv(3,2);
n = mat_inv(4,2);

mat = [phi(0),phi(h),1, i*h;phi(h),phi(0),1,(i+1)*h;1,1,0,0;i*h,(i+1)*h,0,0];
mat_inv = inv(mat);
a = mat_inv(1,1);
c = mat_inv(2,1);
e = mat_inv(3,1);
m = mat_inv(4,1);


if abs(i-j)>=2
    result = 0;
elseif i == j && i == scale
    
    result = integral(@(x)( p(x)*(b*phi_d(x-(i-1)*h)+d*phi_d(x-i*h)+n)*(b*phi_d(x-(i-1)*h)+d*phi_d(x-i*h)+n) ),(i-1)*h,i*h,'ArrayValued',true); 
    
    result = result + integral(@(x)(q(x)*( b*phi(x-(i-1)*h)+d*phi(i*h-x)+f+n*x )*( b*phi(x-(i-1)*h)+d*phi(i*h-x)+f+n*x )),(i-1)*h,i*h,'ArrayValued',true);
elseif abs(i-j) == 1
    index = min(i,j);

    mat = [phi(0),phi(h),1, index*h;phi(h),phi(0),1,(index+1)*h;1,1,0,0;index*h,(index+1)*h,0,0];
    mat_inv = inv(mat);
    b = mat_inv(1,2);
    d = mat_inv(2,2);
    f = mat_inv(3,2);
    n = mat_inv(4,2);
    a = mat_inv(1,1);
    c = mat_inv(2,1);
    e = mat_inv(3,1);
    m = mat_inv(4,1);
    
    result = integral(@(x)( p(x)*(a*phi_d(x-index*h)+c*phi_d(x-(index+1)*h)+m)*(b*phi_d(x-index*h)+d*phi_d(x-(index+1)*h)+n) ),index*h,(index+1)*h,'ArrayValued',true);
    
    result = result + integral(@(x)(q(x)*(a*phi(x-index*h)+c*phi((index+1)*h-x)+e+m*x)*( b*phi(x-index*h)+d*phi((index+1)*h-x)+f+n*x )),index*h,(index+1)*h,'ArrayValued',true);
else
    
    
    result = integral(@(x)( p(x)*(b*phi_d(x-(i-1)*h)+d*phi_d(x-i*h)+n)*(b*phi_d(x-(i-1)*h)+d*phi_d(x-i*h)+n) ),(i-1)*h,i*h,'ArrayValued',true);
    result = result + integral(@(x)( p(x)*(a*phi_d(x-i*h)+c*phi_d(x-(i+1)*h)+m)*(a*phi_d(x-i*h)+c*phi_d(x-(i+1)*h)+m) ),i*h,(i+1)*h,'ArrayValued',true);
    
    result = result + integral(@(x)(q(x)*( b*phi(x-(i-1)*h)+d*phi(i*h-x)+f+n*x )*( b*phi(x-(i-1)*h)+d*phi(i*h-x)+f+n*x )),(i-1)*h,i*h,'ArrayValued',true);
    result = result + integral(@(x)(q(x)*( a*phi(x-i*h)+c*phi((i+1)*h-x)+e+m*x )*( a*phi(x-i*h)+c*phi((i+1)*h-x)+e+m*x )),i*h,(i+1)*h,'ArrayValued',true);
end