function [result] = getInvMQVectorPoly(i,x_0, h,f,x_num)
c = 0.61; %scale=4
%c = 0.6359;%scale=20
%c = 0.63636; %scale=40
%c=0.6364522; %scale=50
c2 = 1.00;
%c = sqrt(c2);

phi = @(x)(sqrt(c.^2+x*x));


mat = [phi(0),phi(h),1;phi(h),phi(0),1;1,1,0];
mat_inv = inv(mat);

a = mat_inv(1,1);
b = mat_inv(1,2);
c = mat_inv(2,1);
d = mat_inv(2,2);
e = mat_inv(3,1);
g = mat_inv(3,2);


if i == x_num
    result = integral(@(x)( f(x)*(b*phi(x-(i-1)*h)+d*phi(i*h-x)+g) ),(i-1)*h,i*h,'ArrayValued',true);
else
    result =  integral(@(x)( f(x)*(  a*phi(x-i*h)+c*phi((i+1)*h-x)+e  ) ),i*h,(i+1)*h,'ArrayValued',true);

    result = result + integral(@(x)( f(x)*(b*phi(x-(i-1)*h)+d*phi(i*h-x)+g) ),(i-1)*h,i*h,'ArrayValued',true);
 
end