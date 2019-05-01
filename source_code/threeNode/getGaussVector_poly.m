function [result] = getGaussVector_poly(i, x_0,h, f, scale)

c = 0.8;

phi = @(x)(exp(-c*x*x));



if i == scale*2
    result = integral(mulFunc(f, getShape_poly(phi, phi, 3, (i/2-1)*h, h)),(i/2-1)*h,i/2*h,'ArrayValued',true);
elseif mod(i,2)==0
    result =  integral(mulFunc(f, getShape_poly(phi, phi, 1, i/2*h, h)),i/2*h,(i/2+1)*h,'ArrayValued',true);

    result = result + integral(mulFunc(f, getShape_poly(phi, phi, 3, (i/2-1)*h, h)) ,(i/2-1)*h,i/2*h,'ArrayValued',true);
else
    result =  integral(mulFunc(f, getShape_poly(phi, phi, 2, (i-1)/2*h, h)),(i-1)/2*h,(i-1)/2*h+h,'ArrayValued',true);
end