function [result] = getMQVector_poly(i, x_0,h, f, scale, epsilon)

pphi = @(x, c)(sqrt(1+c.*x.*x));



if i == scale*2
    phi = @(x)(pphi(x, epsilon(scale)));
       
    result = integral(mulFunc(f, getShape_poly(phi, phi, 3, (i/2-1)*h, h)),(i/2-1)*h,i/2*h,'ArrayValued',true);
elseif mod(i,2)==0
    phi = @(x)(pphi(x, epsilon(i/2+1)));
    result =  integral(mulFunc(f, getShape_poly(phi, phi, 1, i/2*h, h)),i/2*h,(i/2+1)*h,'ArrayValued',true);
    
    phi = @(x)(pphi(x, epsilon(i/2)));
    result = result + integral(mulFunc(f, getShape_poly(phi, phi, 3, (i/2-1)*h, h)) ,(i/2-1)*h,i/2*h,'ArrayValued',true);
else
    phi = @(x)(pphi(x, epsilon((i+1)/2)));
    result =  integral(mulFunc(f, getShape_poly(phi, phi, 2, (i-1)/2*h, h)),(i-1)/2*h,(i-1)/2*h+h,'ArrayValued',true);
end