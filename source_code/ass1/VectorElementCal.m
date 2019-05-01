function [result] = VectorElementCal(j, x_0, h, f, scale)

if j == scale
    result = integral(@(x)(f(x_0+(j-1).*h+h.*x).*x),0,1,'ArrayValued',true).*h;
else
    result = integral(@(x)(f(x_0+(j-1).*h+h.*x).*x),0,1,'ArrayValued',true).*h;
    result = result + integral(@(x)(f(x_0+j.*h+h.*x).*(1-x)),0,1,'ArrayValued',true).*h;
end