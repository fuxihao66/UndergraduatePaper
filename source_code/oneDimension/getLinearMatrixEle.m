function [result] = getLinearMatrixEle(i,j,h,p,q)

if abs(i-j)>=2
    result = 0;
elseif abs(i-j) == 1
    index = min(i,j);
    puv = integral(@(x)(p(x)*(-1/(h.^2))), index*h,(index+1)*h,'ArrayValued',true);
    quv = integral(@(x)(q(x)*( (x-index*h-h)/(-h) )*( (x-index*h)/h )), index*h,(index+1)*h,'ArrayValued',true);
    result = puv + quv;
else
    puv = integral(@(x)(p(x)*(1/(h.^2))), (i-1)*h,(i+1)*h,'ArrayValued',true);
    quv = integral(@(x)(q(x)*(1/(h.^2))*(x-(i-1)*h)*(x-(i-1)*h)), (i-1)*h,i*h,'ArrayValued',true);
    quv = quv + integral(@(x)(q(x)*(1/(h.^2))*(x-(i+1)*h)*(x-(i+1)*h)), i*h,(i+1)*h,'ArrayValued',true);
    result = puv + quv;
end