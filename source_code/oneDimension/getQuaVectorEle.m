function [result] = getQuaVectorEle(i,h,f)


result = integral(@(x)( f(x)*(b*phi(x-(i-1)*h)+d*phi(x-i*h)) ),(i-1)*h,i*h,'ArrayValued',true);
result = result + integral(@(x)( f(x)*(a*phi(x-i*h)+c*phi(x-(i+1)*h)) ),i*h,(i+1)*h,'ArrayValued',true);