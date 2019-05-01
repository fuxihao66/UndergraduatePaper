function [result] = getLinearVectorEle(i,h,f)


result = integral(@(x)( f(x)*(x-(i-1)*h)/h ),(i-1)*h,i*h,'ArrayValued',true);
result = result + integral(@(x)( -f(x)*(x-(i+1)*h)/h ),i*h,(i+1)*h,'ArrayValued',true);