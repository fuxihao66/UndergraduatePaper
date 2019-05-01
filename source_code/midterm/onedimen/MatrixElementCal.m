function [result] = MatrixElementCal(i, j, x_0, h, p, q, scale)
result = 0.0;
 
if i == j && i == scale
    result = result + integral(@(x)(p(x_0+(scale-1).*h+h.*x)),0,1,'ArrayValued',true)/h + integral(@(x)(q(x_0+(scale-1).*h+h.*x).*(x.^2)),0,1,'ArrayValued',true).*h;
elseif i == j
	result = result + integral(@(x)(p(x_0+(i-1).*h+h.*x)),0,1, 'ArrayValued',true)/h + integral(@(x)(q(x_0+(i-1).*h+h.*x).*(x.^2)),0,1, 'ArrayValued',true).*h;
    result = result + integral(@(x)(p(x_0+i.*h+h.*x)),0,1, 'ArrayValued',true)/h + integral(@(x)(q(x_0+i.*h+h.*x).*((1-x).^2)),0,1, 'ArrayValued',true).*h;
elseif abs(i-j) == 1
    minIndex = min(i,j);
	result = result + -1.*integral(@(x)(p(x_0+minIndex.*h+h.*x)), 0, 1,'ArrayValued',true)/h + integral(@(x)(q(x_0+minIndex.*h+h.*x).*x.*(1-x)), 0,1,'ArrayValued',true).*h;
end   