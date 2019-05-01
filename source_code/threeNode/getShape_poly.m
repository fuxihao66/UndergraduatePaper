function [result] = getShape_poly(phi, func, index, startX, h)


f1 = @(x)func(x-startX);
f2 = @(x)func(x-(startX+h/2));
f3 = @(x)func(x- (startX+h));
mat = [phi(0),phi(h/2),phi(h),1;phi(h/2),phi(0),phi(h/2),1;phi(h), phi(h/2), phi(0),1;1,1,1,0];
mat_inv = inv(mat);

result = @(x)(mat_inv(1,index).*f1(x)+mat_inv(2,index).*f2(x)+mat_inv(3,index).*f3(x)+mat_inv(4,index));