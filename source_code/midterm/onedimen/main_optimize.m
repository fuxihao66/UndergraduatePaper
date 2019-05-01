c = 0.6;

for mm = 1:10

scale = 20;  % num of unit
left_b = 0;   % left bound
right_b = 1;  % right bound

x_axis=left_b:(right_b-left_b)/scale:right_b;

matrix = zeros(scale, scale);   % 
b = zeros(scale, 1);


p = @(x)(1);
q = @(x)(pi*pi/4);
f = @(x)(sin(pi*x/2)*pi*pi/2);


% <--------- Gauss -------------->

for i = 1:scale   
    for j = 1:scale
        matrix(i,j) = getGaussMatrixEle(i,j,left_b,(right_b-left_b)/scale, p, q, scale, c);
    end
end
% calculate b
for i = 1:scale
    b(i) = getGaussVector(i, left_b,(right_b-left_b)/scale, f, scale, c);
end

U = matrix\b;   % solve equations
U = [0,U']';

%L1Error(sin(pi/2*x_axis),U)

e = RBFErrorL2_nonpoly(x_axis, U, @(x)(sin(pi/2.*x)), @(x)(exp(-c.*x.*x)) );
c = getOptimizePara(x_axis, U, @(x)(sin(pi/2.*x)), @(a, x)(exp(-a.*x.*x)), c)


end
