

scale = 32;  % num of unit
left_b = 0;   % left bound
right_b = 1;  % right bound
h = (right_b-left_b)/scale;

matrix = zeros(scale*2, scale*2);   % 
b = zeros(scale*2, 1);

% index = 1,3,5,.. means mid point
% index = 2,4,6,.. means node

p = @(x)(1);
q = @(x)(pi.*pi./4);
f = @(x)(sin(pi*x/2)*pi*pi/2);
solu = @(x)(sin(pi/2.*x));
f_d = @(x)(cos(pi.*x/2).*pi.*pi.*pi./4);

epsilon = zeros(scale);

for i=1:scale
    index = i*2-1;
    epsilon(i) = f_d(index.*h./2).*h/(3.*(solu((index+1).*h./2)-solu((index-1).*h./2))) - q(1)/3;
end

%f = @(x)(pi.^2.*(2.*x-x.^2)/4+2);
%solu = @(x)(2.*x-x.^2);

%f = @(x)( 6.*(1-x)+pi.^2./4.*(1-(1-x).^3) );
%solu = @(x)(1-(1-x).^3);

%f = @(x)( (4.*x.^2-8.*x+6).*exp((x-1).^2)+pi.^2./4.*(2.7182818-exp((x-1).^2)) );
%solu = @(x)( 2.7182818-exp((x-1).^2) );

% calculate matrix
for i = 1:scale*2  
    for j = 1:scale*2
        matrix(i,j) = getMQMatrixEle_poly(i,j,left_b,(right_b-left_b)/scale, p, q, scale, epsilon);
        %matrix(i,j) = getGaussMatrixEle(i,j,left_b,(right_b-left_b)/scale, p, q, scale);
        %matrix(i,j) = getGaussMatrixEle_poly(i,j,left_b,(right_b-left_b)/scale, p, q, scale);
    end
end
% calculate b
for i = 1:scale*2
    b(i) = getMQVector_poly(i, left_b,(right_b-left_b)/scale, f, scale, epsilon);
    %b(i) = getGaussVector(i, left_b,(right_b-left_b)/scale, f, scale);
    %b(i) = getGaussVector_poly(i, left_b,(right_b-left_b)/scale, f, scale);
end

U = matrix\b;


U = [0,U']';
x_axis=left_b:(right_b-left_b)/scale/2:right_b;

%c = 0.8;
%phi = @(x)(exp(-c.*x.*x));
%phi = @(x,c)(exp(-c.*x.*x));
pphi = @(x, c)(sqrt(1+c.*x.*x));
%RBF_ERROR(x_axis, U, solu, phi)

RBF_ERROR_varyC(x_axis, U, solu, pphi, epsilon)

figure;
plot(x_axis,U,x_axis,sin(pi/2*x_axis));
legend('Numerical Solution','Accurate Solution');





