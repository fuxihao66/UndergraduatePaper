

scale = 50;  % num of unit
left_b = 0;   % left bound
right_b = 1;  % right bound


matrix = zeros(scale, scale);   % 
b = zeros(scale, 1);


p = @(x)(1);
q = @(x)(pi*pi/4);
f = @(x)(sin(pi*x/2)*pi*pi/2);
%solu = @(x)(sin(pi/2.*x))

%f = @(x)(pi.^2.*(2.*x-x.^2)/4+2);
%solu = @(x)(2.*x-x.^2);

%f = @(x)( 6.*(1-x)+pi.^2./4.*(1-(1-x).^3) );
%solu = @(x)(1-(1-x).^3);

%f = @(x)( (4.*x.^2-8.*x+6).*exp((x-1).^2)+pi.^2./4.*(2.7182818-exp((x-1).^2)) );
%solu = @(x)( 2.7182818-exp((x-1).^2) );

% calculate matrix
for i = 1:scale   
    for j = 1:scale
        %matrix(i,j) = getMQMatrixEle(i,j,left_b,(right_b-left_b)/scale, p, q, scale);
        %matrix(i,j) = getMQMatrixElePoly(i,j,left_b,(right_b-left_b)/scale, p, q, scale);
        %matrix(i,j) = getMQMatrixEleAnother(i,j,left_b,(right_b-left_b)/scale, p, q, scale);
        %matrix(i,j) = getMQMatrixEleAnotherPoly(i,j,left_b,(right_b-left_b)/scale, p, q, scale);
        %matrix(i,j) = getInvMQMatrixEle(i,j,left_b,(right_b-left_b)/scale, p, q, scale);
        matrix(i,j) = getGaussMatrixEle(i,j,left_b,(right_b-left_b)/scale, p, q, scale);
        %matrix(i,j) = getGaussMatrixElePoly(i,j,left_b,(right_b-left_b)/scale, p, q, scale);
        %matrix(i,j) = MatrixElementCal(i,j,left_b,(right_b-left_b)/scale, p, q, scale);
        %matrix(i,j) = getGaussMatrixElePolyOneOrder(i,j,left_b,(right_b-left_b)/scale, p, q, scale);
    end
end
% calculate b
for i = 1:scale
    %b(i) = getMQVector(i, left_b,(right_b-left_b)/scale, f, scale);
    %b(i) = getMQVectorPoly(i, left_b,(right_b-left_b)/scale, f, scale);
    %b(i) = getMQVectorAnother(i, left_b,(right_b-left_b)/scale, f, scale);
    %b(i) = getMQVectorAnotherPoly(i, left_b,(right_b-left_b)/scale, f, scale);
    %b(i) = getInvMQVector(i, left_b,(right_b-left_b)/scale, f, scale);
    b(i) = getGaussVector(i, left_b,(right_b-left_b)/scale, f, scale);
    %b(i) = getGaussVectorPoly(i, left_b,(right_b-left_b)/scale, f, scale);
    %b(i) = VectorElementCal(i, left_b,(right_b-left_b)/scale, f, scale);
    %b(i) = getGaussVectorPolyOneOrder(i, left_b,(right_b-left_b)/scale, f, scale);
end



U = matrix\b;   % solve equations
U = [0,U']';




x_axis=left_b:(right_b-left_b)/scale:right_b;
%L1Error(sin(pi/2*x_axis),U)
%RBFErrorL2_nonpoly(x_axis, U, @(x)(sin(pi/2.*x)),@(x)(sqrt(c.^2+x.*x)) )
%RBFErrorL2(x_axis, U, @(x)(sin(pi/2.*x)), @(x)(sqrt(c.^2+x.*x)) )

c = 0.6;

%RBFErrorL2_OneOrder(x_axis, U, @(x)(sin(pi/2.*x)), @(x)(exp(-c.*x.*x)) )

%RBFErrorL2_nonpoly(x_axis, U, @(x)(sin(pi/2.*x)), @(x)(1./sqrt(1+c.*x.*x)) )

%RBFErrorL2_nonpoly(x_axis, U, @(x)(sin(pi/2.*x)), @(x)(sqrt(1+c.*x.*x)))

%RBFErrorL2(x_axis, U, @(x)(sin(pi/2.*x)), @(x)(exp(-c.*x.*x)) )
RBFErrorL2_nonpoly(x_axis, U, @(x)(sin(pi/2.*x)), @(x)(exp(-c.*x.*x)) )
%ErrorCal(x_axis, U, @(x)(sin(pi/2.*x)) )    % for linear

figure;
plot(x_axis,U,x_axis,sin(pi/2*x_axis));
legend('Numerical Solution','Accurate Solution');




