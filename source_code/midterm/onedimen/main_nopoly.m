
% main for non_poly


methodCell = {'method'};
errorCell = {'error'};
hCell = {'h'};
cCell = {'c'};



for c = 0.5:0.1:3.0
    scale = 2.5;
for sc = 1:5
    

scale = scale*2  % num of unit
left_b = 0;   % left bound
right_b = 1;  % right bound

x_axis=left_b:(right_b-left_b)/scale:right_b;

matrix = zeros(scale, scale);   % 
b = zeros(scale, 1);


p = @(x)(1);
q = @(x)(pi*pi/4);
f = @(x)(sin(pi*x/2)*pi*pi/2);


% <--------- MQ -------------->

% calculate matrix
for i = 1:scale   
    for j = 1:scale
        matrix(i,j) = getMQMatrixEleAnother(i,j,left_b,(right_b-left_b)/scale, p, q, scale, c);
    end
end
% calculate b
for i = 1:scale
    b(i) = getMQVectorAnother(i, left_b,(right_b-left_b)/scale, f, scale, c);
end

U = matrix\b;   % solve equations
U = [0,U']';

%L1Error(sin(pi/2*x_axis),U)

e = RBFErrorL2_nonpoly(x_axis, U, @(x)(sin(pi/2.*x)), @(x)(sqrt(1+c.*x.*x)) )
methodCell{end+1} = 'MQ';
hCell{end+1} = 1/scale;
cCell{end+1} = c;
errorCell{end+1} = e;







% <--------- InvMQ -------------->

for i = 1:scale   
    for j = 1:scale
        matrix(i,j) = getInvMQMatrixEle(i,j,left_b,(right_b-left_b)/scale, p, q, scale, c);
    end
end
% calculate b
for i = 1:scale
    b(i) = getInvMQVector(i, left_b,(right_b-left_b)/scale, f, scale, c);
end

U = matrix\b;   % solve equations
U = [0,U']';

%L1Error(sin(pi/2*x_axis),U)

e = RBFErrorL2_nonpoly(x_axis, U, @(x)(sin(pi/2.*x)), @(x)(1./sqrt(1+c.*x.*x)) )

methodCell{end+1} = 'InvMQ';
hCell{end+1} = 1/scale;
cCell{end+1} = c;
errorCell{end+1} = e;






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

e = RBFErrorL2_nonpoly(x_axis, U, @(x)(sin(pi/2.*x)), @(x)(exp(-c.*x.*x)) )
methodCell{end+1} = 'Gauss';
hCell{end+1} = 1/scale;
cCell{end+1} = c;
errorCell{end+1} = e;

end
end

data = table(methodCell', hCell', errorCell', cCell');
writetable(data, 'nonPoly.csv');
