
% main for linear

methodCell = {'method'};
errorCell = {'error'};
hCell = {'h'};

scale = 2.5;
for sc = 0:5
    

scale = scale*2;  % num of unit
left_b = 0;   % left bound
right_b = 1;  % right bound


matrix = zeros(scale, scale);   % 
b = zeros(scale, 1);


p = @(x)(1);
q = @(x)(pi*pi/4);
f = @(x)(sin(pi*x/2)*pi*pi/2);

% calculate matrix
for i = 1:scale   
    for j = 1:scale
        matrix(i,j) = MatrixElementCal(i,j,left_b,(right_b-left_b)/scale, p, q, scale);
    end
end
% calculate b
for i = 1:scale
    b(i) = VectorElementCal(i, left_b,(right_b-left_b)/scale, f, scale);
end

%cond(matrix)

U = matrix\b;   % solve equations
U = [0,U']';

x_axis=left_b:(right_b-left_b)/scale:right_b;
%L1Error(sin(pi/2*x_axis),U)

e = ErrorCal(x_axis, U, @(x)(sin(pi/2.*x)) )

methodCell{end+1} = 'Linear';
hCell{end+1} = 1/scale;
errorCell{end+1} = e;

end


data = table(methodCell', hCell', errorCell');
writetable(data, 'linear.csv');
