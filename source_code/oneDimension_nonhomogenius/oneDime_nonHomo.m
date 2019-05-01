% u(0) = alpha, u'(1) = beta
% -d(pdu/dx)/dx + qu = f
% virtual work: int_0^1 pu'v'dx  + int_0^1 quvdx = int_0^1 fvdx + beta*v(1)*p(1)
% ??rbf  Hermite

% <------------------------------- one more node
x_start = 0;
x_end = 1;
h = 0.2;
num_of_unit = (x_end-x_start)/h;


p = @(x)(1);
q = @(x)(0);
f = @(x)(sin(x));
% unit base
%phi1 = @(x)(a*phi(abs(x-x_sample(i-1))) + c*phi(abs(x-x_sample(i))) );
%phi2 = @(x)(b*phi(abs(x-x_sample(i-1))) + d*phi(abs(x-x_sample(i))) );



fem_mat = zeros(num_of_unit,num_of_unit);
fem_vec = zeros(num_of_unit,1);
for i=1:num_of_unit-1
    for j = 1:num_of_unit-1
        fem_mat(i,j) = getMatrixEle(i,j,h,p,q);
    end
end

for i = 1:num_of_unit-1
    fem_vec(i) = getVector(i,h,f);
end

U = fem_mat\fem_vec;
U
U = [U',0]';

accSolution = zeros(num_of_unit+1,1);
for i = 1:num_of_unit+1
   accSolution(i,1) = sin((i-1)*h)-sin(1)*((i-1)*h);
end




x_axis=x_start:h:x_end;

figure;
plot(x_axis,U, x_axis, accSolution);
legend('Numerical Solution', 'Accurate Solution');

