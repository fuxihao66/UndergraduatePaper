

% u(0) = 0, u(1) = 0
% -d(pdu/dx)/dx + qu = f
% virtual work: int_0^1 pu'v'dx  + int_0^1 qudx = int_0^1 fvdx

x_start = 0;
x_end = 1;
h = 0.05;
num_of_unit = (x_end-x_start)/h;


p = @(x)(1);
q = @(x)(0);
f = @(x)(sin(x));
% unit base
%phi1 = @(x)(a*phi(abs(x-x_sample(i-1))) + c*phi(abs(x-x_sample(i))) );
%phi2 = @(x)(b*phi(abs(x-x_sample(i-1))) + d*phi(abs(x-x_sample(i))) );



fem_mat = zeros(num_of_unit-1,num_of_unit-1);
fem_vec = zeros(num_of_unit-1,1);
for i=1:num_of_unit-1
    for j = 1:num_of_unit-1
        %fem_mat(i,j) = getGaussMatrixEle(i,j,h,p,q);
        fem_mat(i,j) = getLinearMatrixEle(i,j,h,p,q);
    end
end

for i = 1:num_of_unit-1
    %fem_vec(i) = getGaussVector(i,h,f);
    fem_vec(i) = getLinearVectorEle(i,h,f);
end

U = fem_mat\fem_vec;
U
U = [0,U',0]';

accSolution = zeros(num_of_unit+1,1);
for i = 1:num_of_unit+1
   accSolution(i,1) = sin((i-1)*h)-sin(1)*((i-1)*h);
end


L1Error(accSolution, U)

x_axis=x_start:h:x_end;

figure;
plot(x_axis,U, x_axis, accSolution);
legend('Numerical Solution', 'Accurate Solution');













