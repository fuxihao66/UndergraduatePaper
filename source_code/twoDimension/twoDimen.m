% y = 0 sin(pi/3)x
% p^2u/px^2 + p^2u/py^2 = 0
% virtual work: int (pu/px pv/px)+(pu/py pv/py)
% a(sum phi_i u_i, phi_j) = -a(sum_y=0 phi_k u_k, phi_j)
% analytical: sin(pix/3)sinh(pi(y+1)/3)/sinh(pi/3)

x_start = 0;
x_end = 3;
y_start = -1;
y_end = 0;  
h1 = 0.5;
h2 = 0.25;
N1 = (x_end-x_start)/h1;
N2 = (y_end-y_start)/h2;

boundary = zeros(N1-1,1);
for i=1:N1-1
    boundary(i) = sin(pi*i*h1/3);
end

num_of_undetermined = (N1-1)*(N2-1);

fem_mat = zeros(num_of_undetermined, num_of_undetermined);
fem_vec = zeros(num_of_undetermined, 1);
for i=1:num_of_undetermined
    for j = 1:num_of_undetermined
        %fem_mat(i,j) = InvMQMatrixTwoDimen(i,j,h1,h2,N1);
        %fem_mat(i,j) = GaussMatrixTwoDimen(i,j,h1,h2,N1);
        fem_mat(i,j) = LinearMatrixTwoDimen(i,j,h1,h2,N1);
    end
end

for i = 1:num_of_undetermined
    %fem_vec(i) = InvMQVectorTwoDimen(i,h1,h2,N1,boundary);
    %fem_vec(i) = GaussVectorTwoDimen(i,h1,h2,N1,boundary);
    fem_vec(i) = LinearVectorTwoDimen(i,h1,h2,N1,boundary);
end

U = fem_mat\fem_vec;


accSolution = zeros(N2+1,N1+1);
for i = 1:N1+1
    for j = 1:N2+1
        accSolution(j,i) = sin(pi*(i-1)*h1/3)*sinh(pi*(1-(j-1)*h2)/3)/sinh(pi/3);
    end
end


x=x_start:h1:x_end; % x range
%y=y_start:h2:y_end; % y range
y = y_end:-h2:y_start;
[xx,yy]=meshgrid(x,y); %??????

% row is y and column is x
zz = zeros(N2+1, N1+1);
for i = 1:N1+1
    for j = 1:N2+1
        if j == N2+1 || i == 1 || i == N1+1
            zz(j,i) = 0;
        elseif j == 1
            zz(j,i) = boundary(i-1);
        else
            zz(j,i) = U((i-1)+(j-2)*(N1-1));
        end
    end
end

c = 0.81;

%RBFL2Error_nonpoly(x, y, zz, @(x,y)(sin(pi.*x./3).*sinh(pi.*(y+1)./3)./sinh(pi./3)), @(x,y)( 1./sqrt(1+c.*c.*(x.*x+y.*y)) ) )
%RBFL2Error_nonpoly(x, y, zz, @(x,y)(sin(pi.*x./3).*sinh(pi.*(y+1)./3)./sinh(pi./3)), @(x,y)(exp(-c.*c.*(x.^2+y.^2))) )
LinearL2Error(x, y, zz, @(x,y)(sin(pi.*x./3).*sinh(pi.*(y+1)./3)./sinh(pi./3)))
surf(xx,yy,zz);

%surf(xx,yy,accSolution);
shading interp;


%x_1 = reshape(xx,[1,(N2+1)*(N1+1)])
%y_1 = reshape(yy,[1,(N2+1)*(N1+1)])
%z_1 = reshape(zz,[1,(N2+1)*(N1+1)])
%scatter3(x_1,y_1,z_1);
