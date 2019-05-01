% 
% -(p^2u/px^2 + p^2u/py^2) = f
% virtual work: int (pu/px pv/px)+(pu/py pv/py) = int( fv )dxdy
% a(sum phi_i u_i, phi_j) = -a(sum_(x,y)in \Omega phi_k u_k, phi_j)
% analytical: 

x_start = -1;
x_end = 1;
y_start = -1;
y_end = 1;  
h1 = 0.05;
h2 = 0.05;
N1 = (x_end-x_start)/h1;
N2 = (y_end-y_start)/h2;

%boundary = @(x,y)(sin(pi.*x./3).*sinh(pi.*(y+1)./3)./sinh(pi./3));    % boundary is defined by u
%f = @(x,y)(2.*pi.^2.*sin(pi.*x).*sin(pi.*y));
%boundary = @(x,y)( sin(pi.*x).*sin(pi.*y) );


%f = @(x,y)(-4);
%boundary = @(x,y)(x.^2+y.^2);


%f = @(x,y)( -( 0.75.*exp(-((9.*x-2).^2+(9.*y-2).^2)/4).*(81.*(9.*x-2).^2./4+81.*(9.*y-2).^2/4-81)...
%            +( 0.75.*exp(-(9.*x+1).^2/49-(9.*y+1)/10).*(18.*18./(49.^2).*(9.*x+1).^2+81/100-162/49) ) ...
%            +( 0.5.*exp(-((9.*x-7).^2+(9.*y-3).^2)/4).*( 81.*(9.*x-7).^2./4+81.*(9.*y-3).^2./4-81 ) ) ...
%            -( 0.2.*exp(-(9.*x-4).^2-(9.*y-7).^2).*( 18.^2.*(9.*x-4).^2+18.^2.*(9.*y-7).^2-18.*18 ) ) ));
%boundary = @(x,y)( 0.75.*exp(-((9.*x-2).^2+(9.*y-2).^2)/4)+0.75.*exp(-(9.*x+1).^2/49-(9.*y+1)/10)+...
%                   0.5.*exp(-((9.*x-7).^2+(9.*y-3).^2)/4)-0.2.*exp(-(9.*x-4).^2-(9.*y-7).^2) );

               
%f = @(x,y)( -( 7.5.*(1-x/2).^4.*(1-y/2).^6 + 7.5.*(1-x/2).^6.*(1-y/2).^4 ...
%                + 1000.*(1-x).^3.*(x.^3).*( 6.*y.*(1-y).^3 + 6.*(y.^3).*(1-y)-18.*y.^2.*(1-y).^2 ) ...
%                + 1000.*(1-y).^3.*(y.^3).*( 6.*x.*(1-x).^3 + 6.*(x.^3).*(1-x)-18.*x.^2.*(1-x).^2 ) ...
%                + 7.5.*y.^6.*(1-x/2).^4 + 7.5.*x.^6.*(1-y/2).^4 ...
%                + 30.*x.^4.*(1-y/2).^6 + 30.*y.^4.*(1-x/2).^6 ) );
%boundary = @(x,y)( (1-x/2).^6.*(1-y/2).^6 + 1000.*(1-x).^3.*(x.^3).*(1-y).^3.*(y.^3)+y.^6.*(1-x/2).^6+x.^6.*(1-y/2).^6 );

% inner
u1 = @(x,y)(1./((x+0.5).^2+(y-0.5).^2+0.01));
u2 = @(x,y)(1./((x-0.5).^2+(y+0.5).^2+0.01));
f = @(x,y)(4.*u1(x,y).^2-8.*u1(x,y).^3.*((x+0.5).^2+(y-0.5).^2)-4.*u2(x,y).^2+8.*u2(x,y).^3.*((x-0.5).^2+(y+0.5).^2) );
boundary = @(x,y)(u1(x,y)-u2(x,y));

% boundary singularity
% f = @(x,y)(2./(x.^2+0.01)-8.*x^2./((x.^2+0.01).^3)-2);
% boundary = @(x,y)(1./(x.^2+0.01)+y.^2);


% f = @(x,y)(1.6.*y.^2./(x.^2+y.^2).^1.8+1.6.*x.^2./(x.^2+y.^2).^1.8 - 0.8./(x.^2+y.^2).^0.8);
% boundary = @(x,y)((x.^2+y.^2).^0.2);

num_of_undetermined = (N1-1)*(N2-1);

fem_mat = zeros(num_of_undetermined, num_of_undetermined);
fem_vec = zeros(num_of_undetermined, 1);
for i=1:num_of_undetermined
    for j = 1:num_of_undetermined
        numNodes = N1 - 1;
        i_x = mod(i-1, numNodes)+1;
        i_y = floor((i-1)/numNodes)+1;
        j_x = mod(j-1, numNodes)+1;
        j_y = floor((j-1)/numNodes)+1;
        
        %fem_mat(i,j) = LinearMatrixTwoDimen(x_start, y_start, i_x, i_y, j_x, j_y,h1,h2,N1);
        
        fem_mat(i,j) = InvMQMatrixTwoDimen(x_start, y_start, i_x, i_y, j_x, j_y,h1,h2,N1);
    end
end

for i = 1:num_of_undetermined
    %fem_vec(i) = LinearVectorTwoDimen(x_start, y_start,i,h1,h2,f,N1,N2,boundary);
    
    fem_vec(i) = InvMQVectorTwoDimen(x_start, y_start,i,h1,h2,f,N1,N2,boundary);
end

U = fem_mat\fem_vec;


x = x_start:h1:x_end; % x range
y = y_start:h2:y_end;
[xx,yy]=meshgrid(x,y); %??????

% row is y and column is x
zz = zeros(N2+1, N1+1);
for i = 1:N1+1
    for j = 1:N2+1
        if j == N2+1 || i == 1 || i == N1+1 || j == 1
            zz(j,i) = boundary(x_start+(i-1)*h1, y_start+(j-1)*h2);
        else
            zz(j,i) = U((i-1)+(j-2)*(N1-1));
            %zz(j,i) = boundary(x_start+(i-1)*h1, y_start+(j-1)*h2);
        end
    end
end

c = 1.2;

RBFL2Error_nonpoly(x, y, zz, boundary, @(x,y)( 1./sqrt(1+c.*c.*(x.*x+y.*y)) ) )

%RBFL2Error_nonpoly(x, y, zz, boundary, @(x,y)(exp(-c.*c.*(x.^2+y.^2))) )
%LinearL2Error(x, y, zz, boundary)
surf(xx,yy,zz);

%surf(xx,yy,accSolution);
shading interp;