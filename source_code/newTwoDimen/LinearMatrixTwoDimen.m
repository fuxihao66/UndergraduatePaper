%function [result] = LinearMatrixTwoDimen(i,j,h1,h2,x_num)
function [result] = LinearMatrixTwoDimen(x_start, y_start, i_x, i_y, j_x, j_y, h1,h2,x_num)
% numbering start from (1,1)


h1_sq = h1.*h1;
h2_sq = h2.*h2;

x_real = x_start+i_x*h1;
y_real = y_start+i_y*h2;
if (i_x-j_x)^2 > 1 || (i_y-j_y)^2 > 1
    result = 0;
elseif (i_x-j_x)^2+(i_y-j_y)^2 == 0
    
    func1 = @(x,y)( (1-(y_real-y)./h2).^2./h1_sq + (1-(x-x_real)./h1).^2./h2_sq ); 
    func2 = @(x,y)( (1-(y_real-y)./h2).^2./h1_sq + (1-(x_real-x)./h1).^2./h2_sq );
    func3 = @(x,y)( (1-(y-y_real)./h2).^2./h1_sq + (1-(x-x_real)./h1).^2./h2_sq ); 
    func4 = @(x,y)( (1-(y-y_real)./h2).^2./h1_sq + (1-(x_real-x)./h1).^2./h2_sq );
    
    result = integral2(func1,x_real, x_real+h1, y_real-h2, y_real);
    result = result + integral2(func2, x_real-h1, x_real, y_real-h2, y_real);
    result = result + integral2(func3, x_real, x_real+h1, y_real, y_real+h2);
    result = result + integral2(func4, x_real-h1, x_real, y_real, y_real+h2);
elseif i_x == j_x
    minIy = min(i_y, j_y);
    minIy_real = y_start+minIy*h2;
    func1 = @(x,y)( (1-(minIy_real+h2-y)./h2).*(1-(y-minIy_real)./h2)./h1_sq - (1-(x-x_real)./h1).^2./h2_sq );  % 1, 3
    func2 = @(x,y)( (1-(minIy_real+h2-y)./h2).*(1-(y-minIy_real)./h2)./h1_sq - (1-(x_real-x)./h1).^2./h2_sq );  % 2, 4
    result = integral2(func1, x_real, x_real+h1, minIy_real, minIy_real+h2);   % right
    result = result + integral2(func2, x_real-h1, x_real, minIy_real, minIy_real+h2); % left
elseif i_y == j_y
    minIx = min(i_x, j_x);
    minIx_real = x_start+minIx*h1;
    func1 = @(x,y)( -(1-(y_real-y)./h2).^2./h1_sq + (1-(x-minIx_real)./h1).*(1-(minIx_real+h1-x)./h1)./h2_sq );
    func2 = @(x,y)( -(1-(y-y_real)./h2).^2./h1_sq + (1-(x-minIx_real)./h1).*(1-(minIx_real+h1-x)./h1)./h2_sq );
    result = integral2(func1, minIx_real, minIx_real+h1, y_real-h2, y_real);   % bottom
    result = result + integral2(func2, minIx_real, minIx_real+h1, y_real, y_real+h2); % top
%elseif (i_x-j_x)*(i_y-j_y) < 0   % reverse diag
%    maxIy = max(i_y, j_y);
%    minIx = min(i_x, j_x);
%    f = @(x,y)( -().*()./h1_sq - ().*()./h2_sq );
%    result = integral2(f, minIx*h1, minIx*h1+h1, -maxIy*h2, -maxIy*h2+h2);
else                             % diag
    minIy = min(i_y, j_y);
    minIy_real = y_start+minIy*h2;
    minIx = min(i_x, j_x);
    minIx_real = x_start+minIx*h1;
    f = @(x,y)( -( 1-(minIy_real+h2-y)./h2 ).*( 1-(y-minIy_real)./h2 )./h1_sq - ( 1-(x-minIx_real)./h1 ).*( 1-(minIx_real+h1-x)./h1 )./h2_sq );
    result = integral2(f, minIx_real, minIx_real+h1, minIy_real, minIy_real+h2);    
end





