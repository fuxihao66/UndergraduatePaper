function [result] = LinearMatrixTwoDimen(i,j,h1,h2,x_num)
% numbering start from (1,1)
numNodes = x_num - 1;
i_x = mod(i-1, numNodes)+1;
i_y = floor((i-1)/numNodes)+1;
j_x = mod(j-1, numNodes)+1;
j_y = floor((j-1)/numNodes)+1;

h1_sq = h1.*h1;
h2_sq = h2.*h2;

if (i_x-j_x)^2 > 1 || (i_y-j_y)^2 > 1
    result = 0;
elseif (i_x-j_x)^2+(i_y-j_y)^2 == 0
    xx = i_x*h1;
    yy = -i_y*h2;
    func1 = @(x,y)( (1-(yy-y)./h2).^2./h1_sq + (1-(x-xx)./h1).^2./h2_sq ); 
    func2 = @(x,y)( (1-(yy-y)./h2).^2./h1_sq + (1-(xx-x)./h1).^2./h2_sq );
    func3 = @(x,y)( (1-(y-yy)./h2).^2./h1_sq + (1-(x-xx)./h1).^2./h2_sq ); 
    func4 = @(x,y)( (1-(y-yy)./h2).^2./h1_sq + (1-(xx-x)./h1).^2./h2_sq );
    
    result = integral2(func1,i_x*h1, i_x*h1+h1, -i_y*h2-h2, -i_y*h2);
    result = result + integral2(func2, i_x*h1-h1, i_x*h1, -i_y*h2-h2, -i_y*h2);
    result = result + integral2(func3, i_x*h1, i_x*h1+h1, -i_y*h2, -i_y*h2+h2);
    result = result + integral2(func4, i_x*h1-h1, i_x*h1, -i_y*h2, -i_y*h2+h2);
elseif i_x == j_x
    maxIy = max(i_y, j_y);
    func1 = @(x,y)( (1-(-maxIy*h2+h2-y)./h2).*(1-(y+maxIy*h2)./h2)./h1_sq - (1-(x-i_x*h1)./h1).^2./h2_sq );  % 1, 3
    func2 = @(x,y)( (1-(-maxIy*h2+h2-y)./h2).*(1-(y+maxIy*h2)./h2)./h1_sq - (1-(i_x*h1-x)./h1).^2./h2_sq );  % 2, 4
    result = integral2(func1, i_x*h1, i_x*h1+h1, -maxIy*h2, -maxIy*h2+h2);   % right
    result = result + integral2(func2, i_x*h1-h1, i_x*h1, -maxIy*h2, -maxIy*h2+h2); % left
elseif i_y == j_y
    minIx = min(i_x, j_x);
    func1 = @(x,y)( -(1-(-i_y*h2-y)./h2).^2./h1_sq + (1-(x-minIx*h1)./h1).*(1-(minIx*h1+h1-x)./h1)./h2_sq );
    func2 = @(x,y)( -(1-(y+i_y*h2)./h2).^2./h1_sq + (1-(x-minIx*h1)./h1).*(1-(minIx*h1+h1-x)./h1)./h2_sq );
    result = integral2(func1, minIx*h1, minIx*h1+h1, -i_y*h2-h2, -i_y*h2);   % bottom
    result = result + integral2(func2, minIx*h1, minIx*h1+h1, -i_y*h2, -i_y*h2+h2); % top
%elseif (i_x-j_x)*(i_y-j_y) < 0   % reverse diag
%    maxIy = max(i_y, j_y);
%    minIx = min(i_x, j_x);
%    f = @(x,y)( -().*()./h1_sq - ().*()./h2_sq );
%    result = integral2(f, minIx*h1, minIx*h1+h1, -maxIy*h2, -maxIy*h2+h2);
else                             % diag
    maxIy = max(i_y, j_y);
    minIx = min(i_x, j_x);
    f = @(x,y)( -( 1-(-maxIy*h2+h2-y)./h2 ).*( 1-(y+maxIy*h2)./h2 )./h1_sq - ( 1-(x-minIx*h1)./h1 ).*( 1-(minIx*h1+h1-x)./h1 )./h2_sq );
    result = integral2(f, minIx*h1, minIx*h1+h1, -maxIy*h2, -maxIy*h2+h2);    
end