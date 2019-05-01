%function [result] = InvMQMatrixTwoDimen(x_start, y_start, i,j,h1,h2,x_num)
function [result] = IQMatrixTwoDimen(x_start, y_start, i_x, i_y, j_x, j_y, h1,h2,x_num)
% numbering start from (1,1)
%numNodes = x_num - 1;
%i_x = mod(i-1, numNodes)+1;
%i_y = floor((i-1)/numNodes)+1;
%j_x = mod(j-1, numNodes)+1;
%j_y = floor((j-1)/numNodes)+1;

i_x_real = x_start + i_x*h1;
i_y_real = y_start + i_y*h2;
j_x_real = x_start + j_x*h1;
j_y_real = y_start + j_y*h2;

c = 0.81;
func_dx = @(x,y)( -2.*c.*c.*x./( (1+c.*c.*(x.*x+y.*y)).^2 ) );
func_dy = @(x,y)( -2.*c.*c.*y./( (1+c.*c.*(x.*x+y.*y)).^2 ) );
phi = @(r)( 1./(1+c.*c.*r.*r) );
% relative distance

relaPoint = [0, 0; h1,0; 0, -h2; h1,-h2];

shapeMat = zeros(4,4);
for index_1 = 1:4
    for index_2 = 1:4
        shapeMat(index_1,index_2) = phi(norm( relaPoint([index_1],[1,2])-relaPoint([index_2],[1,2]) ));
    end
end
shapeMat_inv = inv(shapeMat);



if (i_x-j_x)^2 > 1 || (i_y-j_y)^2 > 1
    result = 0;
elseif (i_x-j_x)^2+(i_y-j_y)^2 == 0
    func1 = addFunctionHandle(powerFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 1, i_x_real, i_y_real, h1, h2),2),...
                              powerFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 1, i_x_real, i_y_real, h1, h2),2)); 
    func2 = addFunctionHandle(powerFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 2, i_x_real, i_y_real, h1, h2),2),...
                              powerFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 2, i_x_real, i_y_real, h1, h2),2));
    func3 = addFunctionHandle(powerFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 3, i_x_real, i_y_real, h1, h2),2),...
                              powerFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 3, i_x_real, i_y_real, h1, h2),2));
    func4 = addFunctionHandle(powerFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 4, i_x_real, i_y_real, h1, h2),2),...
                              powerFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 4, i_x_real, i_y_real, h1, h2),2));
    
    
    result = integral2(func1,i_x_real, i_x_real+h1, i_y_real-h2, i_y_real);
    result = result + integral2(func2, i_x_real-h1, i_x_real, i_y_real-h2, i_y_real);
    result = result + integral2(func3, i_x_real, i_x_real+h1, i_y_real, i_y_real+h2);
    result = result + integral2(func4, i_x_real-h1, i_x_real, i_y_real, i_y_real+h2);
elseif i_x == j_x
    minIy = min(i_y, j_y);
    min_iy_real = y_start + minIy*h2;
    func1 = addFunctionHandle(mulFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 1, i_x_real, min_iy_real+h2, h1, h2),getShapeFunc(func_dx, shapeMat_inv, 3, i_x_real, min_iy_real, h1, h2)),...
                              mulFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 1, i_x_real, min_iy_real+h2, h1, h2),getShapeFunc(func_dy, shapeMat_inv, 3, i_x_real, min_iy_real, h1, h2)));
    func2 = addFunctionHandle(mulFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 2, i_x_real, min_iy_real+h2, h1, h2),getShapeFunc(func_dx, shapeMat_inv, 4, i_x_real, min_iy_real, h1, h2)),...
                              mulFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 2, i_x_real, min_iy_real+h2, h1, h2),getShapeFunc(func_dy, shapeMat_inv, 4, i_x_real, min_iy_real, h1, h2)));
    
    result = integral2(func1,i_x_real, i_x_real+h1, min_iy_real, min_iy_real+h2);   % right
    result = result + integral2(func2, i_x_real-h1, i_x_real, min_iy_real, min_iy_real+h2); % left
elseif i_y == j_y
    minIx = min(i_x, j_x);
    min_ix_real = x_start + minIx*h1;
    func1 = addFunctionHandle(mulFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 1, min_ix_real, i_y_real, h1, h2),getShapeFunc(func_dx, shapeMat_inv, 2, min_ix_real+h1, i_y_real, h1, h2)),...
                              mulFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 1, min_ix_real, i_y_real, h1, h2),getShapeFunc(func_dy, shapeMat_inv, 2, min_ix_real+h1, i_y_real, h1, h2)));
    func2 = addFunctionHandle(mulFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 3, min_ix_real, i_y_real, h1, h2),getShapeFunc(func_dx, shapeMat_inv, 4, min_ix_real+h1, i_y_real, h1, h2)),...
                              mulFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 3, min_ix_real, i_y_real, h1, h2),getShapeFunc(func_dy, shapeMat_inv, 4, min_ix_real+h1, i_y_real, h1, h2)));
    
    result = integral2(func1, min_ix_real, min_ix_real+h1, i_y_real-h2, i_y_real);   % bottom
    result = result + integral2(func2, min_ix_real, min_ix_real+h1, i_y_real, i_y_real+h2); % top
elseif (i_x-j_x)*(i_y-j_y) > 0   %  diag
    minIy = min(i_y, j_y);
    minIx = min(i_x, j_x);
    min_ix_real = x_start + minIx*h1;
    min_iy_real = y_start + minIy*h2;
    f = addFunctionHandle(mulFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 2, min_ix_real+h1, min_iy_real+h2, h1, h2),getShapeFunc(func_dx, shapeMat_inv, 3, min_ix_real, min_iy_real, h1, h2)),...
                          mulFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 2, min_ix_real+h1, min_iy_real+h2, h1, h2),getShapeFunc(func_dy, shapeMat_inv, 3, min_ix_real, min_iy_real, h1, h2)));
    
    result = integral2(f, min_ix_real, min_ix_real+h1, min_iy_real, min_iy_real+h2);
    
    
else                             %  reverse diag
    minIy = min(i_y, j_y);
    minIx = min(i_x, j_x);
    min_ix_real = x_start + minIx*h1;
    min_iy_real = y_start + minIy*h2;
    f = addFunctionHandle(mulFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 1, min_ix_real, min_iy_real+h2, h1, h2),getShapeFunc(func_dx, shapeMat_inv, 4, min_ix_real+h1, min_iy_real, h1, h2)),...
                          mulFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 1, min_ix_real, min_iy_real+h2, h1, h2),getShapeFunc(func_dy, shapeMat_inv, 4, min_ix_real+h1, min_iy_real, h1, h2)));
    
    result = integral2(f, min_ix_real, min_ix_real+h1, min_iy_real, min_iy_real+h2);  
   
end



