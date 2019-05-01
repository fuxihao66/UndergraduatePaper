%function [result] = InvMQMatrixTwoDimen(x_start, y_start, i,j,h1,h2,x_num)
function [result] = GaussMatrixTwoDimen_one(x_start, y_start, i_x, i_y, j_x, j_y, h1,h2,x_num)
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
func_dx = @(x,y)(exp(-c.*c.*(x.*x+y.*y) ).*(-2.*c.*c.*x));
func_dy = @(x,y)(exp(-c.*c.*(x.*x+y.*y) ).*(-2.*c.*c.*y));
phi = @(r)(exp(-c*c*r*r));
% relative distance

indicator_x = 1;
indicator_y = 2;


if (i_x-j_x)^2 > 1 || (i_y-j_y)^2 > 1
    result = 0;
elseif (i_x-j_x)^2+(i_y-j_y)^2 == 0
    func1 = addFunctionHandle(powerFunctionHandle(getShapeFunc_one_deri(func_dx, phi, indicator_x, 1, i_x_real, i_y_real, h1, h2),2),...
                              powerFunctionHandle(getShapeFunc_one_deri(func_dy, phi, indicator_y, 1, i_x_real, i_y_real, h1, h2),2)); 
    func2 = addFunctionHandle(powerFunctionHandle(getShapeFunc_one_deri(func_dx, phi, indicator_x, 2, i_x_real, i_y_real, h1, h2),2),...
                              powerFunctionHandle(getShapeFunc_one_deri(func_dy, phi, indicator_y, 2, i_x_real, i_y_real, h1, h2),2));
    func3 = addFunctionHandle(powerFunctionHandle(getShapeFunc_one_deri(func_dx, phi, indicator_x, 3, i_x_real, i_y_real, h1, h2),2),...
                              powerFunctionHandle(getShapeFunc_one_deri(func_dy, phi, indicator_y, 3, i_x_real, i_y_real, h1, h2),2));
    func4 = addFunctionHandle(powerFunctionHandle(getShapeFunc_one_deri(func_dx, phi, indicator_x, 4, i_x_real, i_y_real, h1, h2),2),...
                              powerFunctionHandle(getShapeFunc_one_deri(func_dy, phi, indicator_y, 4, i_x_real, i_y_real, h1, h2),2));
    
    
    result = integral2(func1,i_x_real, i_x_real+h1, i_y_real-h2, i_y_real);
    result = result + integral2(func2, i_x_real-h1, i_x_real, i_y_real-h2, i_y_real);
    result = result + integral2(func3, i_x_real, i_x_real+h1, i_y_real, i_y_real+h2);
    result = result + integral2(func4, i_x_real-h1, i_x_real, i_y_real, i_y_real+h2);
elseif i_x == j_x
    minIy = min(i_y, j_y);
    min_iy_real = y_start + minIy*h2;
    func1 = addFunctionHandle(mulFunctionHandle(getShapeFunc_one_deri(func_dx, phi, indicator_x, 1, i_x_real, min_iy_real+h2, h1, h2),getShapeFunc_one_deri(func_dx, phi, indicator_x, 3, i_x_real, min_iy_real, h1, h2)),...
                              mulFunctionHandle(getShapeFunc_one_deri(func_dy, phi, indicator_y, 1, i_x_real, min_iy_real+h2, h1, h2),getShapeFunc_one_deri(func_dy, phi, indicator_y, 3, i_x_real, min_iy_real, h1, h2)));
    func2 = addFunctionHandle(mulFunctionHandle(getShapeFunc_one_deri(func_dx, phi, indicator_x, 2, i_x_real, min_iy_real+h2, h1, h2),getShapeFunc_one_deri(func_dx, phi, indicator_x, 4, i_x_real, min_iy_real, h1, h2)),...
                              mulFunctionHandle(getShapeFunc_one_deri(func_dy, phi, indicator_y, 2, i_x_real, min_iy_real+h2, h1, h2),getShapeFunc_one_deri(func_dy, phi, indicator_y, 4, i_x_real, min_iy_real, h1, h2)));
    
    result = integral2(func1,i_x_real, i_x_real+h1, min_iy_real, min_iy_real+h2);   % right
    result = result + integral2(func2, i_x_real-h1, i_x_real, min_iy_real, min_iy_real+h2); % left
elseif i_y == j_y
    minIx = min(i_x, j_x);
    min_ix_real = x_start + minIx*h1;
    func1 = addFunctionHandle(mulFunctionHandle(getShapeFunc_one_deri(func_dx, phi, indicator_x, 1, min_ix_real, i_y_real, h1, h2),getShapeFunc_one_deri(func_dx, phi, indicator_x, 2, min_ix_real+h1, i_y_real, h1, h2)),...
                              mulFunctionHandle(getShapeFunc_one_deri(func_dy, phi, indicator_y, 1, min_ix_real, i_y_real, h1, h2),getShapeFunc_one_deri(func_dy, phi, indicator_y, 2, min_ix_real+h1, i_y_real, h1, h2)));
    func2 = addFunctionHandle(mulFunctionHandle(getShapeFunc_one_deri(func_dx, phi, indicator_x, 3, min_ix_real, i_y_real, h1, h2),getShapeFunc_one_deri(func_dx, phi, indicator_x, 4, min_ix_real+h1, i_y_real, h1, h2)),...
                              mulFunctionHandle(getShapeFunc_one_deri(func_dy, phi, indicator_y, 3, min_ix_real, i_y_real, h1, h2),getShapeFunc_one_deri(func_dy, phi, indicator_y, 4, min_ix_real+h1, i_y_real, h1, h2)));
    
    result = integral2(func1, min_ix_real, min_ix_real+h1, i_y_real-h2, i_y_real);   % bottom
    result = result + integral2(func2, min_ix_real, min_ix_real+h1, i_y_real, i_y_real+h2); % top
elseif (i_x-j_x)*(i_y-j_y) > 0   %  diag
    minIy = min(i_y, j_y);
    minIx = min(i_x, j_x);
    min_ix_real = x_start + minIx*h1;
    min_iy_real = y_start + minIy*h2;
    f = addFunctionHandle(mulFunctionHandle(getShapeFunc_one_deri(func_dx, phi, indicator_x, 2, min_ix_real+h1, min_iy_real+h2, h1, h2),getShapeFunc_one_deri(func_dx, phi, indicator_x, 3, min_ix_real, min_iy_real, h1, h2)),...
                          mulFunctionHandle(getShapeFunc_one_deri(func_dy, phi, indicator_y, 2, min_ix_real+h1, min_iy_real+h2, h1, h2),getShapeFunc_one_deri(func_dy, phi, indicator_y, 3, min_ix_real, min_iy_real, h1, h2)));
    
    result = integral2(f, min_ix_real, min_ix_real+h1, min_iy_real, min_iy_real+h2);
    
    
else                             %  reverse diag
    minIy = min(i_y, j_y);
    minIx = min(i_x, j_x);
    min_ix_real = x_start + minIx*h1;
    min_iy_real = y_start + minIy*h2;
    f = addFunctionHandle(mulFunctionHandle(getShapeFunc_one_deri(func_dx, phi, indicator_x, 1, min_ix_real, min_iy_real+h2, h1, h2),getShapeFunc_one_deri(func_dx, phi, indicator_x, 4, min_ix_real+h1, min_iy_real, h1, h2)),...
                          mulFunctionHandle(getShapeFunc_one_deri(func_dy, phi, indicator_y, 1, min_ix_real, min_iy_real+h2, h1, h2),getShapeFunc_one_deri(func_dy, phi, indicator_y, 4, min_ix_real+h1, min_iy_real, h1, h2)));
    
    result = integral2(f, min_ix_real, min_ix_real+h1, min_iy_real, min_iy_real+h2);  
   
end


