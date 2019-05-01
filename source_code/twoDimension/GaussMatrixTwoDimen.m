function [result] = GaussMatrixTwoDimen(i,j,h1,h2,x_num)
% numbering start from (1,1)
numNodes = x_num - 1;
i_x = mod(i-1, numNodes)+1;
i_y = floor((i-1)/numNodes)+1;
j_x = mod(j-1, numNodes)+1;
j_y = floor((j-1)/numNodes)+1;

c = 0.611;
func_dx = @(x,y)(exp(-c.*c.*(x.*x+y.*y) ).*(-2.*c.*c.*x));
func_dy = @(x,y)(exp(-c.*c.*(x.*x+y.*y) ).*(-2.*c.*c.*y));
phi = @(r)(exp(-c*c*r*r));
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
    func1 = addFunctionHandle(powerFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 1, i_x*h1, -i_y*h2, h1, h2),2),...
                              powerFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 1, i_x*h1, -i_y*h2, h1, h2),2)); 
    func2 = addFunctionHandle(powerFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 2, i_x*h1, -i_y*h2, h1, h2),2),...
                              powerFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 2, i_x*h1, -i_y*h2, h1, h2),2));
    func3 = addFunctionHandle(powerFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 3, i_x*h1, -i_y*h2, h1, h2),2),...
                              powerFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 3, i_x*h1, -i_y*h2, h1, h2),2));
    func4 = addFunctionHandle(powerFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 4, i_x*h1, -i_y*h2, h1, h2),2),...
                              powerFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 4, i_x*h1, -i_y*h2, h1, h2),2));
    
    %func2 = @(x,y)( getShapeFunc(func_dx, shapeMat_inv, 2, i_x*h1, -i_y*h2, h1, h2).^2 + getShapeFunc(func_dy, shapeMat_inv, 2, i_x*h1, -i_y*h2, h1, h2).^2 );
    %func3 = @(x,y)( getShapeFunc(func_dx, shapeMat_inv, 3, i_x*h1, -i_y*h2, h1, h2).^2 + getShapeFunc(func_dy, shapeMat_inv, 3, i_x*h1, -i_y*h2, h1, h2).^2 );
    %func4 = @(x,y)( getShapeFunc(func_dx, shapeMat_inv, 4, i_x*h1, -i_y*h2, h1, h2).^2 + getShapeFunc(func_dy, shapeMat_inv, 4, i_x*h1, -i_y*h2, h1, h2).^2 );
    result = integral2(func1,i_x*h1, i_x*h1+h1, -i_y*h2-h2, -i_y*h2);
    result = result + integral2(func2, i_x*h1-h1, i_x*h1, -i_y*h2-h2, -i_y*h2);
    result = result + integral2(func3, i_x*h1, i_x*h1+h1, -i_y*h2, -i_y*h2+h2);
    result = result + integral2(func4, i_x*h1-h1, i_x*h1, -i_y*h2, -i_y*h2+h2);
elseif i_x == j_x
    maxIy = max(i_y, j_y);
    func1 = addFunctionHandle(mulFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 1, i_x*h1, -maxIy*h2+h2, h1, h2),getShapeFunc(func_dx, shapeMat_inv, 3, i_x*h1, -maxIy*h2, h1, h2)),...
                              mulFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 1, i_x*h1, -maxIy*h2+h2, h1, h2),getShapeFunc(func_dy, shapeMat_inv, 3, i_x*h1, -maxIy*h2, h1, h2)));
    func2 = addFunctionHandle(mulFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 2, i_x*h1, -maxIy*h2+h2, h1, h2),getShapeFunc(func_dx, shapeMat_inv, 4, i_x*h1, -maxIy*h2, h1, h2)),...
                              mulFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 2, i_x*h1, -maxIy*h2+h2, h1, h2),getShapeFunc(func_dy, shapeMat_inv, 4, i_x*h1, -maxIy*h2, h1, h2)));
    %func1 = @(x,y)(getShapeFunc(func_dx, shapeMat_inv, 1, i_x*h1, -maxIy*h2+h2, h1, h2)*getShapeFunc(func_dx, shapeMat_inv, 3, i_x*h1, -maxIy*h2, h1, h2)...
    %             + getShapeFunc(func_dy, shapeMat_inv, 1, i_x*h1, -maxIy*h2+h2, h1, h2)*getShapeFunc(func_dy, shapeMat_inv, 3, i_x*h1, -maxIy*h2, h1, h2)  );
    %func2 = @(x,y)(getShapeFunc(func_dx, shapeMat_inv, 2, i_x*h1, -maxIy*h2+h2, h1, h2)*getShapeFunc(func_dx, shapeMat_inv, 4, i_x*h1, -maxIy*h2, h1, h2) ...
    %             + getShapeFunc(func_dy, shapeMat_inv, 2, i_x*h1, -maxIy*h2+h2, h1, h2)*getShapeFunc(func_dy, shapeMat_inv, 4, i_x*h1, -maxIy*h2, h1, h2)  );
    result = integral2(func1,i_x*h1, i_x*h1+h1, -maxIy*h2, -maxIy*h2+h2);   % right
    result = result + integral2(func2, i_x*h1-h1, i_x*h1, -maxIy*h2, -maxIy*h2+h2); % left
elseif i_y == j_y
    minIx = min(i_x, j_x);
    func1 = addFunctionHandle(mulFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 1, minIx*h1, -i_y*h2, h1, h2),getShapeFunc(func_dx, shapeMat_inv, 2, minIx*h1+h1, -i_y*h2, h1, h2)),...
                              mulFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 1, minIx*h1, -i_y*h2, h1, h2),getShapeFunc(func_dy, shapeMat_inv, 2, minIx*h1+h1, -i_y*h2, h1, h2)));
    func2 = addFunctionHandle(mulFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 3, minIx*h1, -i_y*h2, h1, h2),getShapeFunc(func_dx, shapeMat_inv, 4, minIx*h1+h1, -i_y*h2, h1, h2)),...
                              mulFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 3, minIx*h1, -i_y*h2, h1, h2),getShapeFunc(func_dy, shapeMat_inv, 4, minIx*h1+h1, -i_y*h2, h1, h2)));
    %func1 = @(x,y)(getShapeFunc(func_dx, shapeMat_inv, 1, minIx*h1, -i_y*h2, h1, h2)*getShapeFunc(func_dx, shapeMat_inv, 2, minIx*h1+h1, -i_y*h2, h1, h2)...
    %             + getShapeFunc(func_dy, shapeMat_inv, 1, minIx*h1, -i_y*h2, h1, h2)*getShapeFunc(func_dy, shapeMat_inv, 2, minIx*h1+h1, -i_y*h2, h1, h2)  );
    %func2 = @(x,y)(getShapeFunc(func_dx, shapeMat_inv, 3, minIx*h1, -i_y*h2, h1, h2)*getShapeFunc(func_dx, shapeMat_inv, 4, minIx*h1+h1, -i_y*h2, h1, h2) ...
    %             + getShapeFunc(func_dy, shapeMat_inv, 3, minIx*h1, -i_y*h2, h1, h2)*getShapeFunc(func_dy, shapeMat_inv, 4, minIx*h1+h1, -i_y*h2, h1, h2)  );
    result = integral2(func1, minIx*h1, minIx*h1+h1, -i_y*h2-h2, -i_y*h2);   % bottom
    result = result + integral2(func2, minIx*h1, minIx*h1+h1, -i_y*h2, -i_y*h2+h2); % top
elseif (i_x-j_x)*(i_y-j_y) < 0   % reverse diag
    maxIy = max(i_y, j_y);
    minIx = min(i_x, j_x);
    f = addFunctionHandle(mulFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 2, minIx*h1+h1, -maxIy*h2+h2, h1, h2),getShapeFunc(func_dx, shapeMat_inv, 3, minIx*h1, -maxIy*h2, h1, h2)),...
                          mulFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 2, minIx*h1+h1, -maxIy*h2+h2, h1, h2),getShapeFunc(func_dy, shapeMat_inv, 3, minIx*h1, -maxIy*h2, h1, h2)));
    %f = @(x,y)(getShapeFunc(func_dx, shapeMat_inv, 2, minIx*h1+h1, -maxIy*h2+h2, h1, h2)*getShapeFunc(func_dx, shapeMat_inv, 3, minIx*h1, -maxIy*h2, h1, h2)...
    %         + getShapeFunc(func_dy, shapeMat_inv, 2, minIx*h1+h1, -maxIy*h2+h2, h1, h2)*getShapeFunc(func_dy, shapeMat_inv, 3, minIx*h1, -maxIy*h2, h1, h2)  );
    result = integral2(f, minIx*h1, minIx*h1+h1, -maxIy*h2, -maxIy*h2+h2);
else                             % diag
    maxIy = max(i_y, j_y);
    minIx = min(i_x, j_x);
    f = addFunctionHandle(mulFunctionHandle(getShapeFunc(func_dx, shapeMat_inv, 1, minIx*h1, -maxIy*h2+h2, h1, h2),getShapeFunc(func_dx, shapeMat_inv, 4, minIx*h1+h1, -maxIy*h2, h1, h2)),...
                          mulFunctionHandle(getShapeFunc(func_dy, shapeMat_inv, 1, minIx*h1, -maxIy*h2+h2, h1, h2),getShapeFunc(func_dy, shapeMat_inv, 4, minIx*h1+h1, -maxIy*h2, h1, h2)));
    %f = @(x,y)(getShapeFunc(func_dx, shapeMat_inv, 1, minIx*h1, -maxIy*h2+h2, h1, h2)*getShapeFunc(func_dx, shapeMat_inv, 4, minIx*h1+h1, -maxIy*h2, h1, h2)...
    %         + getShapeFunc(func_dy, shapeMat_inv, 1, minIx*h1, -maxIy*h2+h2, h1, h2)*getShapeFunc(func_dy, shapeMat_inv, 4, minIx*h1+h1, -maxIy*h2, h1, h2)  );
    result = integral2(f, minIx*h1, minIx*h1+h1, -maxIy*h2, -maxIy*h2+h2);    
end