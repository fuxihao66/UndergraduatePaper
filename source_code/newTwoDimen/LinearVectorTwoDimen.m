function [result] = LinearVectorTwoDimen(x_start, y_start,i,h1,h2,f,x_num,y_num,boundary)

numNodesX = x_num - 1;
numNodesY = y_num - 1;
i_x = mod(i-1, numNodesX)+1;
i_y = floor((i-1)/numNodesX)+1;
i_x_real = x_start + i_x*h1;
i_y_real = y_start + i_y*h2;


result = 0.0;
if i_x == 1 || i_y == 1 || i_x == numNodesX || i_y == numNodesY
    for yindex=0:y_num
    % x = 0 / x_num
        result = result - LinearMatrixTwoDimen(x_start, y_start, i_x, i_y, 0, yindex, h1,h2,x_num)*boundary(x_start,y_start+yindex*h2);           %(x_start, y_start+yindex*h2)
        result = result - LinearMatrixTwoDimen(x_start, y_start, i_x, i_y, x_num, yindex, h1,h2,x_num)*boundary(x_start+x_num*h1,y_start+yindex*h2);  %(x_end, y_start+yindex*h2)
    end
    for xindex=1:numNodesX       
        result = result - LinearMatrixTwoDimen(x_start, y_start, i_x, i_y, xindex, 0, h1,h2,x_num)* boundary(x_start+xindex*h1 ,y_start);            %(x_start+xindex*h1, y_start)
        result = result - LinearMatrixTwoDimen(x_start, y_start, i_x, i_y, xindex, y_num, h1,h2,x_num)* boundary(x_start+xindex*h1 ,y_start+y_num*h2);   %(x_start+xindex*h1, y_end)
    end
end


func1 = mulFunctionHandle(f, @(x,y)( (1-(x-i_x_real)./(h1)).*(1-(i_y_real-y)./(h2)) ) );
func2 = mulFunctionHandle(f, @(x,y)( (1-(i_x_real-x)./(h1)).*(1-(i_y_real-y)./(h2)) ) );
func3 = mulFunctionHandle(f, @(x,y)( (1-(x-i_x_real)./(h1)).*(1-(y-i_y_real)./(h2)) ) );
func4 = mulFunctionHandle(f, @(x,y)( (1-(i_x_real-x)./(h1)).*(1-(y-i_y_real)./(h2)) ) );
result = result + integral2(func1, i_x_real, i_x_real+h1, i_y_real-h2, i_y_real);
result = result + integral2(func2, i_x_real-h1, i_x_real, i_y_real-h2, i_y_real);
result = result + integral2(func3, i_x_real, i_x_real+h1, i_y_real, i_y_real+h2);
result = result + integral2(func4, i_x_real-h1, i_x_real, i_y_real, i_y_real+h2);