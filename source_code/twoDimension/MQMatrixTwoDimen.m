function [result] = GaussMatrixTwoDimen(i,j,h1,h2,x_num)
% numbering start from (0,0)
numNodes = x_num - 1;
i_x = mod(i-1, numNodes);
i_y = floor((i-1)/numNodes);
j_x = mod(j-1, numNodes);
j_y = floor((j-1)/numNodes);

if (i_x-j_x)^2 > 1 || (i_y-j_y)^2 > 1
    result = 0;
elseif (i_x-j_x)^2+(i_y-j_y)^2 == 0
 
elseif i_x == j_x
    
elseif i_y == j_y
    
else
        
end