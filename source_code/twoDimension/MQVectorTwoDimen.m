function [result] = GaussVectorTwoDimen(i,h1,h2,x_num,boundary)

numNodes = x_num - 1;
i_x = mod(i-1, numNodes);
i_y = floor((i-1)/numNodes);

if i_y ~= 0
    result = 0;
else
    
end