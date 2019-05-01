function [result] = GaussVectorTwoDimen(i,h1,h2,x_num,boundary)

numNodes = x_num - 1;
i_x = mod(i-1, numNodes)+1;
i_y = floor((i-1)/numNodes)+1;

result = 0.0;
if i_y == 1
    for index=1:numNodes
        j = index - numNodes;
        result = result - GaussMatrixTwoDimen(i,j,h1,h2,x_num)*boundary(index);
    end
end