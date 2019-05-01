function [result] = getShapeFunc_const(func, mat, shapeIndex, node_x, node_y, h1, h2)   % shapeIndex = 1,2,3 or 4
% get pu/px
newh2 = -h2;    % 
if shapeIndex == 1
    x1 = [node_x, node_y];
    x2 = [node_x + h1, node_y];
    x3 = [node_x, node_y + newh2];
    x4 = [node_x + h1, node_y + newh2];
elseif shapeIndex == 2
    x1 = [node_x - h1, node_y];
    x2 = [node_x, node_y];
    x3 = [node_x - h1, node_y + newh2];
    x4 = [node_x, node_y + newh2];
elseif shapeIndex == 3
    x1 = [node_x, node_y - newh2];
    x2 = [node_x + h1, node_y - newh2];
    x3 = [node_x, node_y];
    x4 = [node_x + h1, node_y];
else
    x1 = [node_x - h1, node_y - newh2];
    x2 = [node_x, node_y - newh2];
    x3 = [node_x - h1, node_y];
    x4 = [node_x, node_y];    
end

result =  @(x,y)(  func(x-x1(1),y-x1(2)).*mat(1,shapeIndex)+ func(x-x2(1),y-x2(2)).*mat(2,shapeIndex)...
                 + func(x-x3(1),y-x3(2)).*mat(3,shapeIndex)+ func(x-x4(1),y-x4(2)).*mat(4,shapeIndex)...
                 + mat(5,shapeIndex));