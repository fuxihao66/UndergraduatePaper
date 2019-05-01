function [result] = getShapeFunc_one(func, phi, shapeIndex, node_x, node_y, h1, h2)   % shapeIndex = 1,2,3 or 4
% get pu/px
newh2 = -h2;    % 
xx = [[0,0];[0,0];[0,0];[0,0]];
if shapeIndex == 1
    xx(1,:) = [node_x, node_y];
    xx(2,:) = [node_x + h1, node_y];
    xx(3,:) = [node_x, node_y + newh2];
    xx(4,:) = [node_x + h1, node_y + newh2];
elseif shapeIndex == 2
    xx(1,:) = [node_x - h1, node_y];
    xx(2,:) = [node_x, node_y];
    xx(3,:) = [node_x - h1, node_y + newh2];
    xx(4,:) = [node_x, node_y + newh2];
elseif shapeIndex == 3
    xx(1,:) = [node_x, node_y - newh2];
    xx(2,:) = [node_x + h1, node_y - newh2];
    xx(3,:) = [node_x, node_y];
    xx(4,:) = [node_x + h1, node_y];
else
    xx(1,:) = [node_x - h1, node_y - newh2];
    xx(2,:) = [node_x, node_y - newh2];
    xx(3,:) = [node_x - h1, node_y];
    xx(4,:) = [node_x, node_y];    
end

relaPoint = [0, 0; h1,0; 0, -h2; h1,-h2];

shapeMat = zeros(7,7);
for index_1 = 1:4
    for index_2 = 1:4
        shapeMat(index_1,index_2) = phi(norm( relaPoint([index_1],[1,2])-relaPoint([index_2],[1,2]) ));
    end
end
for index = 1:4
    shapeMat(index,5) = 1;
    shapeMat(5,index) = 1;
end
for index = 1:4
    shapeMat(index,6) = xx(index,1);
    shapeMat(6,index) = xx(index,1);
end
for index = 1:4
    shapeMat(index,7) = xx(index,2);
    shapeMat(7,index) = xx(index,2);
end
mat = inv(shapeMat);




result =  @(x,y)(  func(x-xx(1,1),y-xx(1,2)).*mat(1,shapeIndex)+ func(x-xx(2,1),y-xx(2,2)).*mat(2,shapeIndex)...
                 + func(x-xx(3,1),y-xx(3,2)).*mat(3,shapeIndex)+ func(x-xx(4,1),y-xx(4,2)).*mat(4,shapeIndex)...
                 + mat(5,shapeIndex) + mat(6, shapeIndex).*x + mat(7, shapeIndex).*y );
             
             
             