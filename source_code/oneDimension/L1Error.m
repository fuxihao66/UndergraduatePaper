function [result] = L1Error(accurate, numerical)

[w, h ] = size(accurate);
len = 0;
if w == 1
    len = h;
else
    len = w;
end
result = 0;
for i = 1:len
    result = result + abs(accurate(i) - numerical(i));
end

