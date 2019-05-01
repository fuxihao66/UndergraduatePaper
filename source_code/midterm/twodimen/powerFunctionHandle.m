function [result] = powerFunctionHandle(func, index)

result = @(x,y)(func(x,y).^index);