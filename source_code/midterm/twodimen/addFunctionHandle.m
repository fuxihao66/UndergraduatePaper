function [result] = addFunctionHandle(func1, func2)

result = @(x,y)(func1(x,y)+func2(x,y));