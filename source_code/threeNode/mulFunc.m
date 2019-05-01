function [result] = mulFunc(func1,func2)

result = @(x)(func1(x).*func2(x));