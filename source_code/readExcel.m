
method = {'gauss', 'gauss', ..}';
error = {0.1, 0.2, ...}';
    
method{end+1} = 0.11;

data = table(method, error);
writetable(data, 'data.csv');