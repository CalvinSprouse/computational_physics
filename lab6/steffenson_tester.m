func = @(x) 2*exp(-x)-x;

steffenson_finder(func, 0)
fzero(func, 0)