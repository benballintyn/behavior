function [palatibility] = getPalatability(logNaClConcentration)
f = @(x)-27*x.^2 - 67*x + 30;
fbase = f(-2);
palatibility = f(logNaClConcentration)./fbase;
end

