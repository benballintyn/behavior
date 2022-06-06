function [cdfvals,x] = mycdf(vals,rangemin,rangemax,rangeinterval)
x = rangemin:rangeinterval:rangemax;
sortedVals = sort(vals);
N = length(sortedVals);
cdfvals = zeros(1,length(x));
for i=1:length(x)
    cdfvals(i) = sum(sortedVals < x(i))/N;
end
end

