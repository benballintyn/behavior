function [ema] = getEMA(data,alpha)
ema(1) = data(1);
n = length(data);
for i=2:n
    ema(i) = alpha*data(i) + (1-alpha)*ema(i-1);
end

