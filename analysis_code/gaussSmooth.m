function [smoothed_data] = gaussSmooth(data,wsize)
w = gausswin(wsize); w = w./sum(w);
smoothed_data(1:(ceil(wsize/2)-1)) = data(1:(ceil(wsize/2)-1));
for i=ceil(wsize/2):(length(data)-floor(wsize/2))
    smoothed_data(i) = sum(w.*data((i-floor(wsize/2)):(i+floor(wsize/2))));
end
smoothed_data(end+1:end+floor(wsize/2)) = data((end-floor(wsize/2)+1):end);
end

