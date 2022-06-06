function [offset] = findOnset(data,thresh)
%{
% data1 - ema_fast
% data2 - ema_slow
offset = 0;
for i=1:length(data1)
    if (data1(end-i) <= data2(end-i))
        break;
    else
        offset = offset+1;
    end
end
%}
offset = intmax;
for i=1:(length(data)-1)
    if (data(end-i) <= thresh)
        offset = i;
        break;
    end
end
offset = min(offset,length(data));
end

