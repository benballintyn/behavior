function [offset] = findOffset(data,thresh)
%{
% data1 - ema_fast
% data2 - ema_slow
offset = 0;
for i=1:length(data1)
    if (data1(i) <= data2(i))
        break;
    else
        offset = offset+1;
    end
end
%}
offset = intmax;
for i=1:length(data)
    if (data(i) <= thresh)
        offset = i;
        break;
    end
end
offset = min(offset,length(data));
end

