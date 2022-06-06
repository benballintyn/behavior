function [ind] = refineOnset(data)
ind = 0;
if (data(1) < 20)
    return;
end
for i=1:(length(data)-1)
    if (data(end-i) < data(end - i + 1))
        ind = ind+1;
    else
        break;
    end
end
end

