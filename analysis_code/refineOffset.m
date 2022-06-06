function [ind] = refineOffset(data)
ind = 0;
count = 0;
for i=1:(length(data)-1)
    if (data(i) < 20)
        return;
    end
    if (data(i) > data(i + 1))
        ind = ind+1;
    else
        count = count+1;
        if (count == 20)
            return;
        end
    end
end
end

