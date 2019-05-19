function [newlicks,lb,ub] = find_licks3(data,startind,thresh)
inlick = 0;
nlicks = 0;
lb = zeros(1,length(data));
ub = zeros(1,length(data));
if (startind < 1000)
    error('Starting index not late enough to use this method')
end

for i=startind:length(data)
    if (mod(i,1000) == 0)
        disp(['time = ' num2str(i)])
    end
    if ((i+1000) > length(data))
        endind = length(data);
    else
        endind = i+1000;
    end
    dataWindow = data((i-999):endind);
    if (isempty(find(dataWindow > thresh)))
        continue;
    end
    [n,edges] = histcounts(dataWindow,'BinWidth',2);
    
    [pks,locs] = findpeaks(n,'SortStr','descend','NPeaks',2);
    bounds = edges(locs);
    if (length(bounds) == 1)
        lb(i) = bounds;
        ub(i) = bounds;
    else
        lb(i) = bounds(1);
        ub(i) = bounds(2);
    end
    if (sum((data((i-4):i) > max(bounds))) == 5 && sum((data((i-4):i) > thresh)) && ~inlick)
        inlick = 1;
        nlicks = nlicks+1;
        licks(nlicks).onset = i-4;
    end
    if (sum((data((i-4):i) < min(bounds))) == 5 && sum((data((i-4):i) > thresh)) && inlick)
        inlick = 0;
        licks(nlicks).offset = i-4;
        licks(nlicks).duration = licks(nlicks).offset - licks(nlicks).onset;
    end
end
if (inlick)
    newlicks = licks(1:end-1);
else
    newlicks = licks;
end
end

