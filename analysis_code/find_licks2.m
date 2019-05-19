function [licks] = find_licks2(data,startind,thresh)
inlick = 0;
nlicks = 0;
if (startind < 5)
    error('Starting index not late enough to use this method')
end
for i=startind:length(data)
    if (sum((data((i-4):i) > thresh)) == 5 && ~inlick)
        inlick = 1;
        nlicks = nlicks+1;
        licks(nlicks).onset = i-4;
    end
    if (sum((data((i-4):i) < thresh)) == 5 && inlick)
        inlick = 0;
        licks(nlicks).offset = i-4;
        licks(nlicks).duration = licks(nlicks).offset - licks(nlicks).onset;
    end
end

end

