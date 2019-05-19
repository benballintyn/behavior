function [bouts] = get_lick_bouts2(licks)

for i=1:length(licks)
    if (i == 1)
        bouts(1).onset = licks(i).onset;
        bouts(1).nlicks = 1;
        bouts(1).licks(1) = licks(i);
        curBout = 1;
    elseif (licks(i).onset - licks(i-1).offset < 5)
        bouts(curBout).nlicks = bouts(curBout).nlicks + 1;
        bouts(curBout).licks(bouts(curBout).nlicks) = licks(i);
    else
        bouts(curBout).offset = licks(i-1).offset;
        
        curBout = curBout+1;
        bouts(curBout).onset = licks(i).onset;
        bouts(curBout).nlicks = 1;
        bouts(curBout).licks(1) = licks(i);
    end
end
if (isempty(bouts(end).offset))
    bouts(end).offset = licks(end).offset;
end
end

