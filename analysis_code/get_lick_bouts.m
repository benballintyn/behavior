function [bouts] = get_lick_bouts(licks)
nbouts = 0;
inbout = 0;
for i=1:length(licks)
    if (~inbout)
        inbout = 1;
        nbouts = nbouts+1;
        bouts(nbouts).onset = licks(i).onset;
        bouts(nbouts).nlicks = 1;
        bouts(nbouts).licks(1) = licks(i);
    elseif (inbout)
        ili = licks(i).onset - licks(i-1).offset;
        if (ili > 5)
            inbout = 0;
            bouts(nbouts).offset = licks(i-1).offset;
            nbouts = nbouts + 1;
            bouts(nbouts).onset = licks(i).onset;
            inbout = 1;
        else
            bouts(nbouts).nlicks = bouts(nbouts).nlicks+1;
            bouts(nbouts).licks(bouts(nbouts).nlicks) = licks(i);
        end
        if (i == length(licks))
            bouts(nbouts).offset = licks(i).offset;
        end
    end
end
end

