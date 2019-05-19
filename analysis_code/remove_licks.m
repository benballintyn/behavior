function [newlicks] = remove_licks(licks,thresh,ampThresh,lowdur,highdur)
good_licks = 0;
for i=1:length(licks)
    if (licks(i).top < thresh)
        continue;
    elseif (licks(i).duration < lowdur || licks(i).duration > highdur)
        continue;
    elseif ((licks(i).top - licks(i).onset_val) < ampThresh)
        continue;
    else
        good_licks = good_licks + 1;
        newlicks(good_licks) = licks(i);
    end
end
end

