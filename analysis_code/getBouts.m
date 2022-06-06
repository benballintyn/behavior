function [mini_bouts,bouts] = getBouts(licks,thresh1,thresh2)
mini_bouts(1).licks(1) = licks(1);
mini_bouts(1).nlicks = 1;
mini_bouts(1).onset = licks(1).onset;
mini_bouts(1).onset_ind = licks(1).onset_ind;
nmini = 1;
nbouts = 0;
mini_begin_ind = 1;
for i=2:length(licks)
    ili = licks(i).onset - licks(i-1).offset;
    if (ili < thresh1)
        mini_bouts(nmini).nlicks = mini_bouts(nmini).nlicks+1;
        mini_bouts(nmini).licks(mini_bouts(nmini).nlicks) = licks(i);
    elseif (ili > thresh1 && ili <= thresh2)
        mini_bouts(nmini).offset = mini_bouts(nmini).licks(end).offset;
        mini_bouts(nmini).offset_ind = mini_bouts(nmini).licks(end).offset_ind;
        mini_bouts(nmini).duration = mini_bouts(nmini).offset - mini_bouts(nmini).onset;
        mini_bouts(nmini).solution = mini_bouts(nmini).licks(end).solution;
        nmini = nmini + 1;
        mini_bouts(nmini).nlicks = 1;
        mini_bouts(nmini).licks(1) = licks(i);
        mini_bouts(nmini).onset = licks(i).onset;
        mini_bouts(nmini).onset_ind = licks(i).onset_ind;
    else
        mini_bouts(nmini).offset = mini_bouts(nmini).licks(end).offset;
        mini_bouts(nmini).offset_ind = mini_bouts(nmini).licks(end).offset_ind;
        mini_bouts(nmini).duration = mini_bouts(nmini).offset - mini_bouts(nmini).onset;
        mini_bouts(nmini).solution = mini_bouts(nmini).licks(end).solution;
        nbouts = nbouts+1;
        bouts(nbouts).mini_bouts = mini_bouts(mini_begin_ind:nmini);
        bouts(nbouts).onset = bouts(nbouts).mini_bouts(1).onset;
        bouts(nbouts).onset_ind = bouts(nbouts).mini_bouts(1).onset_ind;
        bouts(nbouts).offset = bouts(nbouts).mini_bouts(end).offset;
        bouts(nbouts).offset_ind = bouts(nbouts).mini_bouts(end).offset_ind;
        bouts(nbouts).duration = bouts(nbouts).offset - bouts(nbouts).onset;
        bouts(nbouts).solution = bouts(nbouts).mini_bouts(end).solution;
        mini_begin_ind = nmini;
        nmini = nmini + 1;
        mini_bouts(nmini).nlicks = 1;
        mini_bouts(nmini).licks(1) = licks(i);
        mini_bouts(nmini).onset = licks(i).onset;
        mini_bouts(nmini).onset_ind = licks(i).onset_ind;
    end
end
if (isempty(mini_bouts(end).offset))
    mini_bouts(end).offset = mini_bouts(end).licks(end).offset;
    mini_bouts(end).offset_ind = mini_bouts(end).licks(end).offset_ind;
end
end

