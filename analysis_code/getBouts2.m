function [bouts] = getBouts2(licks,thresh)
% getBouts2
%   Inputs:
%       licks - licks from one channel. e.g. licks{1} if licks is the cell
%               array returned from analysis_wrapper.
%
%       thresh - maximum amount of time, in seconds, that licks can be
%                separated by and still be in the same bout. I have been
%                using 1.5s as the thresh but allowing bouts of separation
%                <2s to be joined.
%
%   Outputs:
%       bouts - struct array containing all the bouts found from the licks
%               input
if (isempty(licks))
    bouts = [];
    return
else
    nbouts = 1;
    bouts(1).nlicks = 1;
    bouts(1).licks(1) = licks(1);
    bouts(1).onset = licks(1).onset;
    bouts(1).onset_ind = licks(1).onset_ind;
    if (length(licks) == 1)
        bouts(1).offset = licks(1).offset;
        bouts(1).offset_ind = licks(1).offset_ind;
        bouts(1).duration = bouts(1).offset - bouts(1).onset;
        bouts(1).solution = licks(1).solution;
        return;
    else
        for i=2:length(licks)
            if ((licks(i).onset - licks(i-1).offset) > thresh)
                bouts(nbouts).offset = bouts(nbouts).licks(end).offset;
                bouts(nbouts).offset_ind = bouts(nbouts).licks(end).offset_ind;
                bouts(nbouts).duration = bouts(nbouts).offset - bouts(nbouts).onset;
                bouts(nbouts).solution = bouts(nbouts).licks(end).solution;

                nbouts = nbouts+1;
                bouts(nbouts).nlicks = 1;
                bouts(nbouts).licks(1) = licks(i);
                bouts(nbouts).onset = licks(i).onset;
                bouts(nbouts).onset_ind = licks(i).onset_ind;
            else
                bouts(nbouts).nlicks = bouts(nbouts).nlicks + 1;
                bouts(nbouts).licks(bouts(nbouts).nlicks) = licks(i);
            end
            if (i == length(licks))
                bouts(nbouts).offset = licks(i).offset;
                bouts(nbouts).offset_ind = licks(i).offset_ind;
                bouts(nbouts).duration = bouts(nbouts).offset - bouts(nbouts).onset;
                bouts(nbouts).solution = bouts(nbouts).licks(end).solution;
            end
        end
    end
end
end

