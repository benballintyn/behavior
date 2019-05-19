function [licks] = find_licks4(signal,tvec,dvlength,endbaseline)
inlick = 0;
nlicks = 0;
timeInLick = 0;
dv = zeros(1,dvlength);
diffSignal = diff(signal(1:endbaseline));
diffTvec = diff(tvec(1:endbaseline));
meanDiff = mean(diffSignal./diffTvec);
stdDiff = std(diffSignal./diffTvec);
for i=(dvlength+1):length(signal)
    dv = [0 dv(1:end-1)];
    dv(1) = (signal(i) - signal(i-1))/(tvec(i) - tvec(i-1));
    upindicator = (dv > (meanDiff + 2*stdDiff));
    downindicator = (dv < -(meanDiff + 2*stdDiff));
    goinUp = (sum(upindicator) == dvlength);
    goinDown = (sum(downindicator) == dvlength);
    if (goinUp && ~inlick)
        inlick = 1;
        nlicks = nlicks+1;
        licks(nlicks).onset_ind = i-dvlength;
        licks(nlicks).onset = tvec(i-dvlength);
        licks(nlicks).onset_val = signal(i-dvlength);
    end
    if (goinDown && inlick)
        inlick = 0;
        timeInLick = 0;
        licks(nlicks).offset_ind = i;
        licks(nlicks).data = signal(licks(nlicks).onset_ind:licks(nlicks).offset_ind);
        licks(nlicks).offset = tvec(i);
        licks(nlicks).offset_val = signal(i);
        licks(nlicks).duration = licks(nlicks).offset - licks(nlicks).onset;
        nTsteps = licks(nlicks).offset_ind - licks(nlicks).onset_ind;
        licks(nlicks).top = signal(licks(nlicks).onset_ind + round(nTsteps/2));
        if (licks(nlicks).duration > .15 || licks(nlicks).duration < .05)
            licks = licks(1:(nlicks-1));
            nlicks = nlicks-1;
        end
    end
end
if (inlick)
    licks = licks(1:end-1);
end
end

