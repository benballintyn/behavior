function [good_licks] = findLicks2(data)
% findLicks2
%   Finds crossings of a voltage threshold in the raw voltage recordings
%   and marks the onset/offset of licks
%   Inputs:
%       data - structure array returned from read_datafiles
%
%   Outputs:
%       good_licks - cell array, of same length as data, of struct arrays
%       with each struct being a lick
b = find(data.tvec < 300);
min_thresh = ceil(mean(data.raw_voltage(b))) + 10*std(data.raw_voltage(b)); % had at 20std
ema_thresh = getEMA(data.smoothed_voltage,.05);
inlick = 0;
nlicks=0;
for i=1:length(data.raw_voltage)
    thresh = max(min_thresh,ema_thresh(i));
    if (~inlick && data.smoothed_voltage(i) > thresh)
        nlicks = nlicks+1;
        licks(nlicks).onset = data.tvec(i-1);
        licks(nlicks).onset_ind = i-1;
        licks(nlicks).onsetVal = data.raw_voltage(i-1);
        inlick = 1;
    elseif (inlick && data.smoothed_voltage(i) < thresh)
        licks(nlicks).offset = data.tvec(i);
        licks(nlicks).offset_ind = i;
        licks(nlicks).offsetVal = data.raw_voltage(i);
        licks(nlicks).duration = licks(nlicks).offset - licks(nlicks).onset;
        inlick = 0;
    end
end
[licks] = refineOnsetsOffsets(licks,data);
for i=1:length(licks)
    licks(i).raw_voltage = data.raw_voltage(licks(i).onset_ind:licks(i).offset_ind);
    licks(i).smoothed_voltage = data.smoothed_voltage(licks(i).onset_ind:licks(i).offset_ind);
    [maxVal,maxInd] = max(licks(i).raw_voltage);
    licks(i).maxInd = licks(i).onset_ind + maxInd;
    licks(i).maxTime = data.tvec(licks(i).maxInd);
    licks(i).maxVal = data.raw_voltage(licks(i).maxInd);
    licks(i).solution = data.solution;
end
ngood = 0;
upper_thresh = min(.12,mean([licks.duration]) + 4*std([licks.duration]));
for i=1:length(licks)
    if (isempty(licks(i).offset))
        continue;
    end
    if (licks(i).duration < .01 || licks(i).duration > upper_thresh)
        continue;
    elseif ((licks(i).maxVal - licks(i).onsetVal) < 20)
        continue;
    else
        ngood=ngood+1;
        good_licks(ngood) = licks(i);
    end
end

if (~exist('good_licks','var'))
    good_licks=[];
else
    ngood = 1;
    temp(ngood) = good_licks(1);
    for i=2:length(good_licks)
        if ((good_licks(i).onset - good_licks(i-1).offset) < .005)
            continue;
        else
            ngood = ngood+1;
            temp(ngood) = good_licks(i);
        end
    end
    good_licks = temp;
end
end

