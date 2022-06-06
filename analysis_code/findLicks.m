function [good_licks] = findLicks(data)
b = find(data.tvec < 300);
thresh = ceil(mean(data.raw_voltage(b))) + 20*std(data.raw_voltage(b))
[pks,locs,w,p] = findpeaks(data.raw_voltage,data.tvec,'MinPeakHeight',thresh,'MinPeakDistance',.06);
[inds,~] = ismember(data.tvec,locs);
loc_inds = find(inds == 1);
licks = [];
for i=1:length(locs)
    start_ind = max(1,loc_inds(i)-100);
    end_ind = min(loc_inds(i)+100,length(data.tvec));
    licks(i).onset_ind = loc_inds(i) - findOnset(data.raw_voltage(start_ind:loc_inds(i)),thresh);
    licks(i).offset_ind = loc_inds(i) + findOffset(data.raw_voltage(loc_inds(i):end_ind),thresh);
    licks(i).onset = data.tvec(licks(i).onset_ind);
    licks(i).offset = data.tvec(licks(i).offset_ind);
    licks(i).onsetVal = data.raw_voltage(licks(i).onset_ind);
    licks(i).offsetVal = data.raw_voltage(licks(i).offset_ind);
    licks(i).duration = licks(i).offset - licks(i).onset;
    licks(i).npts = licks(i).offset_ind - licks(i).onset_ind;
    licks(i).raw_voltage = data.raw_voltage(licks(i).onset_ind:licks(i).offset_ind);
    licks(i).zraw = data.zraw(licks(i).onset_ind:licks(i).offset_ind);
    licks(i).lowpass_voltage = data.lowpass_voltage(licks(i).onset_ind:licks(i).offset_ind);
    licks(i).zlowpass = data.zlowpass(licks(i).onset_ind:licks(i).offset_ind);
    [maxVal,maxInd] = max(licks(i).raw_voltage);
    licks(i).maxInd = licks(i).onset_ind + maxInd - 1;
    licks(i).maxTime = data.tvec(licks(i).maxInd);
    licks(i).maxVal = maxVal;
    licks(i).solution = data.solution;
end
ngood = 0;
nbad = 0;
if (~isempty(licks))
    for i=1:length(licks)
        if (licks(i).duration < .01 || licks(i).duration > .08)
            nbad = nbad+1;
            continue;
        else
            ngood = ngood+1;
            good_licks(ngood) = licks(i);
        end
    end
else
    good_licks = [];
end
length(good_licks)
end

