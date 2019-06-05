function [good_licks] = findLicks(data)
b = find(data.tvec < 300);
thresh = ceil(max(data.zlowpass(b))) + std(data.zlowpass(b));
[pks,locs,w,p] = findpeaks(data.zlowpass,data.tvec,'MinPeakHeight',thresh,'MinPeakDistance',.01);
[inds,~] = ismember(data.tvec,locs);
loc_inds = find(inds == 1);
for i=1:length(locs)
    start_ind = loc_inds(i)-20;
    end_ind = min(loc_inds(i)+20,length(data.tvec));
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
end

ngood = 0;
for i=1:length(licks)
    if (licks(i).duration < 0 || licks(i).duration > .03)
        continue;
    else
        ngood = ngood+1;
        good_licks(ngood) = licks(i);
    end
end
end

