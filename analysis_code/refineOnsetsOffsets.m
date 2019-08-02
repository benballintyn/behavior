function [licks] = refineOnsetsOffsets(licks,data)
for i=1:length(licks)
    licks(i).onset_ind = licks(i).onset_ind - refineOnset(data.smoothed_voltage(licks(i).onset_ind-10:licks(i).onset_ind));
    end_ind = min(licks(i).offset_ind+50,length(data.smoothed_voltage));
    licks(i).offset_ind = licks(i).offset_ind + refineOffset(data.smoothed_voltage(licks(i).offset_ind:end_ind));
    licks(i).onset = data.tvec(licks(i).onset_ind);
    licks(i).offset = data.tvec(licks(i).offset_ind);
    licks(i).onsetVal = data.raw_voltage(licks(i).onset_ind);
    licks(i).offsetVal = data.raw_voltage(licks(i).offset_ind);
end
end

