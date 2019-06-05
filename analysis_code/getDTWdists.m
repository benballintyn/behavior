function [distsFilt,distsRaw] = getDTWdists(licks)
for i=1:length(licks)
    for j=1:length(licks)
        distsFilt(i,j) = dtw(licks(i).lowpass_voltage,licks(j).lowpass_voltage);
        distsRaw(i,j) = dtw(licks(i).raw_voltage,licks(j).raw_voltage);
    end
    disp([num2str((i/length(licks))*100) '%'])
end

end

