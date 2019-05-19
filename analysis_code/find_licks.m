function [licks] = find_licks(data,tracetype,minpeakwidth,minpeakheight)
for n=1:length(data)
    [pks,locs,w] = findpeaks(data(n).(tracetype),'MinPeakWidth',minpeakwidth,'MinPeakHeight',minpeakheight);
    for i=1:length(pks)
        licks{n}(i).time = locs(i);
        licks{n}(i).duration = w(i);
        licks{n}(i).onset = locs(i) - w(i)/2;
        licks{n}(i).offset = locs(i) + w(i)/2;
    end
end
end

