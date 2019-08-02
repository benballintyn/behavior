function [newlicks] = manualLickID(data,licks)
dataMax = max(data.raw_voltage);
nGood = 0;
curLick = 1;
dataOffset = 1000;
good_licks = [];
h=figure;
set(h,'KeyPressFcn',@KeyPressCb);
on = max(1,licks(curLick).onset_ind-dataOffset);
off = min(length(data.tvec),licks(curLick).offset_ind+dataOffset);
lickOn = licks(curLick).onset;
lickOff = licks(curLick).offset;
lickMax = licks(curLick).maxTime;
plot(data.tvec(on:off),data.raw_voltage(on:off))
hold on;
plot(lickOn,licks(curLick).onsetVal,'b*')
plot(lickOff,licks(curLick).offsetVal,'r*')
plot(lickMax,licks(curLick).maxVal,'k*')
hold off;
function KeyPressCb(~,evnt)
    %fprintf('key pressed: %s\n',evnt.Key);
    if strcmpi(evnt.Key,'leftarrow')
        curLick = curLick+1;
        plotNext(curLick)
    elseif strcmpi(evnt.Key,'rightarrow')
        nGood = nGood + 1;
        newlicks(nGood) = licks(curLick);
        curLick = curLick+1;
        good_licks = [good_licks curLick];
        plotNext(curLick)
    elseif strcmpi(evnt.Key,'downarrow')
        curLick = curLick-1;
        if (ismember(curLick,good_licks))
            nGood = nGood - 1;
            newlicks = newlicks(1:end-1);
            good_licks = good_licks(good_licks~=curLick);
        end
        plotNext(curLick)
    end  
end
function plotNext(curLick)
    if (curLick > length(licks))
        close all;
        return;
    end
    on = max(1,licks(curLick).onset_ind-dataOffset);
    off = min(length(data.tvec),licks(curLick).offset_ind+dataOffset);
    lickOn = licks(curLick).onset;
    lickOff = licks(curLick).offset;
    lickMax = licks(curLick).maxTime;
    plot(data.tvec(on:off),data.raw_voltage(on:off))
    hold on;
    plot(lickOn,licks(curLick).onsetVal,'b*')
    plot(lickOff,licks(curLick).offsetVal,'r*')
    plot(lickMax,licks(curLick).maxVal,'k*')
    ylim([0 dataMax])
    hold off;
end
uiwait;
end

