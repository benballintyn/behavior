function [newlicks] = manualLickID(data,licks,savedir,allLicks,chan)
% manualLickID
%   Inputs:
%       data - one entry in the structure array returned from read_datafiles
%              e.g. data(1)
%
%       licks - struct array holding all the licks from one channel
%       lick
%
%   Outputs:
%       newlicks - modified struct array with some licks removed manually
%   Key commands:
%       <- (leftarrow): reject a lick
%       -> (rightarrow): accept a lick
%       downarrow: go back a lick (and remove it from accepted licks if it
%                  was accepted
%       q: shift lick onset back by 1 timestep
%       w: shift lick onset forward by 1 timestep
%       a: shift lick offset back by 1 timestep
%       s: shift lick offset forward by 1 timestep
%
if (isempty(licks))
    newlicks = [];
    return;
end
dataMax = max(data.raw_voltage);
nGood = 0;
curLick = 1;
dataOffset = 1000;
good_licks = [];
offsetAdjustment = 0;
h=figure;
set(h,'KeyPressFcn',@KeyPressCb);
set(h,'Position',[200 200 1600 1000])
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
if (isfield(licks(curLick),'certainty'))
    title(['Certainty = ' num2str(licks(curLick).certainty)])
end
function KeyPressCb(~,evnt)
    %fprintf('key pressed: %s\n',evnt.Key);
    if strcmpi(evnt.Key,'leftarrow')
        curLick = curLick+1;
        plotNext(curLick)
    elseif strcmpi(evnt.Key,'rightarrow')
        licks(curLick).duration = licks(curLick).offset - licks(curLick).onset;
        licks(curLick).raw_voltage = data.raw_voltage(licks(curLick).onset_ind:licks(curLick).offset_ind);
        licks(curLick).smoothed_voltage = data.smoothed_voltage(licks(curLick).onset_ind:licks(curLick).offset_ind);
        [maxVal,maxInd] = max(licks(curLick).raw_voltage);
        licks(curLick).maxInd = licks(curLick).onset_ind + maxInd - 1;
        licks(curLick).maxTime = data.tvec(licks(curLick).maxInd);
        licks(curLick).maxVal = data.raw_voltage(licks(curLick).maxInd);
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
    elseif strcmpi(evnt.Key,'q')
        licks(curLick).onset_ind = licks(curLick).onset_ind-1;
        licks(curLick).onset = data.tvec(licks(curLick).onset_ind);
        licks(curLick).onsetVal = data.raw_voltage(licks(curLick).onset_ind);
        plotNext(curLick)
    elseif strcmpi(evnt.Key,'w')
        licks(curLick).onset_ind = licks(curLick).onset_ind+1;
        licks(curLick).onset = data.tvec(licks(curLick).onset_ind);
        licks(curLick).onsetVal = data.raw_voltage(licks(curLick).onset_ind);
        plotNext(curLick)
    elseif strcmpi(evnt.Key,'z')
        licks(curLick).maxInd = licks(curLick).maxInd-1;
        licks(curLick).maxTime = data.tvec(licks(curLick).maxInd);
        licks(curLick).maxVal = data.raw_voltage(licks(curLick).maxInd);
        plotNext(curLick)
    elseif strcmpi(evnt.Key,'x')
        licks(curLick).maxInd = licks(curLick).maxInd+1;
        licks(curLick).maxTime = data.tvec(licks(curLick).maxInd);
        licks(curLick).maxVal = data.raw_voltage(licks(curLick).maxInd);
        plotNext(curLick)
    elseif strcmpi(evnt.Key,'a')
        licks(curLick).offset_ind = licks(curLick).offset_ind-1;
        licks(curLick).offset = data.tvec(licks(curLick).offset_ind);
        licks(curLick).offsetVal = data.raw_voltage(licks(curLick).offset_ind);
        plotNext(curLick)
    elseif strcmpi(evnt.Key,'s')
        licks(curLick).offset_ind = licks(curLick).offset_ind+1;
        licks(curLick).offset = data.tvec(licks(curLick).offset_ind);
        licks(curLick).offsetVal = data.raw_voltage(licks(curLick).offset_ind);
        plotNext(curLick)
    elseif strcmpi(evnt.Key,'tab')
        disp('saving current progress')
        allLicks{chan} = newlicks;
        oldLicks = licks;
        licks = allLicks;
        save([savedir '/licks.mat'],'licks','-mat')
        licks = oldLicks;
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
    if (isfield(licks(curLick),'certainty'))
        title({['Duration = ' num2str(licks(curLick).duration)], ...
                ['Certainty = ' num2str(licks(curLick).certainty)], ...
                ['Lick ' num2str(curLick) ' / ' num2str(length(licks))]})
    else
        title({['Duration = ' num2str(licks(curLick).duration)],...
               ['Lick ' num2str(curLick) ' / ' num2str(length(licks))]})
    end
end
uiwait;
end

