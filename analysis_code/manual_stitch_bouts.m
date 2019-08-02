function [newbouts] = manual_stitch_bouts(data,bouts)
dataMax = max(data.raw_voltage);
nGood = 0;
curBout = 1;
newbouts = struct('nlicks',{},'licks',{},'onset',{},'onset_ind',{},'offset',{},'offset_ind',{},'duration',{},'solution',{});
stitched = [];
off1 = bouts(curBout).offset_ind;
on2 = bouts(curBout+1).onset_ind;
bout1Off = bouts(curBout).offset;
bout2On = bouts(curBout+1).onset;
tempBout = struct();
oldTempBout = struct();
lastKey = '';
h=figure;
set(h,'KeyPressFcn',@KeyPressCb);
plot(data.tvec(off1-1000:on2+1000),data.raw_voltage(off1-1000:on2+1000))
hold on;
plot([bout1Off bout1Off],[0 1000],'k')
plot([bout2On bout2On],[0 1000],'g')
hold off;
function KeyPressCb(~,evnt)
    %fprintf('key pressed: %s\n',evnt.Key);
    if strcmpi(evnt.Key,'leftarrow')
        nGood=nGood+1;
        oldTempBout = tempBout;
        if (isempty(fieldnames(tempBout)))
            tempBout = bouts(curBout);
        end
        newbouts(nGood) = tempBout;
        tempBout = struct();
        curBout=curBout+1;
        lastKey = 'leftarrow';
        plotNext(curBout)
    elseif strcmpi(evnt.Key,'rightarrow')
        oldTempBout = tempBout;
        if (isempty(fieldnames(tempBout)))
            tempBout.nlicks = bouts(curBout).nlicks;
            tempBout.licks = bouts(curBout).licks;
            tempBout.onset = bouts(curBout).onset;
            tempBout.onset_ind = bouts(curBout).onset_ind;
            tempBout.offset = bouts(curBout).offset;
            tempBout.offset_ind = bouts(curBout).offset_ind;
            tempBout.duration = tempBout.offset - tempBout.onset;
            tempBout.solution = bouts(curBout).solution;
        else
            tempBout.nlicks = tempBout.nlicks + bouts(curBout).nlicks;
            tempBout.licks = [tempBout.licks bouts(curBout).licks];
            tempBout.offset = bouts(curBout).offset;
            tempBout.offset_ind = bouts(curBout).offset_ind;
            tempBout.duration = tempBout.offset - tempBout.onset;
        end
        curBout = curBout+1;
        lastKey = 'rightarrow';
        plotNext(curBout)
    elseif strcmpi(evnt.Key,'downarrow')
        curBout = curBout - 1;
        curBout = max(1,curBout);
        if (strcmp(lastKey,'leftarrow'))
            nGood = nGood - 1;
            nGood = max(0,nGood);
            newbouts = newbouts(1:end-1);
            tempBout = oldTempBout;
        elseif (strcmp(lastKey,'rightarrow'))
            tempBout = oldTempBout;
        end
        lastKey = 'downarrow';
        plotNext(curBout)
    end  
end
function plotNext(curBout)
    if (curBout > (length(bouts)-1))
        close all;
        if (~strcmp(lastKey,'leftarrow'))
            nGood = nGood + 1;
            newbouts(nGood) = tempBout;
        end
        return;
    end
    off1 = bouts(curBout).offset_ind;
    on2 = bouts(curBout+1).onset_ind;
    bout1Off = bouts(curBout).offset;
    bout2On = bouts(curBout+1).onset;
    plot(data.tvec(off1-1000:on2+1000),data.raw_voltage(off1-1000:on2+1000))
    hold on;
    plot([bout1Off bout1Off],[0 1000],'k')
    plot([bout2On bout2On],[0 1000],'g')
    ylim([0 dataMax])
    hold off;
end
uiwait;
end

