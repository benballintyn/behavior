function [licks] = addLicks(data,licks,savedir,allLicks,chan,varargin)
data
% addLicks Manually add licks to an existing set for one channel on one day
%       Inputs:
%           data - 1 entry from the data array for a particular day. e.g.
%                  data(1)
%
%           licks - 1 entry from the licks cell array for a particular day
%                   e.g. licks{1}
%
%           savedir - directory into which to save the modified licks file
%                     e.g. data/190801/bb8
%
%           allLicks - full licks cell array for a particular day
%
%           chan - index (1-3 depending on if it is a 2-bottle or 3-bottle
%                  task) from which the current data and licks were taken.
%                  e.g. if data and licks inputs come from data(1) and
%                  licks{1} then chan should be 1
%
%       Outputs:
%           licks - modified struct array containing prior and new licks
datalength = length(data.raw_voltage); % # of timesteps in raw_voltage
wndwSize = 2000; % size of window to display in # of timesteps
nwndws = ceil(datalength/wndwSize); % total # of windows to display
dataMax = max(data.raw_voltage); % maximum raw_voltage value

onsetLineInds = []; % array to store onset times of new licks
offsetLineInds = []; % array to store offset times of new licks
nNew = 0; % # of new licks
speedFactor = 1; % # of indices to move for left/right commands
offset = 0; % variable that stores how far left/right the user has moved the current window
horz_line_height = nan;
h = figure;
set(h,'KeyPressFcn',@KeyPressCb);
set(h,'Position',[200 200 1600 1000])
wndwNum = 1;
plotWindow(wndwNum); % plot the first window

% callback function that listens for key presses
function KeyPressCb(~,evnt)
    fprintf('key pressed: %s\n',evnt.Key);
    % move to the previous window and eliminate all new licks in this
    % window (that have not yet been accepted)
    if strcmpi(evnt.Key,'leftarrow')
        onsetLineInds = [];
        offsetLineInds = [];
        nNew = 0;
        offset = 0;
        wndwNum = max(1,wndwNum-1);
        horz_line_height = nan;
        plotWindow(wndwNum)
    % Accept all new licks in the current window and add them to the lick
    % array
    elseif strcmpi(evnt.Key,'rightarrow')
        if (length(onsetLineInds) ~= length(offsetLineInds))
            error('Unequal numbers of onset and offset indices for new licks')
        end
        onsets = data.tvec(onsetLineInds);
        onsetVals = data.raw_voltage(onsetLineInds);
        offsets = data.tvec(offsetLineInds);
        offsetVals = data.raw_voltage(offsetLineInds);
        for i=1:length(onsetLineInds)
            tmpLick.onset_ind = onsetLineInds(i);
            tmpLick.onset = onsets(i);
            tmpLick.onsetVal = onsetVals(i);
            tmpLick.offset_ind = offsetLineInds(i);
            tmpLick.offset = offsets(i);
            tmpLick.offsetVal = offsetVals(i);
            tmpLick.duration = tmpLick.offset - tmpLick.onset;
            tmpLick.raw_voltage = data.raw_voltage(tmpLick.onset_ind:tmpLick.offset_ind);
            tmpLick.smoothed_voltage = data.smoothed_voltage(tmpLick.onset_ind:tmpLick.offset_ind);
            [maxVal,maxInd] = max(tmpLick.raw_voltage);
            tmpLick.maxInd = tmpLick.onset_ind + maxInd - 1;
            tmpLick.maxTime = data.tvec(tmpLick.maxInd);
            tmpLick.maxVal = maxVal;
            tmpLick.solution = data.solution;
            tmpLick.certainty = NaN;
            tmpLick.box_side = data.box_side;
            tmpLick.channel = data.channel;
            if (isempty(licks))
                laterLickInds = [];
            else
                laterLickInds = find([licks.onset] > tmpLick.offset);
            end
            if (isempty(laterLickInds))
                licks = [licks tmpLick];
            else
                licks = [licks(1:(laterLickInds(1)-1)) tmpLick licks(laterLickInds)];
            end
        end
        onsetLineInds = [];
        offsetLineInds = [];
        nNew = 0;
        offset = 0;
        wndwNum = min(wndwNum+1,nwndws);
        plotWindow(wndwNum)
    % create the onset and offset lines for a new lick
    elseif strcmpi(evnt.Key,'n')
        nNew = nNew+1;
        if (~isempty(licks))
            meanLickDurs = [licks.offset_ind]-[licks.onset_ind];
            meanLickDur = ceil(mean(meanLickDurs(meanLickDurs < 500)));
        else
            meanLickDur = 60;
        end
        mid_ind=makeOnsetOffsetLines(wndwNum);
        onsetLineInds = [onsetLineInds mid_ind];
        offsetLineInds = [offsetLineInds mid_ind+meanLickDur];
        plotWindow(wndwNum)
    % move the current onset and offset lines speedFactor indices to the left
    elseif strcmpi(evnt.Key,'e')
        if (nNew > 0)
            onsetLineInds(nNew) = onsetLineInds(nNew)-speedFactor;
            offsetLineInds(nNew) = offsetLineInds(nNew)-speedFactor;
            plotWindow(wndwNum)
        end
    % move the current onset and offset lines speedFactor indices to the right
    elseif strcmpi(evnt.Key,'r')
        if (nNew > 0)
            onsetLineInds(nNew) = onsetLineInds(nNew)+speedFactor;
            offsetLineInds(nNew) = offsetLineInds(nNew)+speedFactor;
            plotWindow(wndwNum)
        end
    % move the current onset line speedFactor indices to the left
    elseif strcmpi(evnt.Key,'q')
        if (nNew > 0)
            onsetLineInds(nNew) = onsetLineInds(nNew)-speedFactor;
            plotWindow(wndwNum)
        end
    % move the current onset line speedFactor indices to the right
    elseif strcmpi(evnt.Key,'w')
        if (nNew > 0)
            onsetLineInds(nNew) = onsetLineInds(nNew)+speedFactor;
            plotWindow(wndwNum)
        end
    % move the current offset line speedFactor indices to the left
    elseif strcmpi(evnt.Key,'a')
        if (nNew > 0)
            offsetLineInds(nNew) = offsetLineInds(nNew)-speedFactor;
            plotWindow(wndwNum)
        end
    % move the current offset line speedFactor indices to the right
    elseif strcmpi(evnt.Key,'s')
        if (nNew > 0)
            offsetLineInds(nNew) = offsetLineInds(nNew)+speedFactor;
            plotWindow(wndwNum)
        end
    % accepts any of the number keys as input (0-9) and sets that to be the
    % new speedFactor
    elseif ~isempty(regexp(evnt.Key,'\d','ONCE'))
        speedFactor = str2double(evnt.Key);
    % shift the window left by speedFactor indices
    elseif strcmpi(evnt.Key,'z')
        offset = offset - speedFactor;
        plotWindow(wndwNum)
    % shift the window right by speedFactor indices
    elseif strcmpi(evnt.Key,'x')
        offset = offset + speedFactor;
        plotWindow(wndwNum)
    % takes the current lick array (with any new licks that have been
    % accepted) and adds them to the full lick cell array and saves that.
    % This can be used to save progress as you go along but note that this
    % WILL overwrite the preexisting lick.mat file
    elseif strcmpi(evnt.Key,'tab')
        disp('saving current progress...')
        allLicks{chan} = licks;
        licks = allLicks;
        try 
            save([savedir '/licks.mat'],'licks','-mat')
        catch
            warning('savedir not valid. Saving in ~/temp_lick_analysis/')
            if (~exist('~/temp_lick_analysis','dir'))
                mkdir('~/temp_lick_analysis')
                save('~/temp_lick_analysis/licks.mat','licks','-mat')
            end
        end
        licks = allLicks{chan};
    % Hit enter/return to finish adding licks
    elseif strcmpi(evnt.Key,'return')
        close all;
        return;
    elseif strcmpi(evnt.Key,'t')
        if (~isnan(horz_line_height))
            start_ind = (wndwNum-1)*wndwSize + 1 + offset;
            end_ind = min(datalength,wndwNum*wndwSize+offset);
            inLick = 0;
            for i = 1:(end_ind - start_ind + 1)
                if (~isempty(varargin))
                    if (data.(varargin{1})(start_ind + i - 1) > horz_line_height & ~inLick)
                        inLick = 1;
                        onsetLineInds = [onsetLineInds (start_ind + i - 2)];
                    elseif (data.(varargin{1})(start_ind + i - 1) < horz_line_height & inLick)
                        inLick = 0;
                        offsetLineInds = [offsetLineInds (start_ind + i - 1)];
                    end
                else
                    if (data.raw_voltage(start_ind + i -1) > horz_line_height & ~inLick)
                        inLick = 1;
                        onsetLineInds = [onsetLineInds (start_ind + i - 2)];
                    elseif (data.raw_voltage(start_ind + i - 1) < horz_line_height & inLick)
                        inLick = 0;
                        offsetLineInds = [offsetLineInds (start_ind + i - 1)];
                    end
                end
            end
            if (inLick)
                onsetLineInds = onsetLineInds(1:end-1);
            end
            if (length(onsetLineInds) ~= length(offsetLineInds))
                error('Different number of onsets and offsets')
            end
            lick_durs = data.tvec(offsetLineInds) - data.tvec(onsetLineInds);
            tmpOnsetLineInds = onsetLineInds;
            tmpOffsetLineInds = offsetLineInds;
            for i = 1:length(onsetLineInds)
                if (lick_durs(i) < .025)
                    tmpOnsetLineInds = setdiff(tmpOnsetLineInds,onsetLineInds(i));
                    tmpOffsetLineInds = setdiff(tmpOffsetLineInds,offsetLineInds(i));
                end
            end
            onsetLineInds = tmpOnsetLineInds;
            offsetLineInds = tmpOffsetLineInds;
            horz_line_height = nan;
        else
            horz_line_height = 100;
        end
        plotWindow(wndwNum)
    elseif strcmpi(evnt.Key,'uparrow')
        horz_line_height = horz_line_height + speedFactor;
        plotWindow(wndwNum)
    elseif strcmpi(evnt.Key,'downarrow')
        horz_line_height = horz_line_height - speedFactor;
        plotWindow(wndwNum)
    end
end

function plotWindow(wndwNum)
    start_ind = (wndwNum-1)*wndwSize + 1 + offset;
    end_ind = min(datalength,wndwNum*wndwSize+offset);
    if (~isempty(licks))
        includedLickInds = find([licks.onset_ind] >= start_ind & [licks.onset_ind] <= end_ind);
    else
        includedLickInds = [];
    end
    if (~isempty(varargin))
        plot(data.tvec(start_ind:end_ind),data.(varargin{1})(start_ind:end_ind))
    else
        plot(data.tvec(start_ind:end_ind),data.raw_voltage(start_ind:end_ind))
    end
    hold on;
    for i=1:length(includedLickInds)
        curLick = includedLickInds(i);
        lickOn = licks(curLick).onset;
        lickOff = licks(curLick).offset;
        lickMax = licks(curLick).maxTime;
        hold on;
        plot(lickOn,licks(curLick).onsetVal,'b*')
        plot(lickOff,licks(curLick).offsetVal,'r*')
        plot(lickMax,licks(curLick).maxVal,'k*')
    end
    if (length(onsetLineInds) ~= length(offsetLineInds))
        error('Different number of onsets and offsets')
    end
    tvecVals = data.tvec(onsetLineInds);
    for i=1:length(onsetLineInds)
        plot([tvecVals(i) tvecVals(i)],[0 dataMax],'b')
    end
    tvecVals = data.tvec(offsetLineInds);
    for i=1:length(offsetLineInds)
        plot([tvecVals(i) tvecVals(i)],[0 dataMax],'r')
    end
    if (~isnan(horz_line_height))
        plot([data.tvec(start_ind) data.tvec(end_ind)],[horz_line_height horz_line_height])
    end
    ylim([0 dataMax])
    hold off;
end

function mid_ind=makeOnsetOffsetLines(wndwNum)
    if isempty(offsetLineInds)
        start_ind = (wndwNum-1)*wndwSize + 1;
        end_ind = min(datalength,wndwNum*wndwSize);
        mid_ind = floor((end_ind + start_ind + 1)/2);
    else
        if (~isempty(licks))
            onsets = [licks.onset_ind];
            offsets = [licks.offset_ind];
            ilis = onsets(2:end) - offsets(1:end-1);
            meanILI = ceil(mean(ilis(ilis < 500)));
            if (isnan(meanILI) || meanILI > .5)
                mid_ind = offsetLineInds(end) + 100;
            else
                mid_ind = offsetLineInds(end)+meanILI;
            end
            disp([num2str(mid_ind) ' ' num2str(meanILI)])
        else
            mid_ind = offsetLineInds(end)+200;
        end
    end
end
uiwait;
end

