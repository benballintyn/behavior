function [] = addLicksAndCorrect(basedir,date,animal,varargin)
%   addLicksAndCorrect Go through all channels for a day and manually add
%                      licks
%
%   Inputs:
%       basedir - path to directory containing data files
%                 e.g. 'data' or 'analyzed_data'
%
%       date - date (year,month,day) of data to be analyzed e.g. '190801'
%
%       animal - animal ID to be analyzed e.g. 'bb8'
savedir = [basedir '/' date '/' animal]; % directory that data exists and will be saved
if (~exist([savedir '/licks.mat'],'file'))
    error('No lick file to correct')
else
    % load data and previously identified licks
    licks = load([savedir '/licks.mat']); licks = licks.licks;
    data = load([savedir '/data.mat']); data=data.data;
end
% For each channel manually add licks
for i=1:length(licks)
    if (~isempty(varargin))
        licks{i} = addLicks(data(i),licks{i},savedir,licks,i,varargin{1});
    else
        licks{i} = addLicks(data(i),licks{i},savedir,licks,i);
    end
end
save([savedir '/licks.mat'],'licks','-mat') % save manually modified licks
% For each channel, identify new bouts based on modified licks and manually
% check them
for i=1:length(licks)
    tempBout = getBouts2(licks{i},1.5);
    newbouts = manual_stitch_bouts(data(i),tempBout);
    bouts{i} = newbouts;
end
save([savedir '/bouts.mat'],'bouts','-mat') % save modified bouts
[states] = get_behavior_states(bouts,data(1).tvec);
save([savedir '/states.mat'],'states','-mat')
manual_flag.done = 1;
manual_flag.date = date;
save([savedir '/manual_flag.mat'],'manual_flag','-mat')
end

