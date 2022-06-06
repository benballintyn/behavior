function [data,licks,bouts,states] = lstmLickAnalysis(basedir,dataDate,animal,redo)
% lstmLickAnalysis
%   This function uses a trained LSTM (Long-Short-Term Memory) deep network
%   to identify licks from the raw data.
%
%   Inputs:
%       basedir - top level directory containing all of the animal data
%
%       dataDate - date of the day to be analyzed. e.g. '190808'
%
%       animal - name of the animal to be analyzed. e.g. 'bb8'
%
%       redo - 'yes' or 'no'. 'yes' to reanalyze data.
%
%   Outputs:
%       data - struct array with one entry for each data channel
%
%       licks - cell array of struct arrays. Length of the cell array is
%               the same length as data. Each struct represents 1 lick on
%               that channel
%
%       bouts - cell array of struct arrays. Each struct array contains all
%               of the bouts on that channel
%
%       states - matrix of size(length(data)+1,36000) containing a
%                downsampled sequence of the 'states' of the animal during
%                that session. State == 1 indicates the animal is
%                'wandering'. State == 2 indicates the animal is at
%                solution 1 and so on.
savedir = [basedir '/' dataDate '/' animal];
if (~exist([savedir '/states.mat'],'file') || strcmp(redo,'yes'))
    if (exist([savedir '/licks.mat'],'file'))
        licks = load([savedir '/licks.mat']); licks=licks.licks;
        save([savedir '/oldLicks.mat'],'licks','-mat')
    end
    [data] = read_datafiles(basedir,dataDate,animal);
    save([savedir '/data.mat'],'data','-mat')
    [licks] = findLicksWithLSTM(data,'yes',savedir);
    save([savedir '/licks.mat'],'licks','-mat')
    for i=1:length(data)
        [bouts] = getBouts2(licks{i},1.5);
        newbouts{i} = manual_stitch_bouts(data(i),bouts);
    end
    bouts=newbouts;
    save([savedir '/bouts.mat'],'bouts','-mat')
    [states] = get_behavior_states(bouts,data(1).tvec);
    save([savedir '/states.mat'],'states','-mat')
    manual_flag.done = 1;
    manual_flag.date = date;
    save([savedir '/manual_flag.mat'],'manual_flag','-mat')
    [data] = labelDataWithChannels(basedir,dataDate,animal);
else
    disp(['previous analysis present in ' savedir '. loading...'])
    data = load([savedir '/data.mat']); data=data.data;
    if (~isfield(data,'box_side'))
        labelDataWithChannels(basedir,dataDate,animal);
    end
    licks = load([savedir '/licks.mat']); licks=licks.licks;
    bouts = load([savedir '/bouts.mat']); bouts=bouts.bouts;
    states = load([savedir '/states.mat']); states=states.states;
    save([savedir '/oldLicks.mat'],'licks','-mat')
    for i=1:length(data)
        [nl] = manualLickID(data(i),licks{i});
        newlicks{i} = nl;
        bouts{i} = getBouts2(nl,1.5);
        bouts{i} = manual_stitch_bouts(data(i),bouts{i});
    end
    [states] = get_behavior_states(bouts,data(1).tvec);
    licks = newlicks;
    save([savedir '/data.mat'],'data','-mat')
    save([savedir '/licks.mat'],'licks','-mat')
    save([savedir '/bouts.mat'],'bouts','-mat')
    save([savedir '/states.mat'],'states','-mat')
    labelDataWithChannels(basedir,dataDate,animal);
    manual_flag.done = 1;
    manual_flag.date = date;
    save([savedir '/manual_flag.mat'],'manual_flag','-mat')
end
end

