function [data,newlicks,bouts,states] = manual_wrapper(basedir,date,animal,redo)
% Identify licks using semi-automated process
%   Inputs:
%       basedir - path to directory containing data files
%                 e.g. 'data' or 'analyzed_data'
%
%       date - date (year,month,day) of data to be analyzed e.g. '190801'
%
%       animal - animal ID to be analyzed e.g. 'bb8'
%
%       redo - 'yes', 'no', or 'continue'. 'yes' in order to go through the the full
%       analysis process straight from raw data. 'no' to load previously
%       analyzed data and licks and manually correct them.
%
%   Outputs:
%       data - struct array containing the raw data for each channel for 1
%              day
%
%       newlicks - modified licks stored in a cell array of struct arrays of 
%                  the same length as data (1 entry for each channel)
%
%       bouts - modified bouts stored as a cell array of struct arrays
%
%       states - this is a matrix of size (Nbottles+1)x36000
savedir = [basedir '/' date '/' animal];
if (~exist([savedir '/states.mat'],'file') || strcmp(redo,'yes'))
    [data,licks,bouts,states] = analysis_wrapper(basedir,date,animal,'yes');
    for i=1:length(data)
        [nl] = manualLickID(data(i),licks{i},savedir,licks,i);
        newlicks{i} = nl;
        tempBout = getBouts2(nl,1.5);
        newbouts = manual_stitch_bouts(data(i),tempBout);
        bouts{i} = newbouts;
    end
    [states] = get_behavior_states(bouts,data(1).tvec);
    licks = newlicks;
    save([savedir '/data.mat'],'data','-mat')
    save([savedir '/licks.mat'],'licks','-mat')
    save([savedir '/bouts.mat'],'bouts','-mat')
    save([savedir '/states.mat'],'states','-mat')
    manual_flag.done = 1;
    manual_flag.date = date;
    save([savedir '/manual_flag.mat'],'manual_flag','-mat')
else
    disp(['previous analysis present in ' savedir '. loading...'])
    data = load([savedir '/data.mat']); data=data.data;
    licks = load([savedir '/licks.mat']); licks=licks.licks;
    bouts = load([savedir '/bouts.mat']); bouts=bouts.bouts;
    states = load([savedir '/states.mat']); states=states.states;
    for i=1:length(data)
        [nl] = manualLickID(data(i),licks{i},savedir,licks,i);
        newlicks{i} = nl;
        bouts{i} = getBouts2(nl,1.5);
    end
    [states] = get_behavior_states(bouts,data(1).tvec);
    licks = newlicks;
    save([savedir '/data.mat'],'data','-mat')
    save([savedir '/licks.mat'],'licks','-mat')
    save([savedir '/bouts.mat'],'bouts','-mat')
    save([savedir '/states.mat'],'states','-mat')
    manual_flag.done = 1;
    manual_flag.date = date;
    save([savedir '/manual_flag.mat'],'manual_flag','-mat')
end
end

