function [data,newlicks,bouts,states] = manual_wrapper(basedir,date,animal,redo)
savedir = [basedir '/' date '/' animal];
if (~exist([savedir '/states.mat'],'file') || strcmp(redo,'yes'))
    [data,licks,bouts,states] = analysis_wrapper(basedir,date,animal,'yes');
    for i=1:length(data)
        [nl] = manualLickID(data(i),licks{i});
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
        [nl] = manualLickID(data(i),licks{i});
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

