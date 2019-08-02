function [data,licks,bouts,states] = analysis_wrapper(basedir,date,animal,redo)
savedir = [basedir '/' date '/' animal];
if (~exist([savedir '/states.mat'],'file') || strcmp(redo,'yes'))
    [data] = read_datafiles(basedir,date,animal);
    for i=1:length(data)
        %licks{i} = findLicks(data(i));
        licks{i} = findLicks2(data(i));
        if (~isempty(licks{i}))
            ilis{i} = getILIs(licks{i});
            %bout_thresh = mean(ilis{i}(ilis{i} < 2)) + 5*std(ilis{i}(ilis{i}<2));
            bout_thresh = 1.5;
            [bigBouts] = getBouts2(licks{i},bout_thresh);
            bouts{i} = bigBouts;
        else
            bouts{i} = [];
        end
    end
    [states] = get_behavior_states(bouts,data(1).tvec);
    save([savedir '/data.mat'],'data','-mat')
    save([savedir '/licks.mat'],'licks','-mat')
    save([savedir '/bouts.mat'],'bouts','-mat')
    save([savedir '/states.mat'],'states','-mat')
else
    disp(['previous analysis present in ' savedir '. loading...'])
    data = load([savedir '/data.mat']); data=data.data;
    licks = load([savedir '/licks.mat']); licks=licks.licks;
    bouts = load([savedir '/bouts.mat']); bouts=bouts.bouts;
    states = load([savedir '/states.mat']); states=states.states;
end
end

