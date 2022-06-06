function [] = redoBouts(basedir,dataDate,animal)
dataDir = [basedir '/' dataDate '/' animal];
licks = load([dataDir '/licks.mat']); licks=licks.licks;
data = load([dataDir '/data.mat']); data=data.data;
for i=1:length(licks)
    [bouts{i}] = getBouts2(licks{i},2.0);
end
for i=1:length(bouts)
    newbouts{i} = manual_stitch_bouts(data(i),bouts{i});
end
bouts = newbouts;
save([dataDir '/bouts.mat'],'bouts','-mat')
[states] = get_behavior_states(bouts,data(1).tvec);
save([dataDir '/states.mat'],'states','-mat')
end

