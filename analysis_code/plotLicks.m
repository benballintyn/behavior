function [] = plotLicks(basedir,dataDate,animal)
dataDir = [basedir '/' dataDate '/' animal];
data = load([dataDir '/data.mat']); data=data.data;
licks = load([dataDir '/licks.mat']); licks=licks.licks;
bouts = load([dataDir '/bouts.mat']); bouts=bouts.bouts;
metadata = load([dataDir '/metadata.mat']); metadata=metadata.metadata;
for i=1:length(data)
    h=figure;
    set(h,'Position',[100 100 1400 1000])
    onsets = [licks{i}.onset];
    onsetVals = [licks{i}.onsetVal];
    offsets = [licks{i}.offset];
    offsetVals = [licks{i}.offsetVal];
    maxVals = [licks{i}.maxVal];
    maxTimes = [licks{i}.maxTime];
    maxy = max(data(i).raw_voltage);
    plot(data(i).tvec,data(i).raw_voltage);
    hold on;
    plot(onsets,onsetVals,'r*')
    plot(offsets,offsetVals,'b*')
    plot(maxTimes,maxVals,'k*')
    %{
    for i=1:length(mini_bouts)
        plot([[mini_bouts.onset];[mini_bouts.onset]],[zeros(1,length(mini_bouts));ones(1,length(mini_bouts))*maxy],'m','linew',1)
        plot([[mini_bouts.offset];[mini_bouts.offset]],[zeros(1,length(mini_bouts));ones(1,length(mini_bouts))*maxy],'m','linew',1)
    end
    %}
    for j=1:length(bouts{i})
        plot([bouts{i}(j).onset+.001;bouts{i}(j).onset+.001],[0; maxy],'k')
        plot([bouts{i}(j).offset-.001;bouts{i}(j).offset-.001],[0; maxy],'g')
    end
    title([animal ' soln = ' data(i).solution ' side = ' data(i).box_side])
    hold off;
end
end

