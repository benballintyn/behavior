% volen retreat figures
basedir = 'data';
dates = {'190629','190630','190701','190702','190703','190704','190705',...
         '190808','190809','190810','190811','190812','190813','190814'};
rats = {'bb8','bb9'};
for i=1:length(dates)
    truncDates{i} = [dates{i}(3:4) '/' dates{i}(5:end)];
end
% get rat weights from metadata
for i=1:length(dates)
    dataDir = [basedir '/' dates{i} '/bb8'];
    metadata = load([dataDir '/metadata.mat']); metadata = metadata.metadata;
    weights8(i) = metadata.weight;
    dataDir = [basedir '/' dates{i} '/bb9'];
    metadata = load([dataDir '/metadata.mat']); metadata = metadata.metadata;
    weights9(i) = metadata.weight;
end

% extract bouts and palatability data
[allBouts8,meanAlts8,maxAlts8,altDifs8] = combineBouts(basedir,dates,rats{1});
[allBouts9,meanAlts9,maxAlts9,altDifs9] = combineBouts(basedir,dates,rats{2});
%[meanCorr8,meanPval8,maxCorr8,maxPval8,values8] = getBoutLengthCorrelations(allBouts8);
%[meanCorr9,meanPval9,maxCorr9,maxPval9,values9] = getBoutLengthCorrelations(allBouts9);
[allBoutsBySolution8] = getAllBoutsBySolution(basedir,dates,rats{1}); % H2O, A, B, C
[allBoutsBySolution9] = getAllBoutsBySolution(basedir,dates,rats{2});

% remove bouts with single lick
nlicks8 = [allBouts8.nlicks];
moreThan1Lick8 = find(nlicks8 > 1);
nlicks9 = [allBouts9.nlicks];
moreThan1Lick9 = find(nlicks9 > 1);
allBouts8 = allBouts8(moreThan1Lick8);
meanAlts8 = meanAlts8(moreThan1Lick8);
maxAlts8 = maxAlts8(moreThan1Lick8);
altDifs8 = altDifs8(moreThan1Lick8);
allBouts9 = allBouts9(moreThan1Lick9);
meanAlts9 = meanAlts9(moreThan1Lick9);
maxAlts9 = maxAlts9(moreThan1Lick9);
altDifs9 = altDifs9(moreThan1Lick9);
for i=1:length(allBoutsBySolution8)
    nlicks = [allBoutsBySolution8{i}.nlicks];
    moreThan1 = find(nlicks > 1);
    allBoutsBySolution8{i} = allBoutsBySolution8{i}(moreThan1);
end
for i=1:length(allBoutsBySolution9)
    nlicks = [allBoutsBySolution9{i}.nlicks];
    moreThan1 = find(nlicks > 1);
    allBoutsBySolution9{i} = allBoutsBySolution9{i}(moreThan1);
end

% Plot all bout durations by solution for bb8 and bb9 separately and
% together
for i=1:4
    v8{i} = [allBoutsBySolution8{i}.duration];
end
figure;
violin(v8,'facecolor','b'); set(gca,'xtick',1:4,'xticklabel',{'H2O','.01M','.1M','1M'})
xlabel('Solution'); ylabel('Bout duration')
title('BB8')

for i=1:4
    v9{i} = [allBoutsBySolution9{i}.duration];
end
figure;
violin(v9); set(gca,'xtick',1:4,'xticklabel',{'H2O','.01M','.1M','1M'})
xlabel('Solution'); ylabel('Bout duration')
title('BB9')

for i=1:4
    allBoutsBySolution89{i} = [allBoutsBySolution8{i} allBoutsBySolution9{i}];
    v89{i} = [allBoutsBySolution89{i}.duration];
end
figure;
violin(v89); set(gca,'xtick',1:4,'xticklabel',{'H2O','.01M','.1M','1M'})
xlabel('Solution'); ylabel('Bout duration')
title('BB8 & BB9')

% Plot bout duration vs. mean or max of alternatives and the difference in
% palatability between current and alternative
figure;
subplot(1,2,1); scatter(meanAlts8,[allBouts8.duration],'.');
xlabel('Mean alternative palatability'); ylabel('Bout duration (s)')
title('BB8')
subplot(1,2,2); scatter(maxAlts8,[allBouts8.duration],'.');
xlabel('Max alternative palatability'); ylabel('Bout duration (s)')
title('BB8')
set(gcf,'Position',[100 100 1200 500])

figure;
subplot(1,2,1); scatter(meanAlts9,[allBouts9.duration],'.');
xlabel('Mean alternative palatability'); ylabel('Bout duration (s)')
title('BB9')
subplot(1,2,2); scatter(maxAlts9,[allBouts9.duration],'.');
xlabel('Max alternative palatability'); ylabel('Bout duration (s)')
title('BB9')
set(gcf,'Position',[100 100 1200 500])

figure;
subplot(2,1,1); scatter(altDifs8,[allBouts8.duration],'.')
xlabel('\Delta palatability'); ylabel('Bout duration (s)')
title('BB8')
subplot(2,1,2); scatter(altDifs9,[allBouts9.duration],'.')
xlabel('\Delta palatability'); ylabel('Bout duration (s)')
title('BB9')
set(gcf,'Position',[100 100 1200 500])

altDifs89 = [altDifs8 altDifs9];
durations89 = [[allBouts8.duration] [allBouts9.duration]];
durations89_1 = durations89;
figure;
p1=scatter(altDifs89,durations89,400,'.');
xlabel('\Delta palatability (current - alternative)'); ylabel('Bout duration (s)')
title('BB8 & BB9'); hold on;
u = unique(altDifs89);
for i=1:length(u)
    inds = find(altDifs89 == u(i));
    durs = durations89(inds);
    meanDur(i) = mean(durs);
    medianDur(i) = median(durs);
    muhat(i) = expfit(durs);
    p2=plot([(u(i) - .1) (u(i) + .1)],[meanDur(i) meanDur(i)],'k');
    p3 =plot([(u(i) - .1) (u(i) + .1)],[medianDur(i) medianDur(i)],'r');
    %hold on;
    %plot([(u(i) - .1) (u(i) + .1)],[muhat(i) muhat(i)],'r')
end
legend([p2 p3],{'Mean','Median'})
set(gcf,'Position',[100 100 1200 500])

% plot weights by date
figure;
plot(weights8); hold on; plot(weights9); set(gca,'xtick',1:length(dates),'xticklabel',truncDates,'xticklabelrotation',45);
ylabel('weight (g)'); xlabel('Date'); legend({'bb8','bb9'})

% get total # of licks and total bout duration by solution and use these to
% compute new palatibilities
nLicksBySolution8 = zeros(1,4);
totalDurationBySolution8 = zeros(1,4);
for i=1:length(allBoutsBySolution8)
    for j=1:length(allBoutsBySolution8{i})
        nLicksBySolution8(i) = nLicksBySolution8(i) + allBoutsBySolution8{i}(j).nlicks;
        totalDurationBySolution8(i) = totalDurationBySolution8(i) + allBoutsBySolution8{i}(j).duration;
    end
end
nLicksBySolution9 = zeros(1,4);
totalDurationBySolution9 = zeros(1,4);
for i=1:length(allBoutsBySolution9)
    for j=1:length(allBoutsBySolution9{i})
        nLicksBySolution9(i) = nLicksBySolution9(i) + allBoutsBySolution9{i}(j).nlicks;
        totalDurationBySolution9(i) = totalDurationBySolution9(i) + allBoutsBySolution9{i}(j).duration;
    end
end
% get relative palatibilities for each animal
relativeLicksBySolution8 = nLicksBySolution8;
relativeLicksBySolution8(2:end) = relativeLicksBySolution8(2:end)/2;
relativeLicksBySolution8 = relativeLicksBySolution8./relativeLicksBySolution8(1);
relativeLicksBySolution9 = nLicksBySolution9;
relativeLicksBySolution9(2:end) = relativeLicksBySolution9(2:end)/2;
relativeLicksBySolution9 = relativeLicksBySolution9./relativeLicksBySolution9(1);
relativeTotalDurationBySolution8 = totalDurationBySolution8;
relativeTotalDurationBySolution8(2:end) = relativeTotalDurationBySolution8(2:end)/2;
relativeTotalDurationBySolution8 = relativeTotalDurationBySolution8./relativeTotalDurationBySolution8(1);
relativeTotalDurationBySolution9 = totalDurationBySolution9;
relativeTotalDurationBySolution9(2:end) = relativeTotalDurationBySolution9(2:end)/2;
relativeTotalDurationBySolution9 = relativeTotalDurationBySolution9./relativeTotalDurationBySolution9(1);
% combine total data for both animals
nLicksBySolution89 = nLicksBySolution8 + nLicksBySolution9;
totalDurationBySolution89 = totalDurationBySolution8 + totalDurationBySolution9;
relativeLicksBySolution89 = nLicksBySolution89;
relativeLicksBySolution89(2:end) = relativeLicksBySolution89(2:end)/2;
relativeLicksBySolution89 = relativeLicksBySolution89./relativeLicksBySolution89(1);
relativeTotalDurationBySolution89 = totalDurationBySolution89;
relativeTotalDurationBySolution89(2:end) = relativeTotalDurationBySolution89(2:end)/2;
relativeTotalDurationBySolution89 = relativeTotalDurationBySolution89./relativeTotalDurationBySolution89(1);

[allBouts8,meanAlts8,maxAlts8,altDifs8,totalSampleDur8,totalSampleAltDifs8] =...
    combineBouts(basedir,dates,'bb8',relativeLicksBySolution8);
[allBouts9,meanAlts9,maxAlts9,altDifs9,totalSampleDur9,totalSampleAltDifs9] = combineBouts(basedir,dates,'bb9',relativeLicksBySolution9);
% remove bouts with single lick
nlicks8 = [allBouts8.nlicks];
moreThan1Lick8 = find(nlicks8 > 1);
nlicks9 = [allBouts9.nlicks];
moreThan1Lick9 = find(nlicks9 > 1);
allBouts8 = allBouts8(moreThan1Lick8);
meanAlts8 = meanAlts8(moreThan1Lick8);
maxAlts8 = maxAlts8(moreThan1Lick8);
altDifs8 = altDifs8(moreThan1Lick8);
allBouts9 = allBouts9(moreThan1Lick9);
meanAlts9 = meanAlts9(moreThan1Lick9);
maxAlts9 = maxAlts9(moreThan1Lick9);
altDifs9 = altDifs9(moreThan1Lick9);
for i=1:length(allBoutsBySolution8)
    nlicks = [allBoutsBySolution8{i}.nlicks];
    moreThan1 = find(nlicks > 1);
    allBoutsBySolution8{i} = allBoutsBySolution8{i}(moreThan1);
end
for i=1:length(allBoutsBySolution9)
    nlicks = [allBoutsBySolution9{i}.nlicks];
    moreThan1 = find(nlicks > 1);
    allBoutsBySolution9{i} = allBoutsBySolution9{i}(moreThan1);
end
altDifs89 = [altDifs8 altDifs9];
durations89 = [[allBouts8.duration] [allBouts9.duration]];
figure;
p1=scatter(altDifs89,durations89,400,'.');
xlabel('\Delta palatability (current - alternative)'); ylabel('Bout duration (s)')
title('BB8 & BB9 with palatabilities from data'); hold on;
u = unique(altDifs89);
for i=1:length(u)
    inds = find(altDifs89 == u(i));
    durs = durations89(inds);
    meanDur(i) = mean(durs);
    medianDur(i) = median(durs);
    muhat(i) = expfit(durs);
    p2=plot([(u(i) - .1) (u(i) + .1)],[meanDur(i) meanDur(i)],'k');
    p3 =plot([(u(i) - .1) (u(i) + .1)],[medianDur(i) medianDur(i)],'r');
    %hold on;
    %plot([(u(i) - .1) (u(i) + .1)],[muhat(i) muhat(i)],'r')
end
legend([p2 p3],{'Mean','Median'})
set(gcf,'Position',[100 100 1200 500])

allBouts89 = [allBouts8 allBouts9];
[allLicks] = getAllLicks(basedir,dates,rats);
boutOnsets = [allBouts89.onset];
lickOnsets = [allLicks.onset];
x=1:3600;
fBout = ksdensity(boutOnsets,x);
fLick = ksdensity(lickOnsets,x);
figure; plot(x,fBout(x)); hold on; plot(x,fLick(x))
xlabel('Time (s)'); ylabel('P(bout start OR lick)')
legend({'P(bout start)','P(lick)'})
set(gcf,'Position',[100 100 1200 500])

% plot total sample duration vs. palatability difference
totalSampleDur89 = [totalSampleDur8 totalSampleDur9];
totalSampleAltDifs89 = [totalSampleAltDifs8 totalSampleAltDifs9];
figure;
scatter(totalSampleAltDifs89,totalSampleDur89,400,'.');
xlabel('\Delta palatability (current - alternative)'); ylabel('Total time spent at option (s)')
title('BB8 & BB9'); set(gcf,'Position',[100 100 1000 500])
[rho,pval] = corr(totalSampleAltDifs89',totalSampleDur89');
text(-2,800,sprintf('\\rho = %f \np < 1e-4',[rho]),'FontSize',25,'FontWeight','bold')
f0 = fit(totalSampleAltDifs89',totalSampleDur89','poly1');
minx = min(totalSampleAltDifs89); maxx = max(totalSampleAltDifs89);
hold on; plot(minx:.1:maxx,f0(minx:.1:maxx))
ylim([0 1000])