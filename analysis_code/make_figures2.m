%% make_figures2
clear all; close all;
figFolder = '/home/ben/phd/behavior/figures/08_27_21/';
if (~exist(figFolder,'dir'))
    mkdir(figFolder)
end
do_save = 1; % Boolean for whether to save figures or not
animalInfo = load('analyzed_data/animalInfo.mat'); animalInfo=animalInfo.animalInfo;

if (~exist('animalObjs','var'))
    for i=1:length(animalInfo)
        animalObjs(i) = animalData(animalInfo(i).animal,animalInfo(i).dates,animalInfo(i).sex,'analyzed_data');
        disp(['done with ' animalInfo(i).animal])
    end
end
%B = load(['~/phd/talks/Pizza_Talks/nov_2019/B2.mat']); B = B.B;
%f=@(B,x) B(1).*exp(B(2).*x) + B(3);

% Create table with all bouts and applicable factors
[boutTable] = createBoutDataTable(animalObjs);

% get all bouts
allBouts = [];
for i=1:length(animalObjs)
    for j=1:length(animalObjs(i).dates)
        for k=1:length(animalObjs(i).bouts{j})
            if (~isempty(animalObjs(i).bouts{j}{k}))
                if (any([animalObjs(i).bouts{j}{k}.duration] < 0))
                    disp([animalObjs(i).name ' ' animalObjs(i).dates{j} ' has negative duration bout'])
                end
                allBouts = [allBouts animalObjs(i).bouts{j}{k}];
            end
        end
    end
end


%% FIGURE 1: cdf plot of lick times
cdfs = cell(1,length(animalObjs));
figure;
for i=1:length(animalObjs)
    cdfplot([animalObjs(i).linLicks.onset])
    hold on;
end
xlabel('Time (s)','fontsize',15,'fontweight','bold'); ylabel('P(t_{lick} < x)','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    figName = 'lick_times_cdf';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 2: cdf plot of bout onset times
figure;
for i=1:length(animalObjs)
    cdfplot([animalObjs(i).linBouts.onset])
    hold on;
end
xlabel('Time (s)','fontsize',15,'fontweight','bold'); ylabel('P(t_{bout} < x)','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    figName = 'bout_onset_times_cdf';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 3: bout lengths vs time
figure;
for i=1:length(animalObjs)
    scatter([animalObjs(i).linBouts.onset],[animalObjs(i).linBouts.duration],100,'k.')
    hold on;
end
xlabel('Time(s)','fontsize',15,'fontweight','bold'); ylabel('Bout duration (s)','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    figName = 'bout_lengths_vs_time';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
for i=1:length(animalObjs)
    scatter([animalObjs(i).linDifSolBouts.onset],[animalObjs(i).linDifSolBouts.duration],100,'k.')
    hold on;
end
title('Different solution days only')
xlabel('Time(s)','fontsize',15,'fontweight','bold'); ylabel('Bout duration (s)','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    figName = 'bout_lengths_vs_time_difSolDays';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 4: Bout duration as a function of total prior licks
difSolInds = boutTable.solutionNum ~= boutTable.alternative_solutionNum;
allBoutDurations = boutTable.duration(difSolInds);
allBoutNLicks = boutTable.nlicks(difSolInds);
prevLicks = boutTable.totalPrevLicks(difSolInds);
prevDurations = boutTable.totalPrevDuration(difSolInds);

[rho,pval] = corr(prevLicks,allBoutDurations);
exponentialFit = fit(prevLicks,allBoutDurations,'exp1');
posPrevLicks = find(prevLicks > 0);
noPrevLicks = find(prevLicks == 0);
powerFit = fit(prevLicks(posPrevLicks),allBoutDurations(posPrevLicks),'power1');
adjustedPrevLicks = prevLicks;
adjustedPrevLicks(noPrevLicks) = 1;
fullPowerFit = fit(adjustedPrevLicks,allBoutDurations,'power1');
lf = polyfit(prevLicks,allBoutDurations,1);
meanBoutDuration = mean(allBoutDurations);
SStot = sum((allBoutDurations - meanBoutDuration).^2);
SStotPosPrevLicks = sum((allBoutDurations(posPrevLicks) - meanBoutDuration).^2);
SSres_linear = sum((allBoutDurations - (lf(1)*prevLicks + lf(2))).^2);
SSres_exp = sum((allBoutDurations - exponentialFit.a*exp(exponentialFit.b*prevLicks)).^2);
SSres_power = sum((allBoutDurations(posPrevLicks) - powerFit.a*prevLicks(posPrevLicks).^powerFit.b).^2);
SSres_power_full = sum((allBoutDurations - fullPowerFit.a*adjustedPrevLicks.^fullPowerFit.b).^2);
rsquared_linear = 1 - SSres_linear/SStot;
rsquared_exp = 1 - SSres_exp/SStot;
rsquared_power = 1 - SSres_power/SStotPosPrevLicks;
rsquared_power_full = 1 - SSres_power_full/SStot;
xvals = 0:1:9000;
figure;
subplot(1,2,1)
scatter(prevLicks,allBoutDurations,100,'.')
hold on;
plot(xvals,lf(1)*xvals + lf(2))
plot(xvals,exponentialFit.a*exp(exponentialFit.b*xvals))
plot(xvals,powerFit.a*xvals.^powerFit.b)
plot(xvals,fullPowerFit.a*xvals.^fullPowerFit.b)
xlabel('Total prior licks','fontsize',15,'fontweight','bold')
ylabel('Bout duration','fontsize',15,'fontweight','bold')
text(3000,400,{['\rho = ' num2str(rho)],['p = ' num2str(pval)],['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],...
    ['R^2 (lin): ' num2str(rsquared_linear)],...
    ['R^2 (exp): ' num2str(rsquared_exp)],...
    ['R^2 (pow): ' num2str(rsquared_power)],...
    ['R^2 (pow full): ' num2str(rsquared_power_full)]})
legend({'Data','Linear Fit','Exponential Fit','Power law fit','full Power law fit'})

[rho2,pval2]=corr(prevDurations,allBoutDurations);
lf2 = polyfit(prevDurations,allBoutDurations,1);
subplot(1,2,2)
scatter(prevDurations,allBoutDurations,100,'.')
hold on;
xvals = 0:1:round(max(prevDurations));
plot(xvals,lf2(1)*xvals + lf(2))
xlabel('Total prior bout time','fontsize',15,'fontweight','bold')
ylabel('Bout duration','fontsize',15,'fontweight','bold')
text(max(xvals)/2,400,{['\rho = ' num2str(rho2)],['p = ' num2str(pval2)],['y = ' num2str(lf2(1)) 'x + ' num2str(lf2(2))]})
set(gcf,'Position',[10 10 1400 1200])
if (do_save)
    figName = 'bout_duration_vs_prior_licks';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 5: # of licks vs. bout duration
figure;
scatter(allBoutNLicks,allBoutDurations,100,'.')
xlabel('# of licks in bout','fontsize',15,'fontweight','bold')
ylabel('Bout duration','fontsize',15,'fontweight','bold')
if (do_save)
    figName = 'nlicks_vs_bout_duration';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 6: Average lick rate by solution
allNlicksBySolution = cell(1,4);
allDurationsBySolution = cell(1,4);
allLickRatesBySolution = cell(1,4);
for i=1:length(animalObjs)
    for j=1:length(animalObjs(i).boutsBySolution)
        curNLicks = [animalObjs(i).boutsBySolution{j}.nlicks];
        curDurations = [animalObjs(i).boutsBySolution{j}.duration];
        lickRates = curNLicks./curDurations;
        allNlicksBySolution{j} = [allNlicksBySolution{j} curNLicks];
        allDurationsBySolution{j} = [allDurationsBySolution{j} curDurations];
        allLickRatesBySolution{j} = [allLickRatesBySolution{j} lickRates];
    end
end
allLickRates = [];
solns = [];
for i=1:4
    allLickRates = [allLickRates allLickRatesBySolution{i}];
    solns = [solns ones(1,length(allLickRatesBySolution{i}))*i];
end
solnLabels = {'H2O','.01M','.1M','1M'};
figure;
notBoxPlot(allLickRates,solns)
set(gca,'xtick',1:4,'xticklabels',solnLabels)
ylabel('Lick rate','fontsize',15,'fontweight','bold')
if (do_save)
    figName = 'lickRateBySolution';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

% Do the same but for different solution days only (so no H2O)
allNlicksBySolutionDifSolDays = cell(1,3);
allDurationsBySolutionDifSolDays = cell(1,3);
allLickRatesBySolutionDifSolDays = cell(1,3);
for i=1:length(animalObjs)
    for j=2:length(animalObjs(i).boutsBySolutionDifSolDays)
        curNLicks = [animalObjs(i).boutsBySolutionDifSolDays{j}.nlicks];
        curDurations = [animalObjs(i).boutsBySolutionDifSolDays{j}.duration];
        lickRates = curNLicks./curDurations;
        allNlicksBySolution{j-1} = [allNlicksBySolution{j-1} curNLicks];
        allDurationsBySolution{j-1} = [allDurationsBySolution{j-1} curDurations];
        allLickRatesBySolution{j-1} = [allLickRatesBySolution{j-1} lickRates];
    end
end
allLickRatesDifSolDays = [];
solns = [];
for i=2:4
    allLickRatesDifSolDays = [allLickRatesDifSolDays allLickRatesBySolution{i}];
    solns = [solns ones(1,length(allLickRatesBySolution{i}))*(i-1)];
end
solnLabels = {'.01M','.1M','1M'};
figure;
notBoxPlot(allLickRatesDifSolDays,solns)
set(gca,'xtick',1:3,'xticklabels',solnLabels,'fontsize',15,'fontweight','bold')
ylabel('Lick rate','fontsize',15,'fontweight','bold')
title('Different solution days only','fontsize',15,'fontweight','bold')
if (do_save)
    figName = 'lickRateBySolutionDifSolDays';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 7
% bout duration as a function of previous bout duration
figure;
xs = [];
ys = [];
for i=1:length(animalObjs)
    for j=animalObjs(i).difSolutionDays
        xs = [xs [animalObjs(i).linBoutsByDay{j}(1:end-1).duration]];
        ys = [ys [animalObjs(i).linBoutsByDay{j}(2:end).duration]];
    end
end
[rho,pval] = corr(xs',ys');
lf = polyfit(xs,ys,1);
scatter(xs,ys,'k.'); xlabel('Bout(n) duration'); ylabel('Bout(n+1) duration')
x = 0:600;
hold on; plot(x,lf(1)*x + lf(2))
text(200,550,['\rho = ' num2str(rho)],'fontweight','bold','fontsize',20)
text(200,500,['p = ' num2str(pval)],'fontweight','bold','fontsize',20)
text(200,450,['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],'fontweight','bold','fontsize',20)
xlabel('Current bout duration (s)','fontsize',15,'fontweight','bold')
ylabel('Next bout duration (s)','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    figName = 'bout_duration_vs_prev_bout_duration';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 10
% relative palatability by amount consumed or total licks
figure;
for i=1:length(animalObjs)
    rpalLicks(:,i) = animalObjs(i).relativePalatabilitiesLicks;
    rpalConsumed(:,i) = animalObjs(i).relativePalatabilitiesConsumed;
end
xlabels = {'H2O','.01M','.1M','1M'};
subplot(1,2,1); plot(1:4,rpalConsumed); set(gca,'xtick',1:4,'xticklabels',xlabels); 
ylabel('Relative Palatabilitiy','fontsize',15,'fontweight','bold'); 
title('Relative Palatability based on consumption','fontsize',15,'fontweight','bold')
subplot(1,2,2); plot(1:4,rpalLicks); set(gca,'xtick',1:4,'xticklabels',xlabels); 
ylabel('RelativePalatability','fontsize',15,'fontweight','bold'); 
title('Relative Palatability based on lick totals','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1800 1000])
if (do_save)
    figName = 'relative_palatability_by_amount_consumed_or_licks';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 11: Relative palatability by solution and animal sex
maleRelativePalatabilities = cell(1,4);
femaleRelativePalatabilities = cell(1,4);
allRelativePalatabilities = [];
xinds = [];
for i=1:length(animalObjs)
    for j=1:4
        if (strcmp(animalObjs(i).sex,'M'))
            maleRelativePalatabilities{j} = [maleRelativePalatabilities{j} animalObjs(i).relativePalatabilitiesLicks(j)];
            xinds = [xinds 2*j - .5];
            allRelativePalatabilities = [allRelativePalatabilities animalObjs(i).relativePalatabilitiesLicks(j)];
        elseif (strcmp(animalObjs(i).sex,'F'))
            femaleRelativePalatabilities{j} = [femaleRelativePalatabilities{j} animalObjs(i).relativePalatabilitiesLicks(j)];
            xinds = [xinds 2*j + .5];
            allRelativePalatabilities = [allRelativePalatabilities animalObjs(i).relativePalatabilitiesLicks(j)];
        else
            error(['sex of ' animalObjs(i).name ' is not M/F'])
        end
    end
end
xlabels = {'H2O','.01M','.1M','1M'};
figure;
notBoxPlot(allRelativePalatabilities,xinds)
set(gca,'xtick',2:2:8,'xticklabels',xlabels,'fontsize',15,'fontweight','bold')
ylabel('Relative palatabilitiy','fontsize',15,'fontweight','bold')
if (do_save)
    figName = 'relativePalatabilities_by_solution_and_sex';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end
solnPairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4;];
for i=1:size(solnPairs,1)
    s1 = animalData.solnConverter(solnPairs(i,1));
    s2 = animalData.solnConverter(solnPairs(i,2));
    [h1,p1] = ttest2(maleRelativePalatabilities{solnPairs(i,1)},maleRelativePalatabilities{solnPairs(i,2)},'tail','left');
    [h2,p2] = ttest2(maleRelativePalatabilities{solnPairs(i,2)},maleRelativePalatabilities{solnPairs(i,1)},'tail','left');
    if (p1 < .05)
        disp(['For males, ' s1 ' is significantly less palatable than ' s2 ' | p = ' num2str(p1)])
    elseif (p2 < .05)
        disp(['For males, ' s2 ' is significantly less palatable than ' s1 ' | p = ' num2str(p2)])
    end
end
for i=1:size(solnPairs,1)
    s1 = animalData.solnConverter(solnPairs(i,1));
    s2 = animalData.solnConverter(solnPairs(i,2));
    [h1,p1] = ttest2(femaleRelativePalatabilities{solnPairs(i,1)},femaleRelativePalatabilities{solnPairs(i,2)},'tail','left');
    [h2,p2] = ttest2(femaleRelativePalatabilities{solnPairs(i,2)},femaleRelativePalatabilities{solnPairs(i,1)},'tail','left');
    if (p1 < .05)
        disp(['For females, ' s1 ' is significantly less palatable than ' s2 ' | p = ' num2str(p1)])
    elseif (p2 < .05)
        disp(['For females, ' s2 ' is significantly less palatable than ' s1 ' | p = ' num2str(p2)])
    end
end

%% Figure 12: Bout duration distributions for each solution
allBoutsBySolution = cell(1,4);
allBoutsBySolutionDifSolDays = cell(1,4);
maxDurs = zeros(1,4);
for i=1:length(animalObjs)
    for j=1:4
        allBoutsBySolution{j} = [allBoutsBySolution{j} [animalObjs(i).boutsBySolution{j}.duration]];
        maxDurs(j) = max(maxDurs(j),max(allBoutsBySolution{j}));
    end
end
maxDursDifSolDays = zeros(1,3);
for i=1:length(animalObjs)
    count = 0;
    for j=2:4
        count = count + 1;
        allBoutsBySolutionDifSolDays{count} = [allBoutsBySolutionDifSolDays{count} [animalObjs(i).boutsBySolutionDifSolDays{j}.duration]];
        maxDursDifSolDays(count) = max(maxDursDifSolDays(count),max(allBoutsBySolutionDifSolDays{count}));
    end
end
maxDurs = ceil(maxDurs+1);
maxDursDifSolDays = ceil(maxDursDifSolDays+1);
figure;
subplot(2,2,1); histogram(allBoutsBySolution{1},'Normalization','pdf','binwidth',1); 
xlabel('Bout duration (s)','fontsize',15,'fontweight','bold'); title('H2O','fontsize',15,'fontweight','bold'); xlim([0 maxDurs(1)]); ylim([0 .2])
subplot(2,2,2); histogram(allBoutsBySolution{2},'Normalization','pdf','binwidth',1); 
xlabel('Bout duration (s)','fontsize',15,'fontweight','bold'); title('.01M NaCl','fontsize',15,'fontweight','bold'); xlim([0 maxDurs(2)]); ylim([0 .2])
subplot(2,2,3); histogram(allBoutsBySolution{3},'Normalization','pdf','binwidth',1); 
xlabel('Bout duration (s)','fontsize',15,'fontweight','bold'); title('.1M NaCl','fontsize',15,'fontweight','bold'); xlim([0 maxDurs(3)]); ylim([0 .2])
subplot(2,2,4); histogram(allBoutsBySolution{4},'Normalization','pdf','binwidth',1); 
xlabel('Bout duration (s)','fontsize',15,'fontweight','bold'); title('1M NaCl','fontsize',15,'fontweight','bold'); xlim([0 maxDurs(4)]); ylim([0 .2])
set(gcf,'Position',[10 10 1600 1600])
if (do_save)
    figName = 'bout_duration_distribution_by_solution';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
subplot(2,2,1); histogram(allBoutsBySolutionDifSolDays{1},'Normalization','pdf','binwidth',1); 
xlabel('Bout duration (s)','fontsize',15,'fontweight','bold'); title('.01M','fontsize',15,'fontweight','bold'); xlim([0 maxDursDifSolDays(1)]); ylim([0 .2])
subplot(2,2,2); histogram(allBoutsBySolutionDifSolDays{2},'Normalization','pdf','binwidth',1);
xlabel('Bout duration (s)','fontsize',15,'fontweight','bold'); title('.1M NaCl','fontsize',15,'fontweight','bold'); xlim([0 maxDursDifSolDays(2)]); ylim([0 .2])
subplot(2,2,3); histogram(allBoutsBySolutionDifSolDays{3},'Normalization','pdf','binwidth',1); 
xlabel('Bout duration (s)','fontsize',15,'fontweight','bold'); title('1M NaCl','fontsize',15,'fontweight','bold'); xlim([0 maxDursDifSolDays(3)]); ylim([0 .2])
set(gcf,'Position',[10 10 1600 1600])
if (do_save)
    figName = 'bout_duration_distributions_by_solution_difSolDays';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 13: Exponential fits to bout distributions for each solution
figure;
for i=1:4
    muhat = expfit(allBoutsBySolution{i});
    x = 1:150;
    plot(x,(1/muhat)*exp(-(1/muhat)*x))
    hold on;
end
xlabel('Bout Duration (s)','fontsize',15,'fontweight','bold'); 
ylabel('Probability density','fontsize',15,'fontweight','bold')
legend({'H2O','.01M','.1M','1M'},'fontsize',15,'fontweight','bold'); 
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    figName = 'bout_duration_exp_fits';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
for i=1:3
    muhat = expfit(allBoutsBySolutionDifSolDays{i});
    x = 1:150;
    plot(x,(1/muhat)*exp(-(1/muhat)*x))
    hold on;
end
xlabel('Bout Duration (s)','fontsize',15,'fontweight','bold'); 
ylabel('Probability Density','fontsize',15,'fontweight','bold')
legend({'.01M','.1M','1M'},'fontsize',15,'fontweight','bold'); 
set(gcf,'Position',[10 10 1200 1000])
title('Different solution days only')
if (do_save)
    figName = 'bout_duration_exp_fits_difSolDays';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

solnNames = {'.01M','.1M','1M'};
figure;
for i=1:3
    subplot(2,2,i)
    muhat = expfit(allBoutsBySolutionDifSolDays{i});
    x=1:150;
    plot(x,(1/muhat)*exp(-(1/muhat)*x))
    hold on;
    histogram(allBoutsBySolutionDifSolDays{i},'binwidth',2,'normalization','pdf')
    xlabel('Bout Duration (s)','fontsize',15,'fontweight','bold')
    ylabel('Probability Density','fontsize',15,'fontweight','bold')
    title(solnNames{i},'fontsize',15,'fontweight','bold')
    xlim([0 80])
end
set(gcf,'Position',[10 10 1400 1000])
if (do_save)
    figName = 'bout_duration_exp_fits_difSolDays_w_data';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

muhats = zeros(length(animalObjs),4);
muhatBoxPlotY = [];
muhatBoxPlotX = [];
for i=1:length(animalObjs)
    ratInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    for j=1:4
        solnInds = find(boutTable.solutionNum == j);
        muhats(i,j) = expfit(boutTable.duration(intersect(ratInds,solnInds)));
        muhatBoxPlotY = [muhatBoxPlotY muhats(i,j)];
        muhatBoxPlotX = [muhatBoxPlotX j];
    end
end
figure;
notBoxPlot(muhatBoxPlotY,muhatBoxPlotX)
set(gca,'xtick',1:4,'xticklabels',{'dH2O','.01M','.1M','1M'},'fontsize',15,'fontweight','bold')
ylabel('Mean bout duration','fontsize',15,'fontweight','bold')
title('All days','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1400 1000])
if (do_save)
    figName = 'bout_duration_exp_fits_boxPlot_allDays';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

muhatBoxPlotX = [];
muhatBoxPlotY = [];
muhatsSameSolDays = zeros(length(animalObjs),4);
for i=1:length(animalObjs)
    ratInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    for j=1:4
        sameSolInds = find(boutTable.solutionNum == j & boutTable.alternative_solutionNum == j);
        curInds = intersect(ratInds,sameSolInds);
        muhat = expfit(boutTable.duration(curInds));
        muhatsSameSolDays(i,j) = muhat;
        muhatBoxPlotY = [muhatBoxPlotY muhat];
        muhatBoxPlotX = [muhatBoxPlotX j];
    end
end
figure;
notBoxPlot(muhatBoxPlotY,muhatBoxPlotX)
set(gca,'xtick',1:4,'xticklabels',{'dH2O','.01M','.1M','1M'},'fontsize',15,'fontweight','bold')
ylabel('Mean bout duration','fontsize',15,'fontweight','bold')
title('Same solution days only','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1400 1000])
if (do_save)
    figName = 'boutDuration_exp_fits_boxPlot_sameSolDays';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 14: Bout durations as a function of solution palatability
figure;
bd = [];
relPals = [];
for i=1:length(animalObjs)
    for j=1:4
        b = [animalObjs(i).boutsBySolution{j}.duration];
        bd = [bd b];
        if (any(isnan(animalObjs(i).relativePalatabilitiesLicks(j))))
            disp([animalObjs(i).name ' ' num2str(j) ' has a nan palatability'])
        end
        relPals = [relPals ones(1,length(b))*animalObjs(i).relativePalatabilitiesLicks(j)];
        %scatter(ones(1,length(b))*animalObjs(i).relativePalatabilitiesConsumed(j),b,'k.')
    end
end
[rho,pval] = corr(relPals',bd');
lf = polyfit(relPals,bd,1);
x = min(relPals):.001:max(relPals);
scatter(relPals,bd,'k.'); hold on; plot(x,x*lf(1)+lf(2))
xlabel('Palatability','fontsize',15,'fontweight','bold'); 
ylabel('Bout Duration (s)','fontsize',15,'fontweight','bold')
text(.2,550,['\rho = ' num2str(rho)],'fontweight','bold','fontsize',20); 
text(.2, 500,['p = ' num2str(pval)],'fontweight','bold','fontsize',20);
text(.2,450,['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    figName = 'bout_duration_vs_current_palatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 15: Bout durations as a function of solution palatability, different solution days only
figure;
bd = [];
relPals = [];
for i=1:length(animalObjs)
    for j=1:4
        if (~isempty(animalObjs(i).boutsBySolutionDifSolDays{j}))
            b = [animalObjs(i).boutsBySolutionDifSolDays{j}.duration];
            bd = [bd b];
            if (any(isnan(animalObjs(i).relativePalatabilitiesLicks(j))))
                disp([animalObjs(i).name ' ' num2str(j) ' has a nan palatability'])
            end
            relPals = [relPals ones(1,length(b))*animalObjs(i).relativePalatabilitiesLicks(j)];
        end
    end
end
[rho,pval] = corr(relPals',bd');
lf = polyfit(relPals,bd,1);
x = min(relPals):.001:max(relPals);
scatter(relPals,bd,'k.'); hold on; plot(x,x*lf(1)+lf(2))
xlabel('Palatability','fontsize',15,'fontweight','bold'); 
ylabel('Bout Duration (s)','fontsize',15,'fontweight','bold'); 
title('Different solution days only','fontsize',15,'fontweight','bold')
text(.2,550,['\rho = ' num2str(rho)],'fontweight','bold','fontsize',20); 
text(.2, 500,['p = ' num2str(pval)],'fontweight','bold','fontsize',20);
text(.2,450,['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    figName = 'bout_duration_vs_current_palatability_difSolDays';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 15: Bout durations as a function of alternative palatability
difSolInds = boutTable.solutionNum ~= boutTable.alternative_solutionNum;
meanAltPals = boutTable.alternative_palatability(difSolInds);
bd = boutTable.duration(difSolInds);

[rho,pval] = corr(meanAltPals,bd);
lf = polyfit(meanAltPals,bd,1);
x = min(meanAltPals):.001:max(meanAltPals);
figure;
scatter(meanAltPals,bd,'k.'); hold on; plot(x,x*lf(1)+lf(2));
xlabel('Alternative palatability','fontsize',15,'fontweight','bold'); 
ylabel('Bout duration (s)','fontsize',15,'fontweight','bold')
text(.5,550,['\rho = ' num2str(rho)],'fontsize',20,'fontweight','bold')
text(.5,500,['p = ' num2str(pval)],'fontsize',20,'fontweight','bold')
text(.5,450,['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],'fontsize',20,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    figName = 'bout_duration_vs_alternative_palatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 16: Bout durations as a function of relative palatability difference
difSolInds = boutTable.solutionNum ~= boutTable.alternative_solutionNum;
allPalDifs = boutTable.palatability(difSolInds) - boutTable.alternative_palatability(difSolInds);
allBoutDurs = boutTable.duration(difSolInds);
uPals = unique(allPalDifs);

figure;
for i=1:length(uPals)
    inds = find(allPalDifs == uPals(i));
    durs = allBoutDurs(inds);
    muhats2(i) = expfit(durs);
    medianBoutDurs(i) = median(durs);
    palDifs(i) = allPalDifs(inds(1));
    scatter(allPalDifs(inds),allBoutDurs(inds),'k.'); hold on;
    plot([palDifs(i)-.1 palDifs(i)+.1],[muhats2(i) muhats2(i)],'r')
end
xlabel('Palatability difference','fontsize',15,'fontweight','bold'); 
ylabel('Bout duration (s)','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1400 1000])
if (do_save)
    figName = 'bout_duration_vs_palatability_difference';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 17: Bout duration vs palatability difference with linear fit and correlation
figure;
[rho,pval] = corr(allPalDifs,allBoutDurs);
scatter(allPalDifs,allBoutDurs,'k.')
lf = polyfit(allPalDifs,allBoutDurs,1);
x = min(allPalDifs):.01:max(allPalDifs);
hold on; plot(x,lf(1)*x + lf(2))
xlabel('Palatability difference','fontsize',15,'fontweight','bold');
ylabel('Bout duration (s)','fontsize',15,'fontweight','bold')
text(-1.8,460,['\rho = ' num2str(rho)],'fontweight','bold','fontsize',20)
text(-1.8,430,['p = ' num2str(pval)],'fontweight','bold','fontsize',20)
text(-1.8,400,['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],'fontweight','bold','fontsize',20)
ylim([0 max(allBoutDurs)+10])
set(gcf,'Position',[10 10 1400 1000])
if (do_save)
    figName = 'bout_duration_vs_palatability_difference_linfit_corr';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 18: Bout duration Exponential means and medians as a function of palatability difference
figure;
subplot(1,2,1)
linfit = polyfit(palDifs,muhats2,1);
linfit2 = polyfit(allPalDifs,allBoutDurs,1);
[rho,pval] = corr(palDifs',muhats2');
x = min(palDifs):.01:max(palDifs);
scatter(palDifs,muhats2,'k.'); hold on; h1=plot(x,x*linfit(1) + linfit(2)); 
%h2=plot(x,x*linfit2(1) + linfit2(2));
ylim([0 70])
xlabel('Palatability difference','fontsize',15,'fontweight','bold'); 
ylabel('Exponential mean','fontsize',15,'fontweight','bold'); 
legend([h1],{'Fit to exponential means'},'fontsize',15,'fontweight','bold')
text(-1.8,60,['\rho = ' num2str(rho)],'fontweight','bold','fontsize',20); 
text(-1.8,55,['p = ' num2str(pval)],'fontweight','bold','fontsize',20)
text(-1.8,50,['y = ' num2str(linfit(1)) 'x + ' num2str(linfit(2))],'fontweight','bold','fontsize',20)

subplot(1,2,2)
linfit3 = polyfit(palDifs,medianBoutDurs,1);
[rho2,pval2] = corr(palDifs',medianBoutDurs');
scatter(palDifs,medianBoutDurs,'k.'); hold on; h3 = plot(x,x*linfit3(1) + linfit3(2));
ylim([0 70])
xlabel('Palatability difference','fontsize',15,'fontweight','bold');
ylabel('Median bout length','fontsize',15,'fontweight','bold');
legend(h3,{'Fit to median bout lengths'},'fontsize',15,'fontweight','bold')
text(-1.8,60,['\rho = ' num2str(rho2)],'fontweight','bold','fontsize',20); 
text(-1.8,55,['p = ' num2str(pval2)],'fontweight','bold','fontsize',20)
text(-1.8,50,['y = ' num2str(linfit3(1)) 'x + ' num2str(linfit3(2))],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'bout_duration_vs_palatability_difference_exp_fit_lin_fit';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 19: Bout lengths as a linear model of current and alternative palatability
Y = cell(1,length(animalObjs));
X = cell(1,length(animalObjs));
for i=1:length(animalObjs)
    count = 0;
    for j=animalObjs(i).difSolutionDays
        for k=1:length(animalObjs(i).bouts{j})
            for l=1:length(animalObjs(i).bouts{j}{k})
                count = count+1;
                Y{i}(count) = animalObjs(i).bouts{j}{k}(l).duration;
                X{i}(count,1) = animalObjs(i).Palatabilities{j}{k}(l);
                X{i}(count,2) = animalObjs(i).meanAlternativePalatabilities{j}{k}(l);
                X{i}(count,3) = 1;
            end
        end
    end
    [b,bint,r,rint,stats] = regress(Y{i}',X{i});
    multiLinearBoutDurFits(i).b = b;
    multiLinearBoutDurFits(i).normalizedB = b./mean(Y{i});
    multiLinearBoutDurFits(i).bint = bint;
    multiLinearBoutDurFits(i).r = r;
    multiLinearBoutDurFits(i).rint = rint;
    multiLinearBoutDurFits(i).stats = stats;
    %{
    figure;
    scatter3(X{i}(:,1),X{i}(:,2),Y{i},50,'.')
    hold on;
    %}
end
z = [multiLinearBoutDurFits.b];
normalizedZ = [multiLinearBoutDurFits.normalizedB];
[h1,p1] = ttest(z(1,:));
[h2,p2] = ttest(z(2,:));
[h3,p3] = ttest(z(3,:));
[pr1,hr1] = signrank(z(1,:));
[pr2,hr2] = signrank(z(2,:));
[pr3,hr3] = signrank(z(3,:));
disp(['Result of t-test for Palatability coefficient: H = ' num2str(h1) ', p = ' num2str(p1) ', H = ' num2str(hr1) ', P = ' num2str(pr1)])
disp(['Result of t-test for Alternative Palatability coefficient: H = ' num2str(h2) ', p = ' num2str(p2) ', H = ' num2str(hr2) ', P = ' num2str(pr2)])
disp(['Result of t-test for intercept: H = ' num2str(h3) ', p = ' num2str(p3) ', H = ' num2str(hr3) ', P = ' num2str(pr3)])

cmap = jet;
figure;
for i=1:length(animalObjs)
    Zs{i} = multiLinearBoutDurFits(i).b(2)*(0:.01:3)' + multiLinearBoutDurFits(i).b(1)*(0:.01:3) + multiLinearBoutDurFits(i).b(3);
    surf(0:.01:3,0:.01:3,Zs{i},'FaceColor',cmap(i*floor(size(cmap,1)/length(animalObjs)),:),'FaceAlpha',.3,'EdgeColor',cmap(i*floor(size(cmap,1)/length(animalObjs)),:));
    hold on;
end
xlabel('Current Palatability','fontsize',15,'fontweight','bold');
ylabel('Alternative palatability','fontsize',15,'fontweight','bold');
zlabel('Predicted Bout duration','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1400 1200])
if (do_save)
    figName = 'palatability_alternativePalatability_boutDuration_planes';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
notBoxPlot([z(1,:) z(2,:)],[ones(1,length(z(1,:))) ones(1,length(z(2,:)))*2])
ylabel('Regression coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Palatability','Alternative palatability'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'all_bouts_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
notBoxPlot([normalizedZ(1,:) normalizedZ(2,:)],[ones(1,length(normalizedZ(1,:))) ones(1,length(normalizedZ(2,:)))*2])
ylabel('Normalized regression coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Palatability','Alternative palatability'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'all_bouts_normalized_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Example figure demonstrating multilinear regression
animal2use = 1;
vals = 0:.01:2.5;
figure;
subplot(1,2,1)
plot(vals,z(1,animal2use)*vals + z(2,animal2use)*.5 + z(3,animal2use))
hold on;
plot(vals,z(1,animal2use)*vals + z(2,animal2use)*2 + z(3,animal2use))
legend({'Low alternative palatability','High alternative palatability'},'fontsize',15,'fontweight','bold')
xlabel('Current palatability','fontsize',15,'fontweight','bold')
ylabel('Bout duration (s)','fontsize',15,'fontweight','bold')
ylim([0 25])

subplot(1,2,2)
plot(vals,z(1,animal2use)*.5 + z(2,animal2use)*vals + z(3,animal2use))
hold on;
plot(vals,z(1,animal2use)*2 + z(2,animal2use)*vals + z(3,animal2use))
legend({'Low current palatability','High current palatability'},'fontsize',15,'fontweight','bold')
xlabel('Alternative palatability','fontsize',15,'fontweight','bold')
ylabel('Bout duration (s)','fontsize',15,'fontweight','bold')
ylim([0 25])
set(gcf,'Position',[10 10 1400 600])

if (do_save)
    figName = 'Multilinear_regression_example';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end
%% Bout duration ANOVA
difSolInds = boutTable.solutionNum ~= boutTable.alternative_solutionNum;
y = boutTable.duration(difSolInds);
factors = {boutTable.rat(difSolInds), boutTable.palatability(difSolInds),...
           boutTable.alternative_palatability(difSolInds),boutTable.totalPrevLicks(difSolInds),...
           boutTable.back1_duration(difSolInds), boutTable.back2_duration(difSolInds), ...
           boutTable.back3_duration(difSolInds)};

varnames = {'rat','palatability','alternative_palatability','totalPrevLicks',...
            'back1_duration','back2_duration','back3_duration'};
isContinuous = [2,3,4,5,6,7];
isRandom = [1];
%isContinuous = [2 3];
[p,tbl,stats,terms] = anovan(y,factors,'model','interaction','continuous',...
                             isContinuous,'random',isRandom,'varnames',varnames);
                         
%% Bout duration ANOVA (only palatability factors)
difSolInds = boutTable.solutionNum ~= boutTable.alternative_solutionNum;
y = boutTable.duration(difSolInds);
factors = {boutTable.palatability(difSolInds), boutTable.alternative_palatability(difSolInds)};
varnames = {'palatability','alternative_palatability'};
isContinuous = [1 2];
[p,tbl,stats,terms] = anovan(y,factors,'model','interaction','continuous',...
                             isContinuous,'varnames',varnames);

[p,tbl,stats,terms] = anovan(y,factors,'model','linear','continuous',...
                             isContinuous,'varnames',varnames);
%% Bout duration Scheirer-Ray-Hare test
[p,H,p2] = my_scheirer_ray_hare(y,factors,[1 2])
                         
%% Bout duration ANOVA (only palatability factors + rat random factor)
difSolInds = boutTable.solutionNum ~= boutTable.alternative_solutionNum;
y = boutTable.duration(difSolInds);
factors = {boutTable.palatability(difSolInds), boutTable.alternative_palatability(difSolInds), boutTable.rat(difSolInds)};
varnames = {'palatability','alternative_palatability','rat'};
isContinuous = [1 2];
isRandom = 3;
[p,tbl,stats,terms] = anovan(y,factors,'model','interaction','continuous',...
                             isContinuous,'random',isRandom,'varnames',varnames);
                         
[p,tbl,stats,terms] = anovan(y,factors,'model','linear','continuous',...
                             isContinuous,'random',isRandom,'varnames',varnames);
%% Bout duration Scheirer-Ray-Hare test with rat as a random factor

%% Linear model of bout durations: All animals combined
difSolInds = boutTable.solutionNum ~= boutTable.alternative_solutionNum;
Y = boutTable.duration(difSolInds);
X = [boutTable.palatability(difSolInds) boutTable.alternative_palatability(difSolInds)...
     boutTable.totalPrevLicks(difSolInds) boutTable.back1_duration(difSolInds)...
     boutTable.back2_duration(difSolInds) boutTable.back3_duration(difSolInds) ...
     ones(sum(difSolInds),1)];
[b,bint,r,rint,stats] = regress(Y,X);

%% Figure 27: Animal specific linear models of bout duration following a stay or switch decision
clear switchDurs switchAltPals switchPals stayDurs stayAltPals stayPals
switchDurs = [];
switchAltPals = [];
switchPals = [];
stayDurs = [];
stayAltPals = [];
stayPals = [];
nStay = 0;
nFirstStay = 0;
firstStayDurs = [];
firstStayAltPals = [];
firstStayPals = [];
count = 0;
switchCount = 0;
stayCount = 0;
%for i=setdiff(1:length(animalObjs),5) % exclude bb12
for i=1:length(animalObjs)
    count = count+1;
    switchX = [];
    switchY = [];
    stayX = [];
    stayY = [];
    firstStayX = [];
    firstStayY = [];
    firstStay = true;
    for j=animalObjs(i).difSolutionDays
        for k=2:length(animalObjs(i).linBoutsByDay{j})
            curPal = animalObjs(i).relativePalatabilitiesLicks(animalObjs(i).solnConverter(animalObjs(i).linBoutsByDay{j}(k).solution));
            curDur = animalObjs(i).linBoutsByDay{j}(k).duration;
            if (animalObjs(i).linBoutsByDay{j}(k).channel ~= animalObjs(i).linBoutsByDay{j}(k-1).channel)
                switchCount = switchCount+1;
                switchRats{switchCount} = animalObjs(i).name;
                curAltPal = animalObjs(i).relativePalatabilitiesLicks(animalObjs(i).solnConverter(animalObjs(i).linBoutsByDay{j}(k-1).solution));
                switchDurs = [switchDurs curDur];
                switchAltPals = [switchAltPals curAltPal];
                switchPals = [switchPals curPal];
                switchX = [switchX; [curPal curAltPal]];
                switchY = [switchY curDur];
                firstStay = true;
            else   
                stayCount = stayCount + 1;
                stayRats{stayCount} = animalObjs(i).name;
                curSolnInd = animalObjs(i).solnConverter(animalObjs(i).linBoutsByDay{j}(k).solution);
                curDaySolns = animalObjs(i).solutions(j,:);
                altSolnInd = setdiff(curDaySolns,curSolnInd);
                curAltPal = animalObjs(i).relativePalatabilitiesLicks(altSolnInd);
                stayDurs = [stayDurs curDur];
                stayAltPals = [stayAltPals curAltPal];
                stayPals = [stayPals curPal];
                stayX = [stayX; [curPal curAltPal]];
                stayY = [stayY curDur];
                if (firstStay)
                    firstStayDurs = [firstStayDurs curDur];
                    firstStayAltPals = [firstStayAltPals curAltPal];
                    firstStayPals = [firstStayPals curPal];
                    firstStayX = [firstStayX; [curPal curAltPal]];
                    firstStayY = [firstStayY curDur];
                    nFirstStay = nFirstStay + 1;
                    firstStay = false;
                end
                nStay = nStay + 1;
            end
        end
    end
    switchX = [switchX ones(size(switchX,1),1)];
    [b,bint,r,rint,stats] = regress(switchY',switchX);
    switchBs(count,:) = b;
    switchNormalizedBs(count,:) = b./mean(switchY);
    switchBINTs{count} = bint;
    switchRs{count} = r;
    switchRINTs{count} = rint;
    switchSTATs{count} = stats;
    
    stayX = [stayX ones(size(stayX,1),1)];
    [b,bint,r,rint,stats] = regress(stayY',stayX);
    stayBs(count,:) = b;
    stayNormalizedBs(count,:) = b./mean(stayY);
    stayBINTs{count} = bint;
    stayRs{count} = r;
    stayRINTs{count} = rint;
    staySTATs{count} = stats;
    
    firstStayX = [firstStayX ones(size(firstStayX,1),1)];
    [b,bint,r,rint,stats] = regress(firstStayY',firstStayX);
    firstStayBs(count,:) = b;
    firstStayNormalizedBs(count,:) = b./mean(stayY);
    firstStayBINTs{count} = bint;
    firstStayRs{count} = r;
    firstStayRINTs{count} = rint;
    firstStaySTATs{count} = stats;
end
% Plot bout duration vs alternative palatability after a stay decision
[rho,pval] = corr(stayDurs',stayAltPals');
lf = polyfit(stayAltPals,stayDurs,1);
figure;
scatter(stayAltPals,stayDurs,100,'k.')
x=0:.01:3; hold on; plot(x,lf(1)*x + lf(2))
xlabel('Alternative palatability','fontsize',15,'fontweight','bold')
ylabel('Bout duration','fontsize',15,'fontweight','bold')
title('After stay decision','fontsize',15,'fontweight','bold')
text(1.2,200,['\rho = ' num2str(rho)],'fontweight','bold','fontsize',20)
text(1.2,220,['p = ' num2str(pval)],'fontweight','bold','fontsize',20)
text(1.2,240,['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'bout_duration_after_stay_vs_alternative_palatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

% Plot bout duration vs palatability after a stay decision
[rho,pval] = corr(stayDurs',stayPals');
lf = polyfit(stayPals,stayDurs,1);
figure;
scatter(stayPals,stayDurs,100,'k.')
x=0:.01:3; hold on; plot(x,lf(1)*x + lf(2))
xlabel('Palatability','fontsize',15,'fontweight','bold')
ylabel('Bout duration','fontsize',15,'fontweight','bold')
title('After stay decision','fontsize',15,'fontweight','bold')
text(1.2,200,['\rho = ' num2str(rho)],'fontweight','bold','fontsize',20)
text(1.2,220,['p = ' num2str(pval)],'fontweight','bold','fontsize',20)
text(1.2,240,['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'bout_duration_after_stay_vs_palatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

% Plot bout duration vs alternative palatability after a switch decision
[rho,pval] = corr(switchDurs',switchAltPals');
lf = polyfit(switchAltPals,switchDurs,1);
figure;
scatter(switchAltPals,switchDurs,100,'k.')
x=0:.01:3; hold on; plot(x,lf(1)*x + lf(2))
xlabel('Alternative palatability','fontsize',15,'fontweight','bold')
ylabel('Bout duration','fontsize',15,'fontweight','bold')
title('After switch decision','fontsize',15,'fontweight','bold')
text(1.2,200,['\rho = ' num2str(rho)],'fontweight','bold','fontsize',20)
text(1.2,220,['p = ' num2str(pval)],'fontweight','bold','fontsize',20)
text(1.2,240,['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'bout_duration_after_switch_vs_alternative_palatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

% Plot bout duration vs palatability after a switch decision
[rho,pval] = corr(switchDurs',switchPals');
lf = polyfit(switchPals,switchDurs,1);
figure;
scatter(switchPals,switchDurs,100,'k.')
x=0:.01:3; hold on; plot(x,lf(1)*x + lf(2))
xlabel('Palatability','fontsize',15,'fontweight','bold')
ylabel('Bout duration','fontsize',15,'fontweight','bold')
title('After switch decision','fontsize',15,'fontweight','bold')
text(1.2,200,['\rho = ' num2str(rho)],'fontweight','bold','fontsize',20)
text(1.2,220,['p = ' num2str(pval)],'fontweight','bold','fontsize',20)
text(1.2,240,['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'bout_duration_after_switch_vs_palatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% ANOVA of bout duration (with palatabilities as factors) after stay decision
factors = {stayPals',stayAltPals'};
varnames = {'(stay) palatability','(stay) alternative_palatability'};
isContinuous = [1 2];
[p,tbl,stats,terms] = anovan(stayDurs,factors,'model','interaction','continuous',...
                             isContinuous,'varnames',varnames);

[p,tbl,stats,terms] = anovan(stayDurs,factors,'model','linear','continuous',...
                             isContinuous,'varnames',varnames);
%% ANOVA of bout duration (with palatabilities as factors) after switch decision 
factors = {switchPals',switchAltPals'};
varnames = {'(switch) palatability','(switch) alternative_palatability'};
isContinuous = [1 2];
[p,tbl,stats,terms] = anovan(switchDurs,factors,'model','interaction','continuous',...
                             isContinuous,'varnames',varnames);

[p,tbl,stats,terms] = anovan(switchDurs,factors,'model','linear','continuous',...
                             isContinuous,'varnames',varnames);
%% ANOVA of bout duration (with palatabilities and rat as factors) after stay decision
factors = {stayPals',stayAltPals', stayRats};
varnames = {'(stay) palatability','(stay) alternative_palatability', 'rat'};
isContinuous = [1 2];
isRandom = [3];
[p,tbl,stats,terms] = anovan(stayDurs,factors,'model','interaction','continuous',...
                             isContinuous,'random',isRandom,'varnames',varnames);

%% ANOVA of bout duration (with palatabilities and rat as factors) after switch decision
factors = {switchPals',switchAltPals', switchRats};
varnames = {'(switch) palatability','(switch) alternative_palatability', 'rat'};
isContinuous = [1 2];
isRandom = [3];
[p,tbl,stats,terms] = anovan(switchDurs,factors,'model','interaction','continuous',...
                             isContinuous,'random',isRandom,'varnames',varnames);
%% Comparison of stay/switch bout duration linear model coefficients
[p,h] = signrank(stayBs(:,1));
if (p < .05)
    disp(['Significant effect of palatability on stay bout durations p = ' num2str(p)])
else
    disp(['No significant effect of palatability on stay bout durations p = ' num2str(p)])
end
[p,h] = signrank(stayBs(:,2));
if (p < .05)
    disp(['Significant effect of alternative palatability on stay bout durations p = ' num2str(p)])
else
    disp(['No significant effect of alternative palatability on stay bout durations p = ' num2str(p)])
end
[p,h] = signrank(switchBs(:,1));
if (p < .05)
    disp(['Significant effect of palatability on switch bout durations p = ' num2str(p)])
else
    disp(['No significant effect of palatability on switch bout durations p = ' num2str(p)])
end
[p,h] = signrank(switchBs(:,2));
if (p < .05)
    disp(['Significant effect of alternative palatability on switch bout durations p = ' num2str(p)])
else
    disp(['No significant effect of alternative palatability on switch bout durations p = ' num2str(p)])
end
figure;
notBoxPlot([stayBs(:,1); switchBs(:,1)],[ones(size(stayBs,1),1); ones(size(switchBs,1),1)*2])
ylabel('Palatability coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Stay','Switch'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'Stay_vs_switch_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
notBoxPlot([stayBs(:,2); switchBs(:,2)],[ones(size(stayBs,1),1); ones(size(switchBs,1),1)*2])
ylabel('Alternative palatability coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Stay','Switch'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'Stay_vs_switch_alternative_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

[p,h] = signrank(stayBs(:,1),switchBs(:,1));
disp(['Stay vs. Switch palatability coeff comparison: p = ' num2str(p)])

[p,h] = signrank(stayBs(:,2),switchBs(:,2));
disp(['Stay vs. Switch alternative palatability coeff comparison: p = ' num2str(p)])

figure;
notBoxPlot([stayNormalizedBs(:,1); switchNormalizedBs(:,1)],[ones(size(stayNormalizedBs,1),1); ones(size(switchBs,1),1)*2])
ylabel('Normalized palatability coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Stay','Switch'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'Stay_vs_switch_normalized_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
notBoxPlot([stayNormalizedBs(:,2); switchNormalizedBs(:,2)],[ones(size(stayNormalizedBs,1),1); ones(size(switchNormalizedBs,1),1)*2])
ylabel({'Normalized alternative','palatability coefficient'},'fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Stay','Switch'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'Stay_vs_switch_normalized_alternative_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    
    print([figFolder figName],'-dpng','-r600')
end

figure;
notBoxPlot([firstStayBs(:,1); firstStayBs(:,2)],[ones(size(firstStayBs,1),1); ones(size(firstStayBs,1),1)*2])
ylabel('Regression coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Current','Alternative'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
title('First stay bouts only')
if (do_save)
    figName = 'First_stay_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Animal specific linear models of bout duration for early and late bouts
clear earlyX earlyY earlyBs earlyBINTs earlyRs earlyRINTs earlySTATs ...
      lateX lateY lateBs lateBINTs lateRs lateRINTs lateSTATs
earlyInds = find(boutTable.early);
lateInds = find(boutTable.late);
earlyCount = 0;
lateCount = 0;
count = 0;
for i=1:length(animalObjs)
    count = count+1;
    animalInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    earlyX = [];
    lateX = [];
    earlyY = [];
    lateY = [];
    for j=animalObjs(i).difSolutionDays
        dayInds = find(strcmp(boutTable.date,animalObjs(i).dates{j}));
        earlyTableInds = intersect(intersect(animalInds,dayInds),earlyInds);
        lateTableInds = intersect(intersect(animalInds,dayInds),lateInds);
        earlyX = [earlyX; [boutTable.palatability(earlyTableInds) boutTable.alternative_palatability(earlyTableInds)]];
        lateX = [lateX; [boutTable.palatability(lateTableInds) boutTable.alternative_palatability(lateTableInds)]];
        earlyY = [earlyY; boutTable.duration(earlyTableInds)];
        lateY = [lateY; boutTable.duration(lateTableInds)];
        for k=1:length(earlyTableInds)
            earlyCount = earlyCount + 1;
            earlyRat{earlyCount} = animalObjs(i).name;
        end
        for k=1:length(lateTableInds)
            lateCount = lateCount + 1;
            lateRat{lateCount} = animalObjs(i).name;
        end
    end
    earlyX = [earlyX ones(size(earlyX,1),1)];
    [b,bint,r,rint,stats] = regress(earlyY,earlyX);
    earlyBs(count,:) = b;
    earlyNormalizedBs(count,:) = b./mean(earlyY);
    earlyBINTs{count} = bint;
    earlyRs{count} = r;
    earlyRINTs{count} = rint;
    earlySTATs{count} = stats;
   
    lateX = [lateX ones(size(lateX,1),1)];
    [b,bint,r,rint,stats] = regress(lateY,lateX);
    lateBs(count,:) = b;
    lateNormalizedBs(count,:) = b./mean(lateY);
    lateBINTs{count} = bint;
    lateRs{count} = r;
    lateRINTs{count} = rint;
    lateSTATs{count} = stats;
end
[p,h] = signrank(earlyBs(:,1));
if (p < .05)
    disp(['Significant effect of palatability on early bout durations p = ' num2str(p)])
else
    disp(['No significant effect of palatability on early bout durations p = ' num2str(p)])
end
[p,h] = signrank(earlyBs(:,2));
if (p < .05)
    disp(['Significant effect of alternative palatability on early bout duration p = ' num2str(p)])
else
    disp(['No significant effect of alternative palatability on early bout durations p = ' num2str(p)])
end
[p,h] = signrank(lateBs(:,1));
if (p < .05)
    disp(['Significant effect of palatability on late bout durations p = ' num2str(p)])
else
    disp(['No significant effect of palatability on late bout durations p = ' num2str(p)])
end
[p,h] = signrank(lateBs(:,2));
if (p < .05)
    disp(['Significant effect of alternative palatability on late bout duration p = ' num2str(p)])
else
    disp(['No significant effect of alternative palatability on late bout durations p = ' num2str(p)])
end

% Boxplots of early vs. late palatability coefficients
figure;
notBoxPlot([earlyBs(:,1); lateBs(:,1)],[ones(size(earlyBs,1),1); ones(size(lateBs,1),1)*2])
ylabel('Palatability coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'early','late'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])

[p1,h] = signrank(earlyBs(:,1),lateBs(:,1),'tail','right');
[p2,h] = signrank(lateBs(:,1),earlyBs(:,1),'tail','right');
if (p1 < .05)
    disp(['Early palatability coefficients are significantly higher than late coefficients. p = ' num2str(p1)])
elseif (p2 < .05)
    disp(['Late palatability coefficients are significantly higher than early coefficients. p = ' num2str(p2)])
else
    disp(['Neither early or late palatability coefficients are significantly higher than each other'])
end
if (do_save)
    figName = 'Early_vs_late_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

% Boxplot of early vs late alternative palatability coefficients
figure;
notBoxPlot([earlyBs(:,2); lateBs(:,2)],[ones(size(earlyBs,1),1); ones(size(lateBs,1),1)*2])
ylabel('Alternative palatability coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'early','late'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])

[p1,h] = signrank(earlyBs(:,2),lateBs(:,2),'tail','right');
[p2,h] = signrank(lateBs(:,2),earlyBs(:,2),'tail','right');
if (p1 < .05)
    disp(['Early alternative palatability coefficients are significantly higher than late coefficients. p = ' num2str(p1)])
elseif (p2 < .05)
    disp(['Late alternative palatability coefficients are significantly higher than early coefficients. p = ' num2str(p2)])
else
    disp(['Neither early or late alternative palatability coefficients are significantly higher than each other'])
end
if (do_save)
    figName = 'Early_vs_late_alternative_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

% Boxplot of early vs late palatability coefficient magnitudes
figure;
notBoxPlot([abs(earlyBs(:,1)); abs(lateBs(:,1))],[ones(size(earlyBs,1),1); ones(size(lateBs,1),1)*2])
ylabel('Palatablity coefficient magnitude','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Early','Late'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'early_vs_late_palatability_coefficient_magnitudes';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

% Boxplot of early vs late alternative palatability coefficient magnitudes
figure;
notBoxPlot([abs(earlyBs(:,2)); abs(lateBs(:,2))],[ones(size(earlyBs,1),1); ones(size(lateBs,1),1)*2])
ylabel('Alternative palatablity coefficient magnitude','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Early','Late'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'early_vs_late_alternative_palatability_coefficient_magnitudes';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

% Boxplots of normalized early vs. late palatability coefficients
figure;
notBoxPlot([earlyNormalizedBs(:,1); lateNormalizedBs(:,1)],[ones(size(earlyNormalizedBs,1),1); ones(size(lateNormalizedBs,1),1)*2])
ylabel('Normalized palatability coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'early','late'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
[p1,h] = signrank(earlyNormalizedBs(:,1),lateNormalizedBs(:,1));
if (p1 < .05)
    disp(['Early and late palatability coefficients are significantly different. p = ' num2str(p1)])
else
    disp(['Early and late palatability coefficients are not significantly different.'])
end
if (do_save)
    figName = 'Early_vs_late_normalized_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

% Boxplots of normalized early vs. late alternative palatability coefficients
figure;
notBoxPlot([earlyNormalizedBs(:,2); lateNormalizedBs(:,2)],[ones(size(earlyNormalizedBs,1),1); ones(size(lateNormalizedBs,1),1)*2])
ylabel({'Normalized Alternative','palatability coefficient'},'fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'early','late'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
[p1,h] = signrank(earlyNormalizedBs(:,2),lateNormalizedBs(:,2));
if (p1 < .05)
    disp(['Early and late palatability coefficients are signficantly different. p = ' num2str(p1)])
else
    disp(['Early and late palatability coefficients are significantly different.'])
end
if (do_save)
    figName = 'Early_vs_late_normalized_alternative_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% ANOVA of bout duration (with palatabilities as factors) for early bouts
difSolInds = find(boutTable.solutionNum ~= boutTable.alternative_solutionNum);
earlyInds = find(boutTable.early);
inds = intersect(difSolInds,earlyInds);
earlyDurations = boutTable.duration(inds);
earlyPalatabilities = boutTable.palatability(inds);
earlyAlternativePalatabilities = boutTable.alternative_palatability(inds);

factors = {earlyPalatabilities,earlyAlternativePalatabilities};
varnames = {'(early) palatability', '(early) alternative palatability'};
isContinuous = [1 2];
[p,tbl,stats,terms] = anovan(earlyDurations,factors,'model','interaction','continuous',...
                             isContinuous,'varnames',varnames);

[p,tbl,stats,terms] = anovan(earlyDurations,factors,'model','linear','continuous',...
                             isContinuous,'varnames',varnames);
                         
factors = {earlyPalatabilities, earlyAlternativePalatabilities, earlyRat};
varnames = {'(early) palatability', '(early) alternative palatability', 'rat'};
isContinuous = [1 2];
isRandom = [3];
[p,tbl,stats,terms] = anovan(earlyDurations,factors,'model','interaction','continuous',...
                             isContinuous,'random',isRandom,'varnames',varnames);
%% ANOVA of bout duration (with palatabilities as factors) for late bouts
difSolInds = find(boutTable.solutionNum ~= boutTable.alternative_solutionNum);
lateInds = find(boutTable.late);
inds = intersect(difSolInds,lateInds);
lateDurations = boutTable.duration(inds);
latePalatabilities = boutTable.palatability(inds);
lateAlternativePalatabilities = boutTable.alternative_palatability(inds);

factors = {latePalatabilities,lateAlternativePalatabilities};
varnames = {'(late) palatability', '(late) alternative palatability'};
isContinuous = [1 2];
[p,tbl,stats,terms] = anovan(lateDurations,factors,'model','interaction','continuous',...
                             isContinuous,'varnames',varnames);

[p,tbl,stats,terms] = anovan(lateDurations,factors,'model','linear','continuous',...
                             isContinuous,'varnames',varnames);
                         
factors = {latePalatabilities, lateAlternativePalatabilities, lateRat};
varnames = {'(late) palatability', '(late) alternative palatability', 'rat'};
isContinuous = [1 2];
isRandom = [3];
[p,tbl,stats,terms] = anovan(lateDurations,factors,'model','interaction','continuous',...
                             isContinuous,'random',isRandom,'varnames',varnames);

%% Figure 19: Get solution transition matrices using all bouts
boutTransitionProbs = cell(1,length(animalObjs));
allSolnNswitches = cell(1,length(animalObjs));
allNbouts = cell(1,length(animalObjs));
for i=1:length(animalObjs)
    boutTransitionProbs{i} = cell(4,4);
    allSolnNswitches{i} = cell(4,4);
    allNbouts{i} = cell(4,4);
    for j=1:4
        for k=j:4
            if (j == 1 && k > 1)
                continue;
            end
            dayInds = [];
            for day=1:length(animalObjs(i).dates)
                if (animalObjs(i).solutions(day,1) == j && animalObjs(i).solutions(day,2) == k)
                    dayInds = [dayInds day];
                elseif (animalObjs(i).solutions(day,1) == k && animalObjs(i).solutions(day,2) == j)
                    dayInds = [dayInds day];
                end 
            end
            if (length(dayInds) > 2)
                warning([animalObjs(i).name ' solutions ' num2str(j) ' and ' num2str(k) ' too many days identified'])
            end
            tCount = zeros(size(animalObjs(i).solutions,2));
            totalTrans = zeros(1,size(animalObjs(i).solutions,2));
            for d=1:length(dayInds)
                [T,transitionCount,totalTransitions,nswitches,nbouts] = getBoutTransitionMatrix(animalObjs(i).bouts{dayInds(d)},[j k]);
                tCount = tCount + transitionCount;
                totalTrans = totalTrans + totalTransitions;
                allSolnNswitches{i}{j,k} = [allSolnNswitches{i}{j,k} nswitches];
                allNbouts{i}{j,k} = [allNbouts{i}{j,k} nbouts];
            end
            totalT = tCount./totalTrans;
            boutTransitionProbs{i}{j,k} = totalT;
        end
    end
end
solutionTransitionMatrices = cell(1,length(animalObjs));
for i=1:length(animalObjs)
    solutionTransitionMatrices{i} = getSolutionTransitionMatrix(boutTransitionProbs{i});
end
figure;
transitionProbVals = cell(1,10);
for k=setdiff(1:length(animalObjs),5) % exclude bb12
    count = 0;
    for i=1:4
        for j=1:4
            if ((i == 1 || j == 1) && i~=j)
                continue;
            else
                count = count + 1;
                transitionProbVals{count} = [transitionProbVals{count} solutionTransitionMatrices{k}(i,j)];
            end
        end
    end
end
for i=1:length(transitionProbVals)
    scatter(i*ones(1,length(transitionProbVals{i})),transitionProbVals{i},100,'k.')
    hold on;
    plot([i-.2 i+.2],[mean(transitionProbVals{i}) mean(transitionProbVals{i})],'r')
    plot([i-.2 i+.2],[median(transitionProbVals{i}) median(transitionProbVals{i})],'b')
end
xlabels = {'H2O->H2O','A->A','B->A','C->A','A->B','B->B','C->B','A->C','B->C','C->C'};
set(gca,'xtick',1:10,'xticklabels',xlabels,'xticklabelrotation',20,'fontsize',15,'fontweight','bold'); xlim([0 11]); ylim([0 1])
ylabel('Transition probability','fontsize',15,'fontweight','bold')
hold off;
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    saveas(gcf,[figFolder 'transition_probability_by_solution_pair'],'fig')
    saveas(gcf,[figFolder 'transition_probability_by_solution_pair'],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 20: Get solution transition matrices for early bouts only
earlyBoutTransitionProbs = cell(1,length(animalObjs));
earlyAllSolnNswitches = cell(1,length(animalObjs));
earlyAllNbouts = cell(1,length(animalObjs));
for i=1:length(animalObjs)
    earlyBoutTransitionProbs{i} = cell(4,4);
    earlyAllSolnNswitches{i} = cell(4,4);
    earlyAllNbouts{i} = cell(4,4);
    for j=1:4
        for k=j:4
            if (j == 1 && k > 1)
                continue;
            end
            dayInds = [];
            for day=1:length(animalObjs(i).dates)
                if (animalObjs(i).solutions(day,1) == j && animalObjs(i).solutions(day,2) == k)
                    dayInds = [dayInds day];
                elseif (animalObjs(i).solutions(day,1) == k && animalObjs(i).solutions(day,2) == j)
                    dayInds = [dayInds day];
                end 
            end
            if (length(dayInds) > 2)
                warning([animalObjs(i).name ' solutions ' num2str(j) ' and ' num2str(k) ' too many days identified'])
            end
            earlyTCount = zeros(size(animalObjs(i).solutions,2));
            earlyTotalTrans = zeros(1,size(animalObjs(i).solutions,2));
            for d=1:length(dayInds)
                [earlyT,earlyTransitionCount,earlyTotalTransitions,earlyNswitches,earlyNbouts] = getBoutTransitionMatrix(animalObjs(i).earlyBouts{dayInds(d)},[j k]);
                earlyTCount = earlyTCount + earlyTransitionCount;
                earlyTotalTrans = earlyTotalTrans + earlyTotalTransitions;
                earlyAllSolnNswitches{i}{j,k} = [earlyAllSolnNswitches{i}{j,k} earlyNswitches];
                earlyAllNbouts{i}{j,k} = [earlyAllNbouts{i}{j,k} earlyNbouts];
            end
            earlyTotalT = earlyTCount./earlyTotalTrans;
            earlyBoutTransitionProbs{i}{j,k} = earlyTotalT;
        end
    end
end
earlySolutionTransitionMatrices = cell(1,length(animalObjs));
for i=1:length(animalObjs)
    earlySolutionTransitionMatrices{i} = getSolutionTransitionMatrix(earlyBoutTransitionProbs{i});
end
figure;
earlyTransitionProbVals = cell(1,10);
for k=setdiff(1:length(animalObjs),5) % exclude bb12
    count = 0;
    for i=1:4
        for j=1:4
            if ((i == 1 || j == 1) && i~=j)
                continue;
            else
                count = count + 1;
                earlyTransitionProbVals{count} = [earlyTransitionProbVals{count} earlySolutionTransitionMatrices{k}(i,j)];
            end
        end
    end
end
for i=1:length(earlyTransitionProbVals)
    scatter(i*ones(1,length(earlyTransitionProbVals{i})),earlyTransitionProbVals{i},100,'k.')
    hold on;
    plot([i-.2 i+.2],[mean(earlyTransitionProbVals{i},'omitnan') mean(earlyTransitionProbVals{i},'omitnan')],'r')
    plot([i-.2 i+.2],[median(earlyTransitionProbVals{i},'omitnan') median(earlyTransitionProbVals{i},'omitnan')],'b')
end
xlabels = {'H2O->H2O','A->A','B->A','C->A','A->B','B->B','C->B','A->C','B->C','C->C'};
set(gca,'xtick',1:10,'xticklabels',xlabels,'xticklabelrotation',20,'fontsize',15,'fontweight','bold');
xlim([0 11]); ylim([0 1])
ylabel('Transition probability','fontsize',15,'fontweight','bold')
title('Early bouts only','fontsize',15,'fontweight','bold')
hold off;
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'early_transition_probability_by_solution_pair';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 21: Get solution transition matrices for late bouts only
lateBoutTransitionProbs = cell(1,length(animalObjs));
lateAllSolnNswitches = cell(1,length(animalObjs));
lateAllNbouts = cell(1,length(animalObjs));
for i=1:length(animalObjs)
    lateBoutTransitionProbs{i} = cell(4,4);
    lateAllSolnNswitches{i} = cell(4,4);
    lateAllNbouts{i} = cell(4,4);
    for j=1:4
        for k=j:4
            if (j == 1 && k > 1)
                continue;
            end
            dayInds = [];
            for day=1:length(animalObjs(i).dates)
                if (animalObjs(i).solutions(day,1) == j && animalObjs(i).solutions(day,2) == k)
                    dayInds = [dayInds day];
                elseif (animalObjs(i).solutions(day,1) == k && animalObjs(i).solutions(day,2) == j)
                    dayInds = [dayInds day];
                end 
            end
            if (length(dayInds) > 2)
                warning([animalObjs(i).name ' solutions ' num2str(j) ' and ' num2str(k) ' too many days identified'])
            end
            lateTCount = zeros(size(animalObjs(i).solutions,2));
            lateTotalTrans = zeros(1,size(animalObjs(i).solutions,2));
            for d=1:length(dayInds)
                [lateT,lateTransitionCount,lateTotalTransitions,lateNswitches,lateNbouts] = getBoutTransitionMatrix(animalObjs(i).lateBouts{dayInds(d)},[j k]);
                lateTCount = lateTCount + lateTransitionCount;
                lateTotalTrans = lateTotalTrans + lateTotalTransitions;
                lateAllSolnNswitches{i}{j,k} = [lateAllSolnNswitches{i}{j,k} lateNswitches];
                lateAllNbouts{i}{j,k} = [lateAllNbouts{i}{j,k} lateNbouts];
            end
            lateTotalT = lateTCount./lateTotalTrans;
            lateBoutTransitionProbs{i}{j,k} = lateTotalT;
        end
    end
end
lateSolutionTransitionMatrices = cell(1,length(animalObjs));
for i=1:length(animalObjs)
    lateSolutionTransitionMatrices{i} = getSolutionTransitionMatrix(lateBoutTransitionProbs{i});
end
figure;
lateTransitionProbVals = cell(1,10);
for k=setdiff(1:length(animalObjs),5) % exclude bb12
    count = 0;
    for i=1:4
        for j=1:4
            if ((i == 1 || j == 1) && i~=j)
                continue;
            else
                count = count + 1;
                lateTransitionProbVals{count} = [lateTransitionProbVals{count} lateSolutionTransitionMatrices{k}(i,j)];
            end
        end
    end
end
for i=1:length(lateTransitionProbVals)
    scatter(i*ones(1,length(lateTransitionProbVals{i})),lateTransitionProbVals{i},100,'k.')
    hold on;
    plot([i-.2 i+.2],[mean(lateTransitionProbVals{i},'omitnan') mean(lateTransitionProbVals{i},'omitnan')],'r')
    plot([i-.2 i+.2],[median(lateTransitionProbVals{i},'omitnan') median(lateTransitionProbVals{i},'omitnan')],'b')
end
xlabels = {'H2O->H2O','A->A','B->A','C->A','A->B','B->B','C->B','A->C','B->C','C->C'};
set(gca,'xtick',1:10,'xticklabels',xlabels,'xticklabelrotation',20,'fontsize',15,'fontweight','bold');
xlim([0 11]); ylim([0 1])
ylabel('Transition probability','fontsize',15,'fontweight','bold')
title('Late bouts only','fontsize',15,'fontweight','bold')
hold off;
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'late_transition_probability_by_solution_pair';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 22: Plot transition probability as function of current - alternative palatability
slns = [1 1; 2 2; 3 2; 4 2; 2 3; 3 3; 4 3; 2 4; 3 4; 4 4];
relPalDifs = cell(1,length(transitionProbVals));
altPals = cell(1,size(slns,1));
curPals = cell(1,size(slns,1));
for i=1:size(slns,1)
    count = 0;
    for k=setdiff(1:length(animalObjs),5) % exclude bb12
        count=count+1;
        relPalDifs{i}(count) = animalObjs(k).relativePalatabilitiesLicks(slns(i,1)) - animalObjs(k).relativePalatabilitiesLicks(slns(i,2));
        curPals{i}(count) = animalObjs(k).relativePalatabilitiesLicks(slns(i,1));
        altPals{i}(count) = animalObjs(k).relativePalatabilitiesLicks(slns(i,2));
    end
end
allTransitionProbVals = [];
allRelPalDifs = [];
allCurPals = [];
allAltPals = [];
for i=1:length(transitionProbVals)
    allTransitionProbVals = [allTransitionProbVals transitionProbVals{i}];
    allRelPalDifs = [allRelPalDifs relPalDifs{i}];
    allCurPals = [allCurPals curPals{i}];
    allAltPals = [allAltPals altPals{i}];
end
lf1 = polyfit(allRelPalDifs,allTransitionProbVals,1);
lf2 = polyfit(allCurPals,allTransitionProbVals,1);
lf3 = polyfit(allAltPals,allTransitionProbVals,1);

[rho1,pval1] = corr(allRelPalDifs',allTransitionProbVals');
[rho2,pval2] = corr(allCurPals',allTransitionProbVals');
[rho3,pval3] = corr(allAltPals',allTransitionProbVals');
figure;
scatter(allRelPalDifs,allTransitionProbVals,100,'k.')
x=min(allRelPalDifs):.01:max(allRelPalDifs);
hold on;
plot(x,lf1(1)*x + lf1(2))
xlabel('Current - alternative palatability','fontsize',15,'fontweight','bold');
ylabel('Transition Probability','fontsize',15,'fontweight','bold')
text(0,.95,['\rho = ' num2str(rho1)],'fontweight','bold','fontsize',20)
text(0,.88,['p = ' num2str(pval1)],'fontweight','bold','fontsize',20)
text(0,.81,['y = ' num2str(lf1(1)) 'x + ' num2str(lf1(2))],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'transition_probability_vs_palatability_difference_linfit';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 23: Plot transition probability as a function of palatability and alternative palatability
figure;
scatter(allCurPals,allTransitionProbVals,100,'k.')
x=min(allCurPals):.01:max(allCurPals);
hold on;
plot(x,lf2(1)*x + lf2(2))
xlabel('Current palatability'); ylabel('Transition Probability')
text(1.7,.95,['\rho = ' num2str(rho2)],'fontweight','bold','fontsize',20)
text(1.7,.88,['p = ' num2str(pval2)],'fontweight','bold','fontsize',20)
%text(1.7,.81,['y = ' num2str(lf2(1)) 'x + ' num2str(lf2(2))],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    saveas(gcf,[figFolder 'transition_probability_vs_current_palatability'],'fig')
    saveas(gcf,[figFolder 'transition_probability_vs_current_palatability'],'eps')
    saveas(gcf,[figFolder 'transition_probability_vs_current_palatability'],'png')
end

figure;
scatter(allAltPals,allTransitionProbVals,100,'k.')
x=min(allAltPals):.01:max(allAltPals);
hold on;
plot(x,lf3(1)*x + lf3(2))
xlabel('Alternative palatability','fontsize',15,'fontweight','bold');
ylabel('Transition Probability','fontsize',15,'fontweight','bold')
text(0.3,.95,['\rho = ' num2str(rho3)],'fontweight','bold','fontsize',20)
text(0.3,.88,['p = ' num2str(pval3)],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'transition_probability_vs_alternative_palatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Linear model for transition probabilities
clear YY XX;
solnPairs = [3 2;
             2 3;
             2 4;
             4 2;
             3 4;
             4 3];
count = 0;
for i=1:length(animalObjs)
    for j=1:size(solnPairs,1)
        count = count + 1;
        YY(count) = solutionTransitionMatrices{i}(solnPairs(j,1),solnPairs(j,2));
        XX(count,1) = animalObjs(i).relativePalatabilitiesLicks(solnPairs(j,2));
        XX(count,2) = animalObjs(i).relativePalatabilitiesLicks(solnPairs(j,1));
        rat{count} = animalObjs(i).name;
    end
end
XX = [XX ones(size(XX,1),1)];
[b,bint,r,rint,stats] = regress(YY',XX);

%% Transition probability ANOVA
factors = {XX(:,1), XX(:,2)};
isContinuous = [1 2];
varnames = {'palatability', 'alternative_palatability'};
[p,tbl,stats,terms] = anovan(YY,factors,'model','interaction','continuous',...
                             isContinuous,'varnames',varnames);
                         
%% Transition probability ANOVA with rat as random factor
factors = {XX(:,1), XX(:,2),rat};
isContinuous = [1 2];
isRandom = 3;
varnames = {'palatability', 'alternative_palatability','rat'};
[p,tbl,stats,terms] = anovan(YY,factors,'model','interaction','continuous',...
                             isContinuous,'varnames',varnames);

%% Logistic regression on switch probability
clear wasSwitch pastPalatability alternativePalatability pastDuration numConsecutive totalDuration
wasSwitch = [];
pastPalatability = [];
alternativePalatability = [];
pastDuration = [];
numConsecutive = [];
totalDuration = [];
for i=1:length(animalObjs)
    animalInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    for j=animalObjs(i).difSolutionDays
        dayInds = find(strcmp(boutTable.date,animalObjs(i).dates{j}));
        tableInds = intersect(animalInds,dayInds);
        [~,sortedOrder] = sort(boutTable.boutNumber(tableInds));
        sortedInds = tableInds(sortedOrder);
        nconsecutive = 0;
        curTotalDuration = 0;
        for k=1:(length(sortedInds) - 1)
            curTotalDuration = curTotalDuration + boutTable.duration(sortedInds(k));
            pastPalatability = [pastPalatability boutTable.palatability(sortedInds(k))];
            alternativePalatability = [alternativePalatability boutTable.alternative_palatability(sortedInds(k))];
            pastDuration = [pastDuration boutTable.duration(sortedInds(k))];
            numConsecutive = [numConsecutive nconsecutive];
            totalDuration = [totalDuration curTotalDuration];
            if (boutTable.channel(sortedInds(k)) ~= boutTable.channel(sortedInds(k+1))) % If the animal switches
                wasSwitch = [wasSwitch 1];
                nconsecutive = 0;
            else
                wasSwitch = [wasSwitch 0];
                nconsecutive = nconsecutive + 1;
            end
        end
    end
end

glm = fitglm([pastPalatability' alternativePalatability'],wasSwitch','Distribution','binomial','Link','logit');
%glm = fitglm([pastPalatability' alternativePalatability'],wasSwitch','interactions','Distribution','binomial','Link','logit');
%glm = fitglm([pastPalatability' alternativePalatability' pastDuration'],wasSwitch','interactions','Distribution','binomial','Link','logit');
%glm = fitglm([pastPalatability' alternativePalatability' pastDuration' numConsecutive'],wasSwitch','interactions','Distribution','binomial','Link','logit');
%glm = fitglm([pastPalatability' alternativePalatability' pastDuration' numConsecutive' totalDuration'],wasSwitch','interactions','Distribution','binomial','Link','logit');
glm.plotSlice();
glm.plotResiduals('probability')
glm
glm.Rsquared


palatabilityVals = 0:.01:3;
alternative_palatabilityVals = 0:.01:3;
predictedVals = zeros(length(palatabilityVals),length(alternative_palatabilityVals));
upperCIVals = zeros(length(palatabilityVals),length(alternative_palatabilityVals));
lowerCIVals = zeros(length(palatabilityVals),length(alternative_palatabilityVals));
totalEvals = length(palatabilityVals)*length(alternative_palatabilityVals);
evalCount = 0;
for i=1:length(palatabilityVals)
    for j=1:length(palatabilityVals)
        [y,yci] = glm.predict([palatabilityVals(i),alternative_palatabilityVals(j)]);
        predictedVals(i,j) = y;
        lowerCIVals(i,j) = yci(1);
        upperCIVals(i,j) = yci(2);
        evalCount = evalCount + 1;
        if (mod(evalCount,round(totalEvals/10)) == 0)
            disp([num2str((evalCount/totalEvals)*100) '% done'])
        end
    end
end

%% Surf plot of switch probability as function of palatability and alternative palatability with CIs
figure;
surf(palatabilityVals,alternative_palatabilityVals,predictedVals','FaceColor','g','EdgeColor','b','FaceAlpha',.3)
hold on;
surf(palatabilityVals,alternative_palatabilityVals,lowerCIVals','LineStyle',':','FaceColor','r','EdgeColor','r','FaceAlpha',.2)
surf(palatabilityVals,alternative_palatabilityVals,upperCIVals','LineStyle',':','FaceColor','r','EdgeColor','r','FaceAlpha',.2)
xlabel('Current palatability','fontsize',15,'fontweight','bold'); 
ylabel('Alternative palatability','fontsize',15,'fontweight','bold');
zlabel('Switch probability','fontsize',15,'fontweight','bold')

%% Linear model for early transition probabilities
clear earlyYY earlyXX earlyRat;
count = 0;
for i=1:length(animalObjs)
    for j=1:size(solnPairs,1)
        count = count + 1;
        earlyYY(count) = earlySolutionTransitionMatrices{i}(solnPairs(j,1),solnPairs(j,2));
        earlyXX(count,1) = animalObjs(i).relativePalatabilitiesLicks(solnPairs(j,2));
        earlyXX(count,2) = animalObjs(i).relativePalatabilitiesLicks(solnPairs(j,1));
        earlyRat{count} = animalObjs(i).name;
    end
end
earlyXX = [earlyXX ones(size(earlyXX,1),1)];
[b,bint,r,rint,stats] = regress(earlyYY',earlyXX);

%% Early Transition probability ANOVA
factors = {earlyXX(:,1), earlyXX(:,2)};
isContinuous = [1 2];
varnames = {'(early) palatability', '(early) alternative_palatability'};
[p,tbl,stats,terms] = anovan(earlyYY,factors,'model','interaction','continuous',...
                             isContinuous,'varnames',varnames);
                         
%% Early Transition probability ANOVA with rat as random factor
factors = {earlyXX(:,1), earlyXX(:,2),earlyRat};
isContinuous = [1 2];
isRandom = 3;
varnames = {'(early) palatability', '(early) alternative_palatability','rat'};
[p,tbl,stats,terms] = anovan(earlyYY,factors,'model','interaction','continuous',...
                             isContinuous,'varnames',varnames);

%% Logistic regression on switch probability for early bouts
clear earlyWasSwitch earlyPastPalatability earlyAlternativePalatability earlyPastDuration ...
      earlyNumConsecutive earlyTotalDuration
earlyWasSwitch = [];
earlyPastPalatability = [];
earlyAlternativePalatability = [];
earlyPastDuration = [];
earlyNumConsecutive = [];
earlyTotalDuration = [];
for i=1:length(animalObjs)
    animalInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    for j=animalObjs(i).difSolutionDays
        dayInds = find(strcmp(boutTable.date,animalObjs(i).dates{j}));
        earlyInds = find(boutTable.early);
        tableInds = intersect(intersect(animalInds,dayInds),earlyInds);
        [~,sortedOrder] = sort(boutTable.boutNumber(tableInds));
        sortedInds = tableInds(sortedOrder);
        nconsecutive = 0;
        curTotalDuration = 0;
        for k=1:(length(sortedInds) - 1)
            curTotalDuration = curTotalDuration + boutTable.duration(sortedInds(k));
            earlyPastPalatability = [earlyPastPalatability boutTable.palatability(sortedInds(k))];
            earlyAlternativePalatability = [earlyAlternativePalatability boutTable.alternative_palatability(sortedInds(k))];
            earlyPastDuration = [earlyPastDuration boutTable.duration(sortedInds(k))];
            earlyNumConsecutive = [earlyNumConsecutive nconsecutive];
            earlyTotalDuration = [earlyTotalDuration curTotalDuration];
            if (boutTable.channel(sortedInds(k)) ~= boutTable.channel(sortedInds(k+1))) % If the animal switches
                earlyWasSwitch = [earlyWasSwitch 1];
                nconsecutive = 0;
            else
                earlyWasSwitch = [earlyWasSwitch 0];
                nconsecutive = nconsecutive + 1;
            end
        end
    end
end

glm = fitglm([earlyPastPalatability' earlyAlternativePalatability'],earlyWasSwitch','interactions','Distribution','binomial','Link','logit');
%glm = fitglm([earlyPastPalatability' earlyAlternativePalatability' earlyPastDuration'],earlyWasSwitch','interactions','Distribution','binomial','Link','logit');
%glm = fitglm([earlyPastPalatability' earlyAlternativePalatability' earlyPastDuration' earlyNumConsecutive'],earlyWasSwitch','interactions','Distribution','binomial','Link','logit');
%glm = fitglm([earlyPastPalatability' earlyAlternativePalatability' earlyPastDuration' earlyNumConsecutive' earlyTotalDuration'],earlyWasSwitch','interactions','Distribution','binomial','Link','logit');
glm.plotSlice();
glm.plotResiduals('probability')
glm
glm.Rsquared

palatabilityVals = 0:.01:3;
alternative_palatabilityVals = 0:.01:3;
earlyPredictedVals = zeros(length(palatabilityVals),length(alternative_palatabilityVals));
earlyUpperCIVals = zeros(length(palatabilityVals),length(alternative_palatabilityVals));
earlyLowerCIVals = zeros(length(palatabilityVals),length(alternative_palatabilityVals));
totalEvals = length(palatabilityVals)*length(alternative_palatabilityVals);
evalCount = 0;
for i=1:length(palatabilityVals)
    for j=1:length(palatabilityVals)
        [y,yci] = glm.predict([palatabilityVals(i),alternative_palatabilityVals(j)]);
        earlyPredictedVals(i,j) = y;
        earlyLowerCIVals(i,j) = yci(1);
        earlyUpperCIVals(i,j) = yci(2);
        evalCount = evalCount + 1;
        if (mod(evalCount,round(totalEvals/10)) == 0)
            disp([num2str((evalCount/totalEvals)*100) '% done'])
        end
    end
end
%% Linear model for late transition probabilities
clear lateYY lateXX lateRat;
count = 0;
for i=1:length(animalObjs)
    for j=1:size(solnPairs,1)
        count = count + 1;
        lateYY(count) = lateSolutionTransitionMatrices{i}(solnPairs(j,1),solnPairs(j,2));
        lateXX(count,1) = animalObjs(i).relativePalatabilitiesLicks(solnPairs(j,2));
        lateXX(count,2) = animalObjs(i).relativePalatabilitiesLicks(solnPairs(j,1));
        lateRat{count} = animalObjs(i).name;
    end
end
lateXX = [lateXX ones(size(lateXX,1),1)];
[b,bint,r,rint,stats] = regress(lateYY',lateXX);

%% Late Transition probability ANOVA
factors = {lateXX(:,1), lateXX(:,2)};
isContinuous = [1 2];
varnames = {'(late) palatability', '(late) alternative_palatability'};
[p,tbl,stats,terms] = anovan(lateYY,factors,'model','interaction','continuous',...
                             isContinuous,'varnames',varnames);
                         
%% Late Transition probability ANOVA with rat as random factor
factors = {lateXX(:,1), lateXX(:,2),lateRat};
isContinuous = [1 2];
isRandom = 3;
varnames = {'(late) palatability', '(late) alternative_palatability','rat'};
[p,tbl,stats,terms] = anovan(lateYY,factors,'model','interaction','continuous',...
                             isContinuous,'varnames',varnames);

%% Logistic regression on switch probability for late bouts
clear lateWasSwitch latePastPalatability lateAlternativePalatability latePastDuration ...
      lateNumConsecutive lateTotalDuration
lateWasSwitch = [];
latePastPalatability = [];
lateAlternativePalatability = [];
latePastDuration = [];
lateNumConsecutive = [];
lateTotalDuration = [];
for i=1:length(animalObjs)
    animalInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    for j=animalObjs(i).difSolutionDays
        dayInds = find(strcmp(boutTable.date,animalObjs(i).dates{j}));
        lateInds = find(boutTable.late);
        tableInds = intersect(intersect(animalInds,dayInds),lateInds);
        [~,sortedOrder] = sort(boutTable.boutNumber(tableInds));
        sortedInds = tableInds(sortedOrder);
        nconsecutive = 0;
        curTotalDuration = 0;
        for k=1:(length(sortedInds) - 1)
            curTotalDuration = curTotalDuration + boutTable.duration(sortedInds(k));
            latePastPalatability = [latePastPalatability boutTable.palatability(sortedInds(k))];
            lateAlternativePalatability = [lateAlternativePalatability boutTable.alternative_palatability(sortedInds(k))];
            latePastDuration = [latePastDuration boutTable.duration(sortedInds(k))];
            lateNumConsecutive = [lateNumConsecutive nconsecutive];
            lateTotalDuration = [lateTotalDuration curTotalDuration];
            if (boutTable.channel(sortedInds(k)) ~= boutTable.channel(sortedInds(k+1))) % If the animal switches
                lateWasSwitch = [lateWasSwitch 1];
                nconsecutive = 0;
            else
                lateWasSwitch = [lateWasSwitch 0];
                nconsecutive = nconsecutive + 1;
            end
        end
    end
end

glm = fitglm([latePastPalatability' lateAlternativePalatability'],lateWasSwitch','interactions','Distribution','binomial','Link','logit');
%glm = fitglm([latePastPalatability' lateAlternativePalatability' latePastDuration'],lateWasSwitch','interactions','Distribution','binomial','Link','logit');
%glm = fitglm([latePastPalatability' lateAlternativePalatability' latePastDuration' lateNumConsecutive'],lateWasSwitch','interactions','Distribution','binomial','Link','logit');
%glm = fitglm([latePastPalatability' lateAlternativePalatability' latePastDuration' lateNumConsecutive' lateTotalDuration'],lateWasSwitch','interactions','Distribution','binomial','Link','logit');
glm.plotSlice();
glm.plotResiduals('probability')
glm
glm.Rsquared

palatabilityVals = 0:.01:3;
alternative_palatabilityVals = 0:.01:3;
latePredictedVals = zeros(length(palatabilityVals),length(alternative_palatabilityVals));
lateUpperCIVals = zeros(length(palatabilityVals),length(alternative_palatabilityVals));
lateLowerCIVals = zeros(length(palatabilityVals),length(alternative_palatabilityVals));
totalEvals = length(palatabilityVals)*length(alternative_palatabilityVals);
evalCount = 0;
for i=1:length(palatabilityVals)
    for j=1:length(palatabilityVals)
        [y,yci] = glm.predict([palatabilityVals(i),alternative_palatabilityVals(j)]);
        latePredictedVals(i,j) = y;
        lateLowerCIVals(i,j) = yci(1);
        lateUpperCIVals(i,j) = yci(2);
        evalCount = evalCount + 1;
        if (mod(evalCount,round(totalEvals/10)) == 0)
            disp([num2str((evalCount/totalEvals)*100) '% done'])
        end
    end
end

%% Transition probabilities following stay or switch
allPostStayTransitions = [];
allPostSwitchTransitions = [];
allPostStayCurPals = [];
allPostSwitchCurPals = [];
allPostStayAltPals = [];
allPostSwitchAltPals = [];
stayGLMCurPalCoeffs = [];
stayGLMAltPalCoeffs = [];
switchGLMCurPalCoeffs = [];
switchGLMAltPalCoeffs = [];
%for i=1:length(animalObjs)
for i=setdiff(1:length(animalObjs),5)
    postStayTransition = [];
    postSwitchTransition = [];
    postStayCurPal = [];
    postSwitchCurPal = [];
    postStayAltPal = [];
    postSwitchAltPal = [];
    ratInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    for j=animalObjs(i).difSolutionDays
        dayInds = find(strcmp(boutTable.date,animalObjs(i).dates{j}));
        curTableInds = intersect(ratInds,dayInds);
        [~,sortedOrder] = sort(boutTable.onset(curTableInds));
        sortedInds = curTableInds(sortedOrder);
        for k=2:(length(sortedInds)-1)
            if (boutTable.solutionNum(sortedInds(k)) == boutTable.solutionNum(sortedInds(k-1)))
                % following a stay decision
                if (boutTable.solutionNum(sortedInds(k)) == boutTable.solutionNum(sortedInds(k+1)))
                    postStayTransition = [postStayTransition 0];
                else
                    postStayTransition = [postStayTransition 1];
                end
                postStayCurPal = [postStayCurPal boutTable.palatability(sortedInds(k))];
                postStayAltPal = [postStayAltPal boutTable.alternative_palatability(sortedInds(k))];
            else
                % following a switch decision
                if (boutTable.solutionNum(sortedInds(k)) == boutTable.solutionNum(sortedInds(k+1)))
                    postSwitchTransition = [postSwitchTransition 0];
                else
                    postSwitchTransition = [postSwitchTransition 1];
                end
                postSwitchCurPal = [postSwitchCurPal boutTable.palatability(sortedInds(k))];
                postSwitchAltPal = [postSwitchAltPal boutTable.alternative_palatability(sortedInds(k))];
            end
        end
    end
    stayGLMs{i} = fitglm([postStayCurPal' postStayAltPal'],postStayTransition','interactions','Distribution','binomial','Link','logit');
    switchGLMs{i} = fitglm([postSwitchCurPal' postSwitchAltPal'],postSwitchTransition','interactions','Distribution','binomial','Link','logit');
    allPostStayTransitions = [allPostStayTransitions postStayTransition];
    allPostSwitchTransitions = [allPostSwitchTransitions postSwitchTransition];
    allPostStayCurPals = [allPostStayCurPals postStayCurPal];
    allPostSwitchCurPals = [allPostSwitchCurPals postSwitchCurPal];
    allPostStayAltPals = [allPostStayAltPals postStayAltPal];
    allPostSwitchAltPals = [allPostSwitchAltPals postSwitchAltPal];
    stayGLMCurPalCoeffs = [stayGLMCurPalCoeffs stayGLMs{i}.Coefficients(2,1).Variables];
    stayGLMAltPalCoeffs = [stayGLMAltPalCoeffs stayGLMs{i}.Coefficients(3,1).Variables];
    switchGLMCurPalCoeffs = [switchGLMCurPalCoeffs switchGLMs{i}.Coefficients(2,1).Variables];
    switchGLMAltPalCoeffs = [switchGLMAltPalCoeffs switchGLMs{i}.Coefficients(3,1).Variables];
end
allStayGLM = fitglm([allPostStayCurPals' allPostStayAltPals'],allPostStayTransitions,'interactions','Distribution','binomial','Link','logit');
allSwitchGLM = fitglm([allPostSwitchCurPals' allPostSwitchAltPals'],allPostSwitchTransitions,'interactions','Distribution','binomial','Link','logit');

figure;
subplot(1,2,1)
notBoxPlot([stayGLMCurPalCoeffs switchGLMCurPalCoeffs],[ones(1,length(stayGLMCurPalCoeffs)) ones(1,length(switchGLMCurPalCoeffs))*2])
set(gca,'xtick',[1 2],'xticklabels',{'Stay','Switch'},'fontsize',15,'fontweight','bold')
ylabel('Palatability coefficient','fontsize',15,'fontweight','bold')
subplot(1,2,2)
notBoxPlot([stayGLMAltPalCoeffs switchGLMAltPalCoeffs],[ones(1,length(stayGLMAltPalCoeffs)) ones(1,length(switchGLMAltPalCoeffs))*2])
set(gca,'xtick',[1 2],'xticklabels',{'Stay','Switch'},'fontsize',15,'fontweight','bold')
ylabel('Alternative palatability coefficient','fontsize',15,'fontweight','bold')
if (do_save)
    figName = 'logistic_regression_coefficients_palatability_alternative_palatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end
%% Plot early and late switch probability surfaces
surf(palatabilityVals,alternative_palatabilityVals,earlyPredictedVals','FaceColor','b','EdgeColor','b','FaceAlpha',.3)
hold on;
surf(palatabilityVals,alternative_palatabilityVals,latePredictedVals','FaceColor','r','EdgeColor','r','FaceAlpha',.3)
xlabel('Current palatability','fontsize',15,'fontweight','bold'); 
ylabel('Alternative palatability','fontsize',15,'fontweight','bold'); 
zlabel('Switch probability','fontsize',15,'fontweight','bold')
legend({'Early','Late'},'fontsize',15,'fontweight','bold')

%% Compare common source/target transition probabilities
%
% 3 transition matrices per animal
%   (1 - Pab)     Pba           (1 - Pac)      Pca          (1 - Pbc)     Pcb
%      Pab      (1 - Pba)          Pac      (1 - Pca)          Pbc      (1 - Pcb)
%
% Common source pairs:
%   Pab,Pac
%   Pba,Pbc
%   Pca,Pcb
% Common target pairs:
%   Pba,Pca
%   Pab,Pcb
%   Pac,Pbc
clear commonSourceTransVals
solnPairInds = [2 3; 2 4; 3 4];
commonSourcePairInds = {[2 3 2 1; 2 4 2 1],[2 3 1 2; 3 4 2 1],[2 4 1 2; 3 4 1 2]};
commonTargetPairInds = {[2 3 1 2; 2 4 1 2],[2 3 2 1; 3 4 1 2],[2 4 2 1; 3 4 2 1]};
sourcePairTitles = {'P_{a -> b} vs. P_{a -> c}','P_{b -> a} vs. P_{b -> c}','P_{c -> a} vs. P_{c -> b}'};
targetPairTitles = {'P_{b -> a} vs. P_{c -> a}','P_{a -> b} vs. P_{c -> b}','P_{a -> c} vs. P_{b -> c}'};
sourcePairStrings = {{'P_{a -> b}','P_{a -> c}'},{'P_{b -> a}','P_{b -> c}'},{'P_{c -> a}','P_{c -> b}'}};
targetPairStrings = {{'P_{b -> a}','P_{c -> a}'},{'P_{a -> b}','P_{c -> b}'},{'P_{a -> c}','P_{b -> c}'}};
commonSourceTransVals = cell(length(commonSourcePairInds),2);
earlyCommonSourceTransVals = cell(length(commonSourcePairInds),2);
lateCommonSourceTransVals = cell(length(commonSourcePairInds),2);
for i=1:length(commonSourcePairInds)
    for j=1:length(animalObjs)
        % All bouts
        curTransMat1 = boutTransitionProbs{j}{commonSourcePairInds{i}(1,1),commonSourcePairInds{i}(1,2)};
        curTransMat2 = boutTransitionProbs{j}{commonSourcePairInds{i}(2,1),commonSourcePairInds{i}(2,2)};
        commonSourceTransVals{i,1} = [commonSourceTransVals{i,1} curTransMat1(commonSourcePairInds{i}(1,3),commonSourcePairInds{i}(1,4))];
        commonSourceTransVals{i,2} = [commonSourceTransVals{i,2} curTransMat2(commonSourcePairInds{i}(2,3),commonSourcePairInds{i}(2,4))];
        % Early bouts
        curTransMat1 = earlyBoutTransitionProbs{j}{commonSourcePairInds{i}(1,1),commonSourcePairInds{i}(1,2)};
        curTransMat2 = earlyBoutTransitionProbs{j}{commonSourcePairInds{i}(2,1),commonSourcePairInds{i}(2,2)};
        earlyCommonSourceTransVals{i,1} = [earlyCommonSourceTransVals{i,1} curTransMat1(commonSourcePairInds{i}(1,3),commonSourcePairInds{i}(1,4))];
        earlyCommonSourceTransVals{i,2} = [earlyCommonSourceTransVals{i,2} curTransMat2(commonSourcePairInds{i}(2,3),commonSourcePairInds{i}(2,4))];
        % Late bouts
        curTransMat1 = lateBoutTransitionProbs{j}{commonSourcePairInds{i}(1,1),commonSourcePairInds{i}(1,2)};
        curTransMat2 = lateBoutTransitionProbs{j}{commonSourcePairInds{i}(2,1),commonSourcePairInds{i}(2,2)};
        lateCommonSourceTransVals{i,1} = [lateCommonSourceTransVals{i,1} curTransMat1(commonSourcePairInds{i}(1,3),commonSourcePairInds{i}(1,4))];
        lateCommonSourceTransVals{i,2} = [lateCommonSourceTransVals{i,2} curTransMat2(commonSourcePairInds{i}(2,3),commonSourcePairInds{i}(2,4))];
    end
end
commonTargetTransVals = cell(length(commonTargetPairInds),2);
earlyCommonTargetTransVals = cell(length(commonTargetPairInds),2);
lateCommonTargetTransVals = cell(length(commonTargetPairInds),2);
for i=1:length(commonTargetPairInds)
    for j=1:length(animalObjs)
        % All bouts
        curTransMat1 = boutTransitionProbs{j}{commonTargetPairInds{i}(1,1),commonTargetPairInds{i}(1,2)};
        curTransMat2 = boutTransitionProbs{j}{commonTargetPairInds{i}(2,1),commonTargetPairInds{i}(2,2)}; 
        commonTargetTransVals{i,1} = [commonTargetTransVals{i,1} curTransMat1(commonTargetPairInds{i}(1,3),commonTargetPairInds{i}(1,4))];
        commonTargetTransVals{i,2} = [commonTargetTransVals{i,2} curTransMat2(commonTargetPairInds{i}(2,3),commonTargetPairInds{i}(2,4))];
        % Early bouts
        curTransMat1 = earlyBoutTransitionProbs{j}{commonTargetPairInds{i}(1,1),commonTargetPairInds{i}(1,2)};
        curTransMat2 = earlyBoutTransitionProbs{j}{commonTargetPairInds{i}(2,1),commonTargetPairInds{i}(2,2)}; 
        earlyCommonTargetTransVals{i,1} = [earlyCommonTargetTransVals{i,1} curTransMat1(commonTargetPairInds{i}(1,3),commonTargetPairInds{i}(1,4))];
        earlyCommonTargetTransVals{i,2} = [earlyCommonTargetTransVals{i,2} curTransMat2(commonTargetPairInds{i}(2,3),commonTargetPairInds{i}(2,4))];
        % Late bouts
        curTransMat1 = lateBoutTransitionProbs{j}{commonTargetPairInds{i}(1,1),commonTargetPairInds{i}(1,2)};
        curTransMat2 = lateBoutTransitionProbs{j}{commonTargetPairInds{i}(2,1),commonTargetPairInds{i}(2,2)}; 
        lateCommonTargetTransVals{i,1} = [lateCommonTargetTransVals{i,1} curTransMat1(commonTargetPairInds{i}(1,3),commonTargetPairInds{i}(1,4))];
        lateCommonTargetTransVals{i,2} = [lateCommonTargetTransVals{i,2} curTransMat2(commonTargetPairInds{i}(2,3),commonTargetPairInds{i}(2,4))];
    end
end
figure;
for i=1:length(commonSourcePairInds)
    subplot(1,3,i)
    [p1,h1] = signrank(commonSourceTransVals{i,1},commonSourceTransVals{i,2},'tail','right');
    [p2,h2] = signrank(commonSourceTransVals{i,1},commonSourceTransVals{i,2},'tail','left');
    [p3,h3] = signrank(commonSourceTransVals{i,1},commonSourceTransVals{i,2});
    p1 = p1*6; % Bonferroni correction
    p2 = p2*6; % Bonferroni correction
    if (p1 < .05)
        curText = [sourcePairStrings{i}{1} ' > ' sourcePairStrings{i}{2} ' p = ' num2str(p1)];
    elseif (p2 < .05)
        curText = [sourcePairStrings{i}{1} ' < ' sourcePairStrings{i}{2} ' p = ' num2str(p2)];
    else
        curText = [sourcePairStrings{i}{1} ' = ' sourcePairStrings{i}{2} ' p = ' num2str(p3)];
    end
    notBoxPlot([commonSourceTransVals{i,1} commonSourceTransVals{i,2}],[ones(1,length(commonSourceTransVals{i,1})) ones(1,length(commonSourceTransVals{i,2}))*2]);
    ylim([0 1])
    ylabel('Transition probability','fontsize',15,'fontweight','bold')
    set(gca,'xtick',[1 2],'xticklabels',{sourcePairStrings{i}{1},sourcePairStrings{i}{2}},'fontsize',15,'fontweight','bold')
    title(curText,'fontsize',15,'fontweight','bold')
end
suptitle('Common source')
set(gcf,'Position',[10 10 1600 1000])
if (do_save)
    figName = 'common_source';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
for i=1:length(commonTargetPairInds)
    subplot(1,3,i)
    [p1,h1] = signrank(commonTargetTransVals{i,1},commonTargetTransVals{i,2},'tail','right');
    [p2,h2] = signrank(commonTargetTransVals{i,1},commonTargetTransVals{i,2},'tail','left');
    [p3,h3] = signrank(commonTargetTransVals{i,1},commonTargetTransVals{i,2});
    p1 = p1*6; % Bonferroni correction
    p2 = p2*6; % Bonferroni correction
    if (p1 < .05)
        curText = [targetPairStrings{i}{1} ' > ' targetPairStrings{i}{2} ' p = ' num2str(p1)];
    elseif (p2 < .05)
        curText = [targetPairStrings{i}{1} ' < ' targetPairStrings{i}{2} ' p = ' num2str(p2)];
    else
        curText = [targetPairStrings{i}{1} ' = ' targetPairStrings{i}{2} ' p = ' num2str(p3)];
    end
    notBoxPlot([commonTargetTransVals{i,1} commonTargetTransVals{i,2}],[ones(1,length(commonTargetTransVals{i,1})) ones(1,length(commonTargetTransVals{i,2}))*2]);
    ylim([0 1])
    ylabel('Transition probability','fontsize',15,'fontweight','bold')
    set(gca,'xtick',[1 2],'xticklabels',{targetPairStrings{i}{1},targetPairStrings{i}{2}},'fontsize',15,'fontweight','bold')
    title(curText,'fontsize',15,'fontweight','bold')
end
suptitle('Common target')
set(gcf,'Position',[10 10 1600 1000])
if (do_save)
    figName = 'common_target';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
for i=1:length(commonSourcePairInds)
    subplot(1,3,i)
    [p1,h1] = signrank(earlyCommonSourceTransVals{i,1},earlyCommonSourceTransVals{i,2},'tail','right');
    [p2,h2] = signrank(earlyCommonSourceTransVals{i,1},earlyCommonSourceTransVals{i,2},'tail','left');
    [p3,h3] = signrank(earlyCommonSourceTransVals{i,1},earlyCommonSourceTransVals{i,2});
    p1 = p1*6; % Bonferroni correction
    p2 = p2*6; % Bonferroni correction
    if (p1 < .05)
        curText = [sourcePairStrings{i}{1} ' > ' sourcePairStrings{i}{2} ' p = ' num2str(p1)];
    elseif (p2 < .05)
        curText = [sourcePairStrings{i}{1} ' < ' sourcePairStrings{i}{2} ' p = ' num2str(p2)];
    else
        curText = [sourcePairStrings{i}{1} ' = ' sourcePairStrings{i}{2} ' p = ' num2str(p3)];
    end
    notBoxPlot([earlyCommonSourceTransVals{i,1} earlyCommonSourceTransVals{i,2}],[ones(1,length(earlyCommonSourceTransVals{i,1})) ones(1,length(earlyCommonSourceTransVals{i,2}))*2]);
    ylim([0 1])
    ylabel('Transition probability','fontsize',15,'fontweight','bold')
    set(gca,'xtick',[1 2],'xticklabels',{sourcePairStrings{i}{1},sourcePairStrings{i}{2}},'fontsize',15,'fontweight','bold')
    title(curText,'fontsize',15,'fontweight','bold')
end
suptitle('Early common source')
set(gcf,'Position',[10 10 1600 1000])
if (do_save)
    figName = 'common_source_early';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
for i=1:length(commonTargetPairInds)
    subplot(1,3,i)
    [p1,h1] = signrank(earlyCommonTargetTransVals{i,1},earlyCommonTargetTransVals{i,2},'tail','right');
    [p2,h2] = signrank(earlyCommonTargetTransVals{i,1},earlyCommonTargetTransVals{i,2},'tail','left');
    [p3,h3] = signrank(earlyCommonTargetTransVals{i,1},earlyCommonTargetTransVals{i,2});
    p1 = p1*6; % Bonferroni correction
    p2 = p2*6; % Bonferroni correction
    if (p1 < .05)
        curText = [targetPairStrings{i}{1} ' > ' targetPairStrings{i}{2} ' p = ' num2str(p1)];
    elseif (p2 < .05)
        curText = [targetPairStrings{i}{1} ' < ' targetPairStrings{i}{2} ' p = ' num2str(p2)];
    else
        curText = [targetPairStrings{i}{1} ' = ' targetPairStrings{i}{2} ' p = ' num2str(p3)];
    end
    notBoxPlot([earlyCommonTargetTransVals{i,1} earlyCommonTargetTransVals{i,2}],[ones(1,length(earlyCommonTargetTransVals{i,1})) ones(1,length(earlyCommonTargetTransVals{i,2}))*2]);
    ylim([0 1])
    ylabel('Transition probability','fontsize',15,'fontweight','bold')
    set(gca,'xtick',[1 2],'xticklabels',{targetPairStrings{i}{1},targetPairStrings{i}{2}},'fontsize',15,'fontweight','bold')
    title(curText,'fontsize',15,'fontweight','bold')
end
suptitle('Early common target')
set(gcf,'Position',[10 10 1600 1000])
if (do_save)
    figName = 'common_target_early';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
for i=1:length(commonSourcePairInds)
    subplot(1,3,i)
    [p1,h1] = signrank(lateCommonSourceTransVals{i,1},lateCommonSourceTransVals{i,2},'tail','right');
    [p2,h2] = signrank(lateCommonSourceTransVals{i,1},lateCommonSourceTransVals{i,2},'tail','left');
    [p3,h3] = signrank(lateCommonSourceTransVals{i,1},lateCommonSourceTransVals{i,2});
    p1 = p1*6; % Bonferroni correction
    p2 = p2*6; % Bonferroni correction
    if (p1 < .05)
        curText = [sourcePairStrings{i}{1} ' > ' sourcePairStrings{i}{2} ' p = ' num2str(p1)];
    elseif (p2 < .05)
        curText = [sourcePairStrings{i}{1} ' < ' sourcePairStrings{i}{2} ' p = ' num2str(p2)];
    else
        curText = [sourcePairStrings{i}{1} ' = ' sourcePairStrings{i}{2} ' p = ' num2str(p3)];
    end
    notBoxPlot([lateCommonSourceTransVals{i,1} lateCommonSourceTransVals{i,2}],[ones(1,length(lateCommonSourceTransVals{i,1})) ones(1,length(lateCommonSourceTransVals{i,2}))*2]);
    ylim([0 1])
    ylabel('Transition probability','fontsize',15,'fontweight','bold')
    set(gca,'xtick',[1 2],'xticklabels',{sourcePairStrings{i}{1},sourcePairStrings{i}{2}},'fontsize',15,'fontweight','bold')
    title(curText,'fontsize',15,'fontweight','bold')
end
suptitle('Late common source')
set(gcf,'Position',[10 10 1600 1000])
if (do_save)
    figName = 'common_source_late';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
for i=1:length(commonTargetPairInds)
    subplot(1,3,i)
    [p1,h1] = signrank(lateCommonTargetTransVals{i,1},lateCommonTargetTransVals{i,2},'tail','right');
    [p2,h2] = signrank(lateCommonTargetTransVals{i,1},lateCommonTargetTransVals{i,2},'tail','left');
    [p3,h3] = signrank(lateCommonTargetTransVals{i,1},lateCommonTargetTransVals{i,2});
    p1 = p1*6; % Bonferroni correction
    p2 = p2*6; % Bonferroni correction
    if (p1 < .05)
        curText = [targetPairStrings{i}{1} ' > ' targetPairStrings{i}{2} ' p = ' num2str(p1)];
    elseif (p2 < .05)
        curText = [targetPairStrings{i}{1} ' < ' targetPairStrings{i}{2} ' p = ' num2str(p2)];
    else
        curText = [targetPairStrings{i}{1} ' = ' targetPairStrings{i}{2} ' p = ' num2str(p3)];
    end
    notBoxPlot([lateCommonTargetTransVals{i,1} lateCommonTargetTransVals{i,2}],[ones(1,length(lateCommonTargetTransVals{i,1})) ones(1,length(lateCommonTargetTransVals{i,2}))*2]);
    ylim([0 1])
    ylabel('Transition probability','fontsize',15,'fontweight','bold')
    set(gca,'xtick',[1 2],'xticklabels',{targetPairStrings{i}{1},targetPairStrings{i}{2}},'fontsize',15,'fontweight','bold')
    title(curText,'fontsize',15,'fontweight','bold')
end
suptitle('Late common target')
set(gcf,'Position',[10 10 1600 1000])
if (do_save)
    figName = 'common_target_late';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Look at # of licks and amount consumed at A or B on AB pair days
Alicks = [];
Blicks = [];
Aconsumed = [];
Bconsumed = [];
for i=1:length(animalObjs)
    for j=1:size(animalObjs(i).solutions)
        if (any(ismember(animalObjs(i).solutions(j,:),2)) & any(ismember(animalObjs(i).solutions(j,:),3)))
            Aind = find(animalObjs(i).solutions(j,:) == 2);
            Bind = find(animalObjs(i).solutions(j,:) == 3);
            if (Aind == Bind)
                error('Solution A and B should not have the same index')
            end
            Alicks = [Alicks animalObjs(i).nLicksByBoxSide(j,Aind)];
            Blicks = [Blicks animalObjs(i).nLicksByBoxSide(j,Bind)];
            Aconsumed = [Aconsumed animalObjs(i).amountsConsumed(j,Aind)];
            Bconsumed = [Bconsumed animalObjs(i).amountsConsumed(j,Bind)];
        end
    end
end

figure;
notBoxPlot([Alicks Blicks],[ones(1,length(Alicks)) ones(1,length(Blicks))*2])
set(gca,'xtick',[1 2],'xticklabels',{'.01M','.1M'},'fontsize',15,'fontweight','bold')
title('Licks','fontsize',15,'fontweight','bold')
[p,h] = ranksum(Alicks,Blicks,'tail','right')

figure;
notBoxPlot([Aconsumed Bconsumed],[ones(1,length(Aconsumed)) ones(1,length(Bconsumed))*2])
set(gca,'xtick',[1 2],'xticklabels',{'.01M','.1M'},'fontsize',15,'fontweight','bold')
title('Amount consumed','fontsize',15,'fontweight','bold')
[p,h] = ranksum(Aconsumed,Bconsumed,'tail','right')

%% Look at # of licks and amount consumed at A or C on AC pair days
Alicks = [];
Clicks = [];
Aconsumed = [];
Cconsumed = [];
for i=1:length(animalObjs)
    for j=1:size(animalObjs(i).solutions)
        if (any(ismember(animalObjs(i).solutions(j,:),2)) & any(ismember(animalObjs(i).solutions(j,:),4)))
            Aind = find(animalObjs(i).solutions(j,:) == 2);
            Cind = find(animalObjs(i).solutions(j,:) == 4);
            if (Aind == Cind)
                error('Solution A and B should not have the same index')
            end
            Alicks = [Alicks animalObjs(i).nLicksByBoxSide(j,Aind)];
            Clicks = [Clicks animalObjs(i).nLicksByBoxSide(j,Cind)];
            Aconsumed = [Aconsumed animalObjs(i).amountsConsumed(j,Aind)];
            Cconsumed = [Cconsumed animalObjs(i).amountsConsumed(j,Cind)];
        end
    end
end

figure;
notBoxPlot([Alicks Clicks],[ones(1,length(Alicks)) ones(1,length(Clicks))*2])
set(gca,'xtick',[1 2],'xticklabels',{'.01M','1M'},'fontsize',15,'fontweight','bold')
title('Licks','fontsize',15,'fontweight','bold')
[p,h] = ranksum(Alicks,Clicks,'tail','right')

figure;
notBoxPlot([Aconsumed Cconsumed],[ones(1,length(Aconsumed)) ones(1,length(Cconsumed))*2])
set(gca,'xtick',[1 2],'xticklabels',{'.01M','1M'},'fontsize',15,'fontweight','bold')
title('Amount consumed','fontsize',15,'fontweight','bold')
[p,h] = ranksum(Aconsumed,Cconsumed,'tail','right')
%% Look at # of licks and amount consumed at B or C on BC pair days
Blicks = [];
Clicks = [];
Bconsumed = [];
Cconsumed = [];
for i=1:length(animalObjs)
    for j=1:size(animalObjs(i).solutions)
        if (any(ismember(animalObjs(i).solutions(j,:),3)) & any(ismember(animalObjs(i).solutions(j,:),4)))
            Bind = find(animalObjs(i).solutions(j,:) == 3);
            Cind = find(animalObjs(i).solutions(j,:) == 4);
            if (Bind == Cind)
                error('Solution A and B should not have the same index')
            end
            Blicks = [Blicks animalObjs(i).nLicksByBoxSide(j,Bind)];
            Clicks = [Clicks animalObjs(i).nLicksByBoxSide(j,Cind)];
            Bconsumed = [Bconsumed animalObjs(i).amountsConsumed(j,Bind)];
            Cconsumed = [Cconsumed animalObjs(i).amountsConsumed(j,Cind)];
        end
    end
end

figure;
notBoxPlot([Blicks Clicks],[ones(1,length(Blicks)) ones(1,length(Clicks))*2])
set(gca,'xtick',[1 2],'xticklabels',{'.1M','1M'},'fontsize',15,'fontweight','bold')
title('Licks','fontsize',15,'fontweight','bold')
[p,h] = ranksum(Blicks,Clicks,'tail','right')

figure;
notBoxPlot([Bconsumed Cconsumed],[ones(1,length(Bconsumed)) ones(1,length(Cconsumed))*2])
set(gca,'xtick',[1 2],'xticklabels',{'.1M','1M'},'fontsize',15,'fontweight','bold')
title('Amount consumed','fontsize',15,'fontweight','bold')
[p,h] = ranksum(Bconsumed,Cconsumed,'tail','right')

%% Figure 22: # of bouts for each solution pair
figure;
xlabels = {'H2O/H2O','A/A','A/B','A/C','B/B','B/C','C/C'};
nBoutsByPair = cell(1,7);
for i=setdiff(1:length(animalObjs),5) % exclude bb12
    count = 0;
    for j=1:4
        for k=j:4
            if (j == 1 && k > 1)
                continue;
            end
            count = count+1;
            nBoutsByPair{count} = [nBoutsByPair{count} allNbouts{i}{j,k}];
        end
    end
end
for i=1:length(nBoutsByPair)
    scatter(i*ones(1,length(nBoutsByPair{i})),nBoutsByPair{i},100,'k.')
    hold on;
    plot([i-.2 i+.2],[mean(nBoutsByPair{i}) mean(nBoutsByPair{i})],'r')
    plot([i-.2 i+.2],[median(nBoutsByPair{i}) median(nBoutsByPair{i})],'b')
end
set(gca,'xtick',1:7,'xticklabels',xlabels,'xticklabelrotation',20,'fontsize',15,'fontweight','bold');
xlim([0 8])
ylabel('# of bouts','fontsize',15,'fontweight','bold')
hold off;
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'number_of_bouts_by_solution';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 23: # of switches for each solution pair
figure;
xlabels = {'H2O/H2O','A/A','A/B','A/C','B/B','B/C','C/C'};
nSwitchVals = cell(1,7);
for i=setdiff(1:length(animalObjs),5) % exclude bb12
    count = 0;
    for j=1:4
        for k=j:4
            if (j == 1 && k > 1)
                continue;
            end
            count = count+1;
            nSwitchVals{count} = [nSwitchVals{count} allSolnNswitches{i}{j,k}];
        end
    end
end
for i=1:length(nSwitchVals)
    scatter(i*ones(1,length(nSwitchVals{i})),nSwitchVals{i},100,'k.')
    hold on;
    plot([i-.2 i+.2],[mean(nSwitchVals{i}) mean(nSwitchVals{i})],'r')
    plot([i-.2 i+.2],[median(nSwitchVals{i}) median(nSwitchVals{i})],'b')
end
set(gca,'xtick',1:7,'xticklabels',xlabels,'xticklabelrotation',20); xlim([0 8])
ylabel('# of switches','fontsize',15,'fontweight','bold')
hold off;
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'number_of_switches_by_solution_pair';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

linNSwitchVals = [];
linXcoords = [];
for i=1:length(nSwitchVals)
    linNSwitchVals = [linNSwitchVals nSwitchVals{i}];
    linXcoords = [linXcoords i*ones(1,length(nSwitchVals{i}))];
end
figure;
notBoxPlot(linNSwitchVals,linXcoords)
ylabel('# of switches','fontsize',15,'fontweight','bold')
set(gca,'xtick',1:7,'xticklabels',{'H2O/H2O','.01M/.01M','.01M/.1M','.01M/1M','.1M/.1M','.1M/1M','1M/1M'},...
    'xticklabelrotation',45,'FontSize',20)
if (do_save)
    figName = 'number_of_switches_by_solution_pair_boxplot';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 24: probability of choosing preferred side from day before first by animal
nPriorPreferred = zeros(1,length(setdiff(1:length(animalObjs),5)));
nPriorPreferredFrac = zeros(1,length(setdiff(1:length(animalObjs),5)));
count = 0;
totalBinomTestCount = 0;
for i=setdiff(1:length(animalObjs),5) % exclude bb12
    count = count + 1;
    for j=2:length(animalObjs(i).dates)
        totalBinomTestCount = totalBinomTestCount + 1;
        disp(['Fig 24: ' animalObjs(i).name ' ' animalObjs(i).dates{j}])
        licks = load(['analyzed_data/' animalObjs(i).dates{j-1} '/' animalObjs(i).name '/licks.mat']); licks=licks.licks;
        if (isempty(licks{1}))
            if (strcmp(licks{2}(1).box_side,'left'))
                nLicksLeft = length(licks{2});
                nLicksRight = 0;
            elseif (strcmp(licks{2}(1).box_side,'right'))
                nLicksRight = length(licks{2});
                nLicksLeft = 0;
            else
                error('neither channel had any licks')
            end
        elseif (isempty(licks{2}))
            if (strcmp(licks{1}(1).box_side,'left'))
                nLicksLeft = length(licks{1});
                nLicksRight = 0;
            elseif (strcmp(licks{1}(1).box_side,'right'))
                nLicksRight = length(licks{1});
                nLicksLeft = 0;
            end
        else
            if (strcmp(licks{1}(1).box_side,'left'))
                nLicksLeft = length(licks{1});
                nLicksright = length(licks{2});
            else
                nLicksRight = length(licks{1});
                nLicksLeft = length(licks{2});
            end
        end
        if (nLicksLeft > nLicksRight)
            priorPreferred = 'left';
        else
            priorPreferred = 'right';
        end
        firstBout = animalObjs(i).linBoutsByDay{j}(1);
        if (strcmp(firstBout.box_side,priorPreferred))
            nPriorPreferred(count) = nPriorPreferred(count) + 1;
        end 
    end
    nPriorPreferredFrac(count) = nPriorPreferred(count)/(length(animalObjs(i).dates)-1);
end
binomPDF = binopdf(0:1:totalBinomTestCount,totalBinomTestCount,.5);
figure;
histogram(nPriorPreferredFrac,'binwidth',.05)
xlabel('Probability of choosing prior preferred side','fontsize',15,'fontweight','bold')
ylabel('# of animals','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'probability_of_choosing_prior_day_preferred_side';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
plot(0:totalBinomTestCount,binomPDF); hold on; plot([sum(nPriorPreferred) sum(nPriorPreferred)],[0 .1],'r')
xlabel('# of prior preferred chosen','fontsize',15,'fontweight','bold');
ylabel('Probability density','fontsize',15,'fontweight','bold')
if (do_save)
    figName = 'prior_day_preferred_binomial_test';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end 

%% Figure 25: Time to first lick
xlabels = {};
labelcount = 0;
for i=1:length(animalObjs)
    if (~strcmp(animalObjs(i).name,'bb12'))
        labelcount = labelcount + 1;
        xlabels{labelcount} = animalObjs(i).name;
    end
end
count = 0;
time2firstLick = cell(1,length(setdiff(1:length(animalObjs),5)));
allFirstLickTimes = [];
for i=setdiff(1:length(animalObjs),5) % exclude bb12
    count = count+1;
    for j=1:length(animalObjs(i).dates)
        time2firstLick{count} = [time2firstLick{count} animalObjs(i).linBoutsByDay{j}(1).onset - 300];
        allFirstLickTimes = [allFirstLickTimes animalObjs(i).linBoutsByDay{j}(1).onset - 300];
    end
end
figure;
for i=1:length(time2firstLick)
    scatter(i*ones(1,length(time2firstLick{i})),time2firstLick{i},100,'k.')
    hold on;
    plot([i-.2 i+.2],[mean(time2firstLick{i}) mean(time2firstLick{i})],'r')
    plot([i-.2 i+.2],[median(time2firstLick{i}) median(time2firstLick{i})],'b')
end
set(gca,'xtick',1:length(setdiff(1:length(animalObjs),5)),'xticklabels',xlabels,'xticklabelrotation',20,'fontsize',15,'fontweight','bold'); 
xlim([0 length(xlabels)+1])
ylabel('Time to first lick (s)','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'time_to_first_lick_by_animal';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 26: Distribution of all times to first lick
figure;
histogram(allFirstLickTimes,'binwidth',2)
xlabel('Time to first lick (s)','fontsize',15,'fontweight','bold')
ylabel('# of occurrences','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'time_to_first_lick_distribution';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end


%% Figure 28: # of switches by solution pair 1st vs. 2nd session
allFirstVisitInds = cell(1,length(animalObjs));
allSecondVisitInds = cell(1,length(animalObjs));
for i=1:length(animalObjs)
    firstVisitInds = nan(1,7);
    secondVisitInds = nan(1,7);
    for j=1:length(animalObjs(i).dates)
        if (animalObjs(i).solutions(j,1) == 1 && animalObjs(i).solutions(j,2) == 1)
            solnPairID = 1;
        elseif (any(ismember(animalObjs(i).solutions(j,:),2)) && any(ismember(animalObjs(i).solutions(j,:),3)))
            solnPairID = 2;
        elseif (any(ismember(animalObjs(i).solutions(j,:),2)) && any(ismember(animalObjs(i).solutions(j,:),4)))
            solnPairID = 3;
        elseif (any(ismember(animalObjs(i).solutions(j,:),3)) && any(ismember(animalObjs(i).solutions(j,:),4)))
            solnPairID = 4;
        elseif (animalObjs(i).solutions(j,1) == 3 && animalObjs(i).solutions(j,2) == 3)
            solnPairID = 5;
        elseif (animalObjs(i).solutions(j,1) == 4 && animalObjs(i).solutions(j,2) == 4)
            solnPairID = 6;
        elseif (animalObjs(i).solutions(j,1) == 2 && animalObjs(i).solutions(j,2) == 2)
            solnPairID = 7;
        else
            error('Cannot identify solution pair')
        end
        if (isnan(firstVisitInds(solnPairID)))
            firstVisitInds(solnPairID) = j;
        else
            secondVisitInds(solnPairID) = j;
        end
    end
    allFirstVisitInds{i} = firstVisitInds;
    allSecondVisitInds{i} = secondVisitInds;
end

firstVisitSwitches = zeros(length(animalObjs),7);
secondVisitSwitches = zeros(length(animalObjs),7);
for i=1:length(animalObjs)
    for j=1:7
        if (~isnan(allFirstVisitInds{i}(j)))
            [T,transitionCount,totalTransitions,nswitches,nbouts] = getBoutTransitionMatrix(animalObjs(i).bouts{allFirstVisitInds{i}(j)},animalObjs(i).solutions(allFirstVisitInds{i}(j),:));
            firstVisitSwitches(i,j) = nswitches;
        else
            firstVisitSwitches(i,j) = nan;
        end
        if (~isnan(allSecondVisitInds{i}(j)))
            [T,transitionCount,totalTransitions,nswitches,nbouts] = getBoutTransitionMatrix(animalObjs(i).bouts{allSecondVisitInds{i}(j)},animalObjs(i).solutions(allSecondVisitInds{i}(j),:));
            secondVisitSwitches(i,j) = nswitches;
        else
            secondVisitSwitches(i,j) = nan;
        end
    end
end
meanFirstVisitSwitches = mean(firstVisitSwitches,1,'omitnan');
stdFirstVisitSwitches = std(firstVisitSwitches,[],1,'omitnan');
meanSecondVisitSwitches = mean(secondVisitSwitches,1,'omitnan');
stdSecondVisitSwitches = std(secondVisitSwitches,[],1,'omitnan');
%xlabels = {'H2O/H2O','A/A','A/B','A/C','B/B','B/C','C/C'};
xlabels = {'H2O/H2O','A/B','A/C','B/C','B/B','C/C','A/A'};
figure;
shadedErrorBar(1:7,meanFirstVisitSwitches,stdFirstVisitSwitches,'lineprops','-r')
hold on;
shadedErrorBar(1:7,meanSecondVisitSwitches,stdSecondVisitSwitches,'lineprops','-b')
set(gca,'xtick',1:7,'xticklabels',xlabels,'fontsize',15,'fontweight','bold')
ylabel('# of switches','fontsize',15,'fontweight','bold')
if (do_save)
    figName = 'nSwitches_by_solutionPair_1st_vs_2nd_exposure';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Switch frequency vs. time. Aggregate and by solution pair
solnPairs = [2 3; 2 4; 3 4];
switchTimesPrevBoutEnd = [];
switchTimesNextBoutStart = [];
switchTimesMidpoints = [];
solnPairInds = [];
switchInterBoutIntervals = [];
difSolInds = find(boutTable.solutionNum ~= boutTable.alternative_solutionNum);
for i=1:length(animalObjs)
    ratInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    for j=animalObjs(i).difSolutionDays
        dayInds = find(strcmp(boutTable.date,animalObjs(i).dates{j}));
        curInds = intersect(intersect(difSolInds,ratInds),dayInds);
        [~,sortedOrder] = sort(boutTable.boutNumber(curInds));
        sortedInds = curInds(sortedOrder);
        curSolns = unique([boutTable.solutionNum(sortedInds); boutTable.alternative_solutionNum(sortedInds)]);
        if (any(curSolns == solnPairs(1,1)) & any(curSolns == solnPairs(1,2)))
            curSolnPairInd = 1;
        elseif (any(curSolns == solnPairs(2,1)) & any(curSolns == solnPairs(2,2)))
            curSolnPairInd = 2;
        elseif (any(curSolns == solnPairs(3,1)) & any(curSolns == solnPairs(3,2)))
            curSolnPairInd = 3;
        else
            error('Solution pair not recognized')
        end
        if (length(sortedInds) > 1)
            for k=2:length(sortedInds)
                if (~strcmp(boutTable.box_side(sortedInds(k)),boutTable.box_side(sortedInds(k-1))))
                    switchTimesPrevBoutEnd = [switchTimesPrevBoutEnd boutTable.offset(sortedInds(k-1))];
                    switchTimesNextBoutStart = [switchTimesNextBoutStart boutTable.onset(sortedInds(k))];
                    switchTimesMidpoints = [switchTimesMidpoints (boutTable.offset(sortedInds(k-1)) + boutTable.onset(sortedInds(k)))/2];
                    solnPairInds = [solnPairInds curSolnPairInd];
                    switchInterBoutIntervals = [switchInterBoutIntervals boutTable.onset(sortedInds(k)) - boutTable.offset(sortedInds(k-1))];
                end
            end
        end
    end
end

figure;
subplot(1,2,1)
histogram(switchTimesMidpoints)
title('Switch time midpoints','fontsize',15,'fontweight','bold')
subplot(1,2,2)
histogram(switchInterBoutIntervals)
title('Switch inter-bout-intervals','fontsize',15,'fontweight','bold')

figure;
[rho,pval] = corr(switchTimesMidpoints',switchInterBoutIntervals');
maxX = max(switchTimesMidpoints);
maxY = max(switchInterBoutIntervals);
scatter(switchTimesMidpoints,switchInterBoutIntervals,100,'k.')
text(.1*maxX,.9*maxY,['\rho = ' num2str(rho)])
xlabel('Switch time midpoint','fontsize',15,'fontweight','bold')
ylabel('Inter-bout-interval','fontsize',15,'fontweight','bold')

figure;
for i=1:size(solnPairs,1)
    subplot(2,2,i)
    inds = find(solnPairInds == i);
    histogram(switchTimesMidpoints(inds),'binwidth',100,'normalization','pdf')
    title(['Solution pair ' num2str(i)],'fontsize',15,'fontweight','bold')
end

%% Estimate 'desired' molarity as weighted average of solution consumed
count = 0;
weightedMolarity = [];
for i=1:length(animalObjs)
    for j=animalObjs(i).difSolutionDays
        count = count+1;
        for k=1:2
            solNum(k) = animalObjs(i).solutions(j,k);
            molarity(k) = .01*10^(solNum(k)-2);
            consumed(k) = animalObjs(i).amountsConsumed(j,k);
        end
        weightedMolarity(count) = sum(molarity.*consumed)/sum(consumed);
    end
end

figure;
histogram(weightedMolarity,'binwidth',.01,'normalization','pdf')
xlabel('Consumption weighted molarity','fontsize',15,'fontweight','bold')
ylabel('Probability density','fontsize',15,'fontweight','bold')

%% Estimate 'desired' molarity as weighted average of solution consumed for solution pairs
solnPairs = [2 3; 2 4; 3 4];
weightedMolarityBySolution = cell(1,size(solnPairs,1));

for i=1:length(animalObjs)
    for j=animalObjs(i).difSolutionDays
        if (any(animalObjs(i).solutions(j,:) == solnPairs(1,1)) & any(animalObjs(i).solutions(j,:) == solnPairs(1,2)))
            solnPairInd = 1;
        elseif (any(animalObjs(i).solutions(j,:) == solnPairs(2,1)) & any(animalObjs(i).solutions(j,:) == solnPairs(2,2)))
            solnPairInd = 2;
        elseif (any(animalObjs(i).solutions(j,:) == solnPairs(3,1)) & any(animalObjs(i).solutions(j,:) == solnPairs(3,2)))
            solnPairInd = 3;
        else
            disp(animalObjs(i).solutions(j,:))
            error('Solution pair not recognized')
        end
        for k=1:2
            solNum(k) = animalObjs(i).solutions(j,k);
            molarity(k) = .01*10^(solNum(k)-2);
            consumed(k) = animalObjs(i).amountsConsumed(j,k);
        end
        weightedMolarityBySolution{solnPairInd} = [weightedMolarityBySolution{solnPairInd} sum(molarity.*consumed)/sum(consumed)];
    end
end

xticklabels = {'A/B','A/C','B/C'};
figure;
for i=1:size(solnPairs,1)
    notBoxPlot(weightedMolarityBySolution{i},ones(1,length(weightedMolarityBySolution{i}))*i)
end
set(gca,'xtick',[1 2 3],'xticklabels',xticklabels,'fontsize',15,'fontweight','bold')
ylabel('Consumption weighted molarity','fontsize',15,'fontweight','bold')

%% Figure 30: Fitting different distributions to bout durations
doubleExponential = @(x,w,lambda1,lambda2) w*expdf(x,lambda1) + (1-w)*exppdf(x,lambda2);
doubleNormal = @(x,w,mu1,mu2,std1,std2) w*normpdf(x,mu1,std1) + (1-w)*normpdf(x,mu2,std2);
distributionNames = {'exponential','poisson','normal','lognormal','gamma','doubleExponential','doubleNormal'};
distributionFuncs = {@exppdf,@poisspdf,@normpdf,@lognpdf,@gampdf,@doubleExponential,@doubleNormal};
bestFitParams = cell(length(animalObjs),4,length(distributionNames));
BIC_fits = zeros(length(animalObjs),4,length(distributionNames));
AIC_fits = zeros(length(animalObjs),4,length(distributionNames));
neglogliks = nan(length(animalObjs),4,length(distributionNames));
bestBICInds = nan(length(animalObjs),4);
bestAICInds = nan(length(animalObjs),4);
x = 0:1:600;
for i=1:length(animalObjs)
    for j=1:4
        curDurations = [animalObjs(i).boutsBySolution{j}.duration];
        curDurationsLength = length(curDurations);
        bestBIC = inf;
        bestAIC = inf;
        for k=1:length(distributionNames)
            if (k <= 5)
                pd = fitdist(curDurations',distributionNames{k});
                if (k == 2)
                    disp(num2str(pd.negloglik))
                end
                bestFitParams{i,j,k} = pd.ParameterValues;
                neglogliks(i,j,k) = pd.negloglik;
                bic = pd.NumParameters*log(length(pd.InputData.data)) + 2*pd.negloglik;
                aic = 2*pd.NumParameters + 2*pd.negloglik;
            elseif (k == 6)
                lambdaStart = bestFitParams{i,j,1};
                startVec = [.9 lambdaStart 10*lambdaStart];
                startVec = sort(startVec,'descend');
                try
                    [phat,pic] = mle(curDurations,'pdf',doubleExponential,'start',startVec);
                    negloglik = -sum(log(doubleExponential(curDurations,phat(1),phat(2),phat(3))));
                    bestFitParams{i,j,k} = phat;
                    neglogliks(i,j,k) = negloglik;
                    bic = 3*log(curDurationsLength) + 2*negloglik;
                    aic = 2*3 + 2*negloglik;
                catch
                    bestFitParams{i,j,k} = nan(1,length(startVec));
                    neglogliks(i,j,k) = nan;
                    bic = nan;
                    aic = nan;
                end
            elseif (k == 7)
                muStart = bestFitParams{i,j,3}(1);
                stdStart = bestFitParams{i,j,3}(2);
                muStartVec = [muStart 10*muStart];
                stdStartVec = [.5*stdStart 5*stdStart];
                startVec = [.9 muStartVec stdStartVec];
                try
                    [phat,pic] = mle(curDurations,'pdf',doubleNormal,'start',startVec);
                    negloglik = -sum(log(doubleNormal(curDurations,phat(1),phat(2),phat(3),phat(4),phat(5))));
                    bestFitParams{i,j,k} = phat;
                    neglogliks(i,j,k) = negloglik;
                    bic = 5*log(curDurationsLength) + 2*negloglik;
                    aic = 2*5 + 2*negloglik;
                catch
                    bestFitParams{i,j,k} = nan(1,length(startVec));
                    neglogliks(i,j,k) = nan;
                    bic = nan;
                    aic = nan;
                end
            end
            
            BIC_fits(i,j,k) = bic;
            AIC_fits(i,j,k) = aic;
            if (bic < bestBIC)
                bestBIC = bic;
                bestBICInds(i,j) = k;
            end
            if (aic < bestAIC)
                bestAIC = aic;
                bestAICInds(i,j) = k;
            end
        end
    end
end

%% Figure 31
% Fit distributions to ILI distributions
%doubleNormal = @(x,w,mu1,mu2,std1,std2) w*normpdf(x,mu1,std1) + (1-w)*normpdf(x,mu2,std2);
doubleNormal = @(x,w,mu1,mu2,std1,std2) max(1e-100,(w/(std1*sqrt(2*pi)))*exp(-.5*((x - mu1)./std1).^2) + ((1-w)/(std2*sqrt(2*pi)))*exp(-.5*((x - mu2)./std2).^2));
doubleNormalPhat = cell(length(animalObjs),4);
expPhat = cell(4,length(animalObjs));
iliNegLogLiks = nan(length(animalObjs),4,2);
allILIsBySolution = cell(1,4);
for i=1:length(animalObjs)
    figure;
    for j=1:4
        if (any(animalObjs(i).ilisBySolution{j} < 0))
            disp([num2str(i) ' ' num2str(j) ' has negative ilis'])
            badinds = find(animalObjs(i).ilisBySolution{j} < 0);
            goodILIs = animalObjs(i).ilisBySolution{j}(setdiff(1:length(animalObjs(i).ilisBySolution{j}),badinds));
            disp(['num bad inds ' num2str(length(badinds))])
        else
            goodILIs = animalObjs(i).ilisBySolution{j};
        end
        allILIsBySolution{j} = [allILIsBySolution{j} goodILIs(goodILIs <= 2)];
        if (length(goodILIs) < 100)
            disp([num2str(i) ' ' num2str(j)])
        end
        curMean = mean(goodILIs);
        curStd = std(goodILIs);
        pd = fitdist(goodILIs','exponential');
        expPhat{j,i} = pd.ParameterValues;
        iliNegLogLiks(i,j,1) = pd.negloglik;
        isDone = 0;
        tryNumber = 0;
        bestNegLogLik = Inf;
        bestParams = [];
        for k=1:10
            lb = [.001 .001 .001 .001 .001];
            ub = [.999 2 2 2 2];
            startVec = [(rand*(ub(1)-lb(1)) + lb(1)) (rand*(ub(2)-lb(2)) + lb(2)) (rand*(ub(3)-lb(3)) + lb(3)) (rand*(ub(4)-lb(4)) + lb(4)) (rand*(ub(5)-lb(5)) + lb(5))];
            [phat,pic] = mle(goodILIs,'pdf',doubleNormal,'start',startVec,'LowerBound',lb,'UpperBound',ub,'optimfun','fmincon');
            negloglik = -sum(log(doubleNormal(goodILIs,phat(1),phat(2),phat(3),phat(4),phat(5))));
            if (negloglik < bestNegLogLik)
                bestNegLogLik = negloglik;
                bestParams = phat;
            end
        end
        iliNegLogLiks(i,j,2) = bestNegLogLik;
        doubleNormalPhat{i,j} = bestParams;
        subplot(2,2,j)
        histogram(goodILIs,'normalization','pdf','binwidth',.005);
        hold on;
        xvec = 0:.001:max(goodILIs);
        plot(xvec,doubleNormal(xvec,bestParams(1),bestParams(2),bestParams(3),bestParams(4),bestParams(5)))
        xlabel('ILI (s)','fontsize',15,'fontweight','bold')
        ylabel('Probability Density','fontsize',15,'fontweight','bold')
    end
    disp(['Done fitting ILI distributions for animal ' num2str(i)])
end

%% Fit ILI distributions
bestParams = cell(1,4);
for i=1:4
    bestNegLogLik = Inf;
    lb = [.001 .001 .001 .001 .001];
    ub = [.999 2 2 2 2];
    for j=1:50
        startVec = [(rand*(ub(1)-lb(1)) + lb(1)) (rand*(ub(2)-lb(2)) + lb(2)) (rand*(ub(3)-lb(3)) + lb(3)) (rand*(ub(4)-lb(4)) + lb(4)) (rand*(ub(5)-lb(5)) + lb(5))];
        [phat,pic] = mle(goodILIs,'pdf',doubleNormal,'start',startVec,'LowerBound',lb,'UpperBound',ub,'optimfun','fmincon');
        negloglik = -sum(log(doubleNormal(goodILIs,phat(1),phat(2),phat(3),phat(4),phat(5))));
        if (negloglik < bestNegLogLik)
            bestNegLogLik = negloglik;
            bestParams{i} = phat;
        end
    end
end
figure;
for i=1:4
    subplot(2,2,i)
    histogram(allILIsBySolution{i},'normalization','pdf','binwidth',.001)
    hold on;
    xvec = 0:.001:2;
    plot(xvec,doubleNormal(xvec,bestParams{i}(1),bestParams{i}(2),bestParams{i}(3),bestParams{i}(4),bestParams{i}(5)))
    xlim([0 2])
end
suptitle('All animal ILIs combined')

