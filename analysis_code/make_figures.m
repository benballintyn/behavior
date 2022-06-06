%% make_figures
clear all; close all;
figFolder = '/home/ben/phd/behavior/figures/03_08_21/';
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


%% FIGURE 1
% cdf plot of lick times
figure;
for i=1:length(animalObjs)
    cdfplot([animalObjs(i).linLicks.onset])
    hold on;
end
xlabel('Time (s)'); ylabel('P(t_{lick} < x)')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    saveas(gcf,[figFolder 'lick_times_cdf'],'fig')
    saveas(gcf,[figFolder 'lick_times_cdf'],'eps')
    saveas(gcf,[figFolder 'lick_times_cdf'],'png')
end

%% Figure 2
% cdf plot of bout onset times
figure;
for i=1:length(animalObjs)
    cdfplot([animalObjs(i).linBouts.onset])
    hold on;
end
xlabel('Time (s)'); ylabel('P(t_{bout} < x)')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    saveas(gcf,[figFolder 'bout_onset_times_cdf'],'fig')
    saveas(gcf,[figFolder 'bout_onset_times_cdf'],'eps')
    saveas(gcf,[figFolder 'bout_onset_times_cdf'],'png')
end

%% Figure 3
% bout lengths vs time
figure;
for i=1:length(animalObjs)
    scatter([animalObjs(i).linBouts.onset],[animalObjs(i).linBouts.duration],100,'k.')
    hold on;
end
xlabel('Time(s)'); ylabel('Bout duration (s)')
%x = 300:1:3900;
%plot(x,f(B,x-300))
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    saveas(gcf,[figFolder 'bout_lengths_vs_time'],'fig')
    saveas(gcf,[figFolder 'bout_lengths_vs_time'],'eps')
    saveas(gcf,[figFolder 'bout_lengths_vs_time'],'png')
end

%% Figure 4
% Bout duration as a function of total prior licks
allBoutDurations = [];
allBoutNlicks = [];
prevLicks = [];
prevDurations = [];
for i=1:length(animalObjs)
    for j=1:length(animalObjs(i).dates)
        curDayBoutOnsets = [];
        curDayBoutOffsets = [];
        curDayBoutNlicks = [];
        curDayBoutDurations = [];
        curDayPrevLicks = [];
        curDayPrevDurations = [];
        for k=1:length(animalObjs(i).bouts{j})
            if (length(animalObjs(i).bouts{j}{k}) > 0)
                curDayBoutOnsets = [curDayBoutOnsets [animalObjs(i).bouts{j}{k}.onset]];
                curDayBoutOffsets = [curDayBoutOffsets [animalObjs(i).bouts{j}{k}.offset]];
                curDayBoutNlicks = [curDayBoutNlicks [animalObjs(i).bouts{j}{k}.nlicks]];
                curDayBoutDurations = [curDayBoutDurations [animalObjs(i).bouts{j}{k}.duration]];
            end
        end
        for k=1:length(curDayBoutDurations)
            prevInds = find(curDayBoutOffsets < curDayBoutOnsets(k));
            curDayPrevLicks = [curDayPrevLicks sum(curDayBoutNlicks(prevInds))];
            curDayPrevDurations = [curDayPrevDurations sum(curDayBoutDurations(prevInds))];
        end
        allBoutDurations = [allBoutDurations curDayBoutDurations];
        allBoutNlicks = [allBoutNlicks curDayBoutNlicks];
        prevLicks = [prevLicks curDayPrevLicks];
        prevDurations = [prevDurations curDayPrevDurations];
    end
end
[rho,pval] = corr(prevLicks',allBoutDurations');
exponentialFit = fit(prevLicks',allBoutDurations','exp1');
posPrevLicks = find(prevLicks > 0);
noPrevLicks = find(prevLicks == 0);
powerFit = fit(prevLicks(posPrevLicks)',allBoutDurations(posPrevLicks)','power1');
adjustedPrevLicks = prevLicks;
adjustedPrevLicks(noPrevLicks) = 1;
fullPowerFit = fit(adjustedPrevLicks',allBoutDurations','power1');
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
xlabel('Total prior licks')
ylabel('Bout duration')
text(3000,400,{['\rho = ' num2str(rho)],['p = ' num2str(pval)],['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],...
    ['R^2 (lin): ' num2str(rsquared_linear)],...
    ['R^2 (exp): ' num2str(rsquared_exp)],...
    ['R^2 (pow): ' num2str(rsquared_power)],...
    ['R^2 (pow full): ' num2str(rsquared_power_full)]})
legend({'Data','Linear Fit','Exponential Fit','Power law fit','full Power law fit'})

[rho2,pval2]=corr(prevDurations',allBoutDurations');
lf2 = polyfit(prevDurations,allBoutDurations,1);
subplot(1,2,2)
scatter(prevDurations,allBoutDurations,100,'.')
hold on;
xvals = 0:1:round(max(prevDurations));
plot(xvals,lf2(1)*xvals + lf(2))
xlabel('Total prior bout time')
ylabel('Bout duration')
text(max(xvals)/2,400,{['\rho = ' num2str(rho2)],['p = ' num2str(pval2)],['y = ' num2str(lf2(1)) 'x + ' num2str(lf2(2))]})
set(gcf,'Position',[10 10 1400 1200])
if (do_save)
    saveas(gcf,[figFolder 'bout_duration_vs_prior_licks'],'fig')
    saveas(gcf,[figFolder 'bout_duration_vs_prior_licks'],'eps')
    saveas(gcf,[figFolder 'bout_duration_vs_prior_licks'],'png')
end

%% Figure 5
% # of licks vs. bout duration
figure;
scatter(allBoutNlicks,allBoutDurations,100,'.')
xlabel('# of licks in bout')
ylabel('Bout duration')
if (do_save)
    saveas(gcf,[figFolder 'nlicks_vs_bout_duration'],'fig')
    saveas(gcf,[figFolder 'nlicks_vs_bout_duration'],'eps')
    saveas(gcf,[figFolder 'nlicks_vs_bout_duration'],'png')
end

%% Figure 6
% Average lick rate by solution
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
ylabel('Lick rate')
if (do_save)
    figName = 'lickRateBySolution';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    saveas(gcf,[figFolder figName],'png')
end

%% Figure 7
% bout duration as a function of previous bout duration
figure;
xs = [];
ys = [];
for i=1:length(animalObjs)
    for j=1:length(animalObjs(i).dates)
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
xlabel('Current bout duration (s)')
ylabel('Next bout duration (s)')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    saveas(gcf,[figFolder 'bout_duration_vs_prev_bout_duration'],'fig')
    saveas(gcf,[figFolder 'bout_duration_vs_prev_bout_duration'],'eps')
    saveas(gcf,[figFolder 'bout_duration_vs_prev_bout_duration'],'png')
end

%% Figure 8
% total consumed for all animals
xlabels = {'H2O','.01M','.1M','1M'};
figure;
for i=1:length(animalObjs)
    totalConsumed = animalObjs(i).totalConsumed;
    scales = animalObjs(i).totalBottles./animalObjs(i).totalBottles(1);
    scaledConsumed = totalConsumed./scales;
    normalizedScaledConsumed = scaledConsumed./scaledConsumed(1);
    plot(1:4,normalizedScaledConsumed)
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',xlabels)
ylabel('Relative palatability')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    saveas(gcf,[figFolder 'total_consumed_all_animals'],'fig')
    saveas(gcf,[figFolder 'total_consumed_all_animals'],'eps')
    saveas(gcf,[figFolder 'total_consumed_all_animals'],'png')
end

%% Figure 9
% total licks for all animals
figure;
for i=1:length(animalObjs)
    totalLicks = animalObjs(i).totalLicks;
    scales = animalObjs(i).totalBottles./animalObjs(i).totalBottles(1);
    scaledLicks = totalLicks./scales;
    normalizedScaledLicks = scaledLicks./scaledLicks(1);
    allNormalizedScaledLicks(i,:) = normalizedScaledLicks;
    plot(1:4,normalizedScaledLicks)
    hold on;
end
set(gca,'xtick',1:4,'xticklabels',xlabels)
ylabel('Relative palatability')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    saveas(gcf,[figFolder 'total_licks_all_animals'],'fig')
    saveas(gcf,[figFolder 'total_licks_all_animals'],'eps')
    saveas(gcf,[figFolder 'total_licks_all_animals'],'png')
end

%% Figure 10
% relative palatability by amount consumed or total licks
figure;
for i=1:length(animalObjs)
    rpalLicks(:,i) = animalObjs(i).relativePalatabilitiesLicks;
    rpalConsumed(:,i) = animalObjs(i).relativePalatabilitiesConsumed;
end
subplot(1,2,1); plot(1:4,rpalConsumed); set(gca,'xticklabels',xlabels); ylabel('Relative Palatabilitiy'); title('Relative Palatability based on consumption')
subplot(1,2,2); plot(1:4,rpalLicks); set(gca,'xticklabels',xlabels); ylabel('RelativePalatability'); title('Relative Palatability based on lick totals')
set(gcf,'Position',[10 10 1800 1000])
if (do_save)
    saveas(gcf,[figFolder 'relative_palatability_by_amount_consumed_or_licks'],'fig')
    saveas(gcf,[figFolder 'relative_palatability_by_amount_consumed_or_licks'],'eps')
    saveas(gcf,[figFolder 'relative_palatability_by_amount_consumed_or_licks'],'png')
end

%% Figure 11
% Relative palatability by solution and animal sex
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
set(gca,'xtick',2:2:8,'xticklabels',xlabels)
ylabel('Relative palatabilitiy')
if (do_save)
    figName = 'relativePalatabilities_by_solution_and_sex';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    saveas(gcf,[figFolder figName],'png')
end

%% Figure 12
% Bout duration distributions for each solution
figure;
allBoutsBySolution = cell(1,4);
maxDur = 0;
for i=1:length(animalObjs)
    for j=1:4
        allBoutsBySolution{j} = [allBoutsBySolution{j} [animalObjs(i).boutsBySolution{j}.duration]];
        maxDur = max(maxDur,max(allBoutsBySolution{j}));
    end
end
subplot(2,2,1); histogram(allBoutsBySolution{1},'Normalization','pdf','binwidth',1); xlabel('Bout duration (s)'); title('H2O'); xlim([0 maxDur]); ylim([0 .2])
subplot(2,2,2); histogram(allBoutsBySolution{2},'Normalization','pdf','binwidth',1); xlabel('Bout duration (s)'); title('.01M NaCl'); xlim([0 maxDur]); ylim([0 .2])
subplot(2,2,3); histogram(allBoutsBySolution{3},'Normalization','pdf','binwidth',1); xlabel('Bout duration (s)'); title('.1M NaCl'); xlim([0 maxDur]); ylim([0 .2])
subplot(2,2,4); histogram(allBoutsBySolution{4},'Normalization','pdf','binwidth',1); xlabel('Bout duration (s)'); title('1M NaCl'); xlim([0 maxDur]); ylim([0 .2])
set(gcf,'Position',[10 10 1600 1600])
if (do_save)
    saveas(gcf,[figFolder 'bout_duration_distributions_by_solution'],'fig')
    saveas(gcf,[figFolder 'bout_duration_distributions_by_solution'],'eps')
    saveas(gcf,[figFolder 'bout_duration_distributions_by_solution'],'png')
end

%% Figure 13
% Exponential fits to bout distributions for each solution
figure;
for i=1:4
    muhat = expfit(allBoutsBySolution{i});
    x = 1:150;
    plot(x,(1/muhat)*exp(-(1/muhat)*x))
    hold on;
end
xlabel('Bout Duration (s)'); legend({'H2O','.01M','.1M','1M'}); ylabel('Probability density')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    saveas(gcf,[figFolder 'bout_duration_exp_fits'],'fig')
    saveas(gcf,[figFolder 'bout_duration_exp_fits'],'eps')
    saveas(gcf,[figFolder 'bout_duration_exp_fits'],'png')
end

%% Figure 14
% Bout durations as a function of solution palatability
figure;
bd = [];
relPals = [];
for i=1:length(animalObjs)
    for j=1:4
        b = [animalObjs(i).boutsBySolution{j}.duration];
        bd = [bd b];
        if (any(isnan(animalObjs(i).relativePalatabilitiesLicks(j))))
            disp([animalObjs(i).name ' ' num2str(j)])
        end
        relPals = [relPals ones(1,length(b))*animalObjs(i).relativePalatabilitiesLicks(j)];
        %scatter(ones(1,length(b))*animalObjs(i).relativePalatabilitiesConsumed(j),b,'k.')
    end
end
[rho,pval] = corr(relPals',bd');
lf = polyfit(relPals,bd,1);
x = min(relPals):.001:max(relPals);
scatter(relPals,bd,'k.'); hold on; plot(x,x*lf(1)+lf(2))
xlabel('Palatability'); ylabel('Bout Duration (s)')
text(.2,550,['\rho = ' num2str(rho)],'fontweight','bold','fontsize',20); 
text(.2, 500,['p = ' num2str(pval)],'fontweight','bold','fontsize',20);
text(.2,450,['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    saveas(gcf,[figFolder 'bout_duration_vs_current_palatability'],'fig')
    saveas(gcf,[figFolder 'bout_duration_vs_current_palatability'],'eps')
    saveas(gcf,[figFolder 'bout_duration_vs_current_palatability'],'png')
end

%% Figure 15
% Bout durations as a function of alternative palatability
figure;
meanAltPals = [];
bd = [];
for i=1:length(animalObjs)
    for j=1:length(animalObjs(i).dates)
        for k=1:length(animalObjs(i).bouts{j})
            %disp([num2str(i) ' ' num2str(j) ' ' num2str(k)])
            if (~isempty(animalObjs(i).bouts{j}{k}))
                meanAltPals = [meanAltPals animalObjs(i).meanAlternativePalatabilities{j}{k}];
                %bd = [bd [animalObjs(i).bouts{j}{k}.duration]];
                bd = [bd [animalObjs(i).bouts{j}{k}.duration]];
            end
        end
    end
end
[rho,pval] = corr(meanAltPals',bd');
lf = polyfit(meanAltPals,bd,1);
x = min(meanAltPals):.001:max(meanAltPals);
scatter(meanAltPals,bd,'k.'); hold on; plot(x,x*lf(1)+lf(2));
xlabel('Alternative palatability'); ylabel('Bout duration (s)')
text(.5,550,['\rho = ' num2str(rho)],'fontsize',20,'fontweight','bold')
text(.5,500,['p = ' num2str(pval)],'fontsize',20,'fontweight','bold')
text(.5,450,['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],'fontsize',20,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    saveas(gcf,[figFolder 'bout_duration_vs_alternative_palatability'],'fig')
    saveas(gcf,[figFolder 'bout_duration_vs_alternative_palatability'],'eps')
    saveas(gcf,[figFolder 'bout_duration_vs_alternative_palatability'],'png')
end

%% Figure 16
% Bout durations as a function of relative palatability difference
figure;
allPalDifs = [];
allBoutDurs = [];
for i=1:length(animalObjs)
    for j=1:length(animalObjs(i).dates)
        for k=1:length(animalObjs(i).bouts{j})
            if (~isempty(animalObjs(i).bouts{j}{k}))
                palDifs = [animalObjs(i).Palatabilities{j}{k}] - [animalObjs(i).meanAlternativePalatabilities{j}{k}];
                boutDurs = [animalObjs(i).bouts{j}{k}.duration];
                allPalDifs = [allPalDifs palDifs];
                allBoutDurs = [allBoutDurs boutDurs];
            end
        end
    end
end
uPals = unique(allPalDifs);
for i=1:length(uPals)
    inds = find(allPalDifs == uPals(i));
    durs = allBoutDurs(inds);
    muhats(i) = expfit(durs);
    medianBoutDurs(i) = median(durs);
    palDifs(i) = allPalDifs(inds(1));
    scatter(allPalDifs(inds),allBoutDurs(inds),'k.'); hold on;
    plot([palDifs(i)-.1 palDifs(i)+.1],[muhats(i) muhats(i)],'r')
end
xlabel('Palatability difference'); ylabel('Bout duration (s)')
set(gcf,'Position',[10 10 1400 1000])
if (do_save)
    saveas(gcf,[figFolder 'bout_duration_vs_palatability_difference'],'fig')
    saveas(gcf,[figFolder 'bout_duration_vs_palatability_difference'],'eps')
    saveas(gcf,[figFolder 'bout_duration_vs_palatability_difference'],'png')
end

%% Figure 17
figure;
[rho,pval] = corr(allPalDifs',allBoutDurs');
scatter(allPalDifs,allBoutDurs,'k.')
lf = polyfit(allPalDifs,allBoutDurs,1);
x = min(allPalDifs):.01:max(allPalDifs);
hold on; plot(x,lf(1)*x + lf(2))
xlabel('Palatability difference'); ylabel('Bout duration (s)')
text(-2.2,460,['\rho = ' num2str(rho)],'fontweight','bold','fontsize',20)
text(-2.2,430,['p = ' num2str(pval)],'fontweight','bold','fontsize',20)
text(-2.2,400,['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],'fontweight','bold','fontsize',20)
ylim([0 max(allBoutDurs)+10])
set(gcf,'Position',[10 10 1400 1000])
if (do_save)
    saveas(gcf,[figFolder 'bout_duration_vs_palatability_difference_linfit_corr'],'fig')
    saveas(gcf,[figFolder 'bout_duration_vs_palatability_difference_linfit_corr'],'eps')
    saveas(gcf,[figFolder 'bout_duration_vs_palatability_difference_linfit_corr'],'png')
end

%% Figure 18
figure;
subplot(1,2,1)
linfit = polyfit(palDifs,muhats,1);
linfit2 = polyfit(allPalDifs,allBoutDurs,1);
[rho,pval] = corr(allPalDifs',allBoutDurs');
x = min(palDifs):.01:max(palDifs);
scatter(palDifs,muhats,'k.'); hold on; h1=plot(x,x*linfit(1) + linfit(2)); 
%h2=plot(x,x*linfit2(1) + linfit2(2));
ylim([0 70])
xlabel('Palatability difference'); ylabel('Exponential mean'); legend([h1],{'Fit to exponential means'})
text(-2.2,60,['\rho = ' num2str(rho)],'fontweight','bold','fontsize',20); 
text(-2.2,55,['p = ' num2str(pval)],'fontweight','bold','fontsize',20)
text(-2.2,50,['y = ' num2str(linfit(1)) 'x + ' num2str(linfit(2))],'fontweight','bold','fontsize',20)

subplot(1,2,2)
linfit3 = polyfit(palDifs,medianBoutDurs,1);
[rho2,pval2] = corr(palDifs',medianBoutDurs');
scatter(palDifs,medianBoutDurs,'k.'); hold on; h3 = plot(x,x*linfit3(1) + linfit3(2));
ylim([0 70])
xlabel('Palatability difference'); ylabel('Median bout length'); legend(h3,{'Fit to median bout lengths'})
text(-2.2,60,['\rho = ' num2str(rho2)],'fontweight','bold','fontsize',20); 
text(-2.2,55,['p = ' num2str(pval2)],'fontweight','bold','fontsize',20)
text(-2.2,50,['y = ' num2str(linfit3(1)) 'x + ' num2str(linfit3(2))],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    saveas(gcf,[figFolder 'bout_duration_vs_palatability_difference_exp_fit_lin_fit'],'fig')
    saveas(gcf,[figFolder 'bout_duration_vs_palatability_difference_exp_fit_lin_fit'],'eps')
    saveas(gcf,[figFolder 'bout_duration_vs_palatability_difference_exp_fit_lin_fit'],'png')
end

%% Figure 19
% bout transition matrix
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
set(gca,'xtick',1:10,'xticklabels',xlabels,'xticklabelrotation',20); xlim([0 11]); ylim([0 1])
ylabel('Transition probability')
hold off;
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    saveas(gcf,[figFolder 'transition_probability_by_solution_pair'],'fig')
    saveas(gcf,[figFolder 'transition_probability_by_solution_pair'],'eps')
    saveas(gcf,[figFolder 'transition_probability_by_solution_pair'],'png')
end

%% Figure 20
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
xlabel('Current - alternative palatability'); ylabel('Transition Probability')
text(0,.95,['\rho = ' num2str(rho1)],'fontweight','bold','fontsize',20)
text(0,.88,['p = ' num2str(pval1)],'fontweight','bold','fontsize',20)
text(0,.81,['y = ' num2str(lf1(1)) 'x + ' num2str(lf1(2))],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    saveas(gcf,[figFolder 'transition_probability_vs_palatability_difference_linfit'],'fig')
    saveas(gcf,[figFolder 'transition_probability_vs_palatability_difference_linfit'],'eps')
    saveas(gcf,[figFolder 'transition_probability_vs_palatability_difference_linfit'],'png')
end

%% Figure 21
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
xlabel('Alternative palatability'); ylabel('Transition Probability')
text(0.3,.95,['\rho = ' num2str(rho3)],'fontweight','bold','fontsize',20)
text(0.3,.88,['p = ' num2str(pval3)],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    saveas(gcf,[figFolder 'transition_probability_vs_alternative_palatability'],'fig')
    saveas(gcf,[figFolder 'transition_probability_vs_alternative_palatability'],'eps')
    saveas(gcf,[figFolder 'transition_probability_vs_alternative_palatability'],'png')
end
%text(0.3,.81,['y = ' num2str(lf3(1)) 'x + ' num2str(lf3(2))],'fontweight','bold','fontsize',20)
%{
figure;
X = [];
GRP = [];
for i=1:length(transitionProbVals)
    X = [X; transitionProbVals{i}'];
    GRP = [GRP; repmat(i,length(transitionProbVals{i}),1)];
end
boxplot(X,GRP)
set(gca,'xtick',1:10,'xticklabels',xlabels,'xticklabelrotation',20); xlim([0 11])
ylabel('Transition probability')
%}

%% Figure 22
% # of bouts for each solution pair
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
set(gca,'xtick',1:7,'xticklabels',xlabels,'xticklabelrotation',20); xlim([0 8])
ylabel('# of bouts')
hold off;
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    saveas(gcf,[figFolder 'number_of_bouts_by_solution'],'fig')
    saveas(gcf,[figFolder 'number_of_bouts_by_solution'],'eps')
    saveas(gcf,[figFolder 'number_of_bouts_by_solution'],'png')
end

%{
figure;
X = [];
GRP = [];
for i=1:length(nBoutsByPair)
    X = [X; nBoutsByPair{i}'];
    GRP = [GRP; repmat(i,length(nBoutsByPair{i}),1)];
end
boxplot(X,GRP)
set(gca,'xtick',1:7,'xticklabels',xlabels,'xticklabelrotation',20); xlim([0 8])
ylabel('# of bouts')
%}

%% Figure 23
% # of switches for each solution pair
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
ylabel('# of switches')
hold off;
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    saveas(gcf,[figFolder 'number_of_switches_by_solution_pair'],'fig')
    saveas(gcf,[figFolder 'number_of_switches_by_solution_pair'],'eps')
    saveas(gcf,[figFolder 'number_of_switches_by_solution_pair'],'png')
end

%{
figure;
X = [];
GRP = [];
for i=1:length(nSwitchVals)
    X = [X; nSwitchVals{i}'];
    GRP = [GRP; repmat(i,length(nSwitchVals{i}),1)];
end
boxplot(X,GRP)
set(gca,'xtick',1:7,'xticklabels',xlabels,'xticklabelrotation',20); xlim([0 8])
ylabel('# of switches')
%}

%% Figure 24
% probability of choosing preferred side from day before first by animal
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
xlabel('Probability of choosing prior preferred side')
ylabel('# of animals')
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'probability_of_choosing_prior_day_preferred_side';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    saveas(gcf,[figFolder figName],'png')
end

figure;
plot(binomPDF); hold on; plot([sum(nPriorPreferred) sum(nPriorPreferred)],[0 .1],'r')
xlabel('# of prior preferred chosen'); ylabel('Probability density')
if (do_save)
    figName = 'prior_day_preferred_binomial_test';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    saveas(gcf,[figFolder figName],'png')
end 

%% Figure 25
% Time to first lick
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
set(gca,'xtick',1:length(setdiff(1:length(animalObjs),5)),'xticklabels',xlabels,'xticklabelrotation',20); xlim([0 length(xlabels)+1])
ylabel('Time to first lick (s)')
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    saveas(gcf,[figFolder 'time_to_first_lick_by_animal'],'fig')
    saveas(gcf,[figFolder 'time_to_first_lick_by_animal'],'eps')
    saveas(gcf,[figFolder 'time_to_first_lick_by_animal'],'png')
end

%% Figure 26
% Distribution of all times to first lick
figure;
histogram(allFirstLickTimes,'binwidth',2)
xlabel('Time to first lick (s)')
ylabel('# of occurrences')
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    saveas(gcf,[figFolder 'time_to_first_lick_distribution'],'fig')
    saveas(gcf,[figFolder 'time_to_first_lick_distribution'],'eps')
    saveas(gcf,[figFolder 'time_to_first_lick_distribution'],'png')
end

%{
% Figure 21
% Two-way anova of current and alternative palatability on bout duration
count = 0;
allMeanDurs = [];
allCurPals = [];
allAltPals = [];
for i=[1 2 3 4 6 7 8]
    for j=1:length(animalObjs(i).dates)
        if (isempty(animalObjs(i).bouts{j}{1}))
            meanDur2 = mean([animalObjs(i).bouts{j}{2}.duration]);
            curpal2 = animalObjs(i).Palatabilities{j}{2}(1);
            altpal2 = animalObjs(i).meanAlternativePalatabilities{j}{2}(1);
            allMeanDurs = [allMeanDurs meanDur2];
            allCurPals = [allCurPals curpal2];
            allAltPals = [allAltPals altpal2];
        elseif (isempty(animalObjs(i).bouts{j}{2}))
            meanDur1 = mean([animalObjs(i).bouts{j}{1}.duration]);
            curpal1 = animalObjs(i).Palatabilities{j}{1}(1);
            altpal1 = animalObjs(i).meanAlternativePalatabilities{j}{1}(1);
            allMeanDurs = [allMeanDurs meanDur1];
            allCurPals = [allCurPals curpal1];
            allAltPals = [allAltPals altpal1];
        else
            meanDur1 = mean([animalObjs(i).bouts{j}{1}.duration]);
            curpal1 = animalObjs(i).Palatabilities{j}{1}(1);
            altpal1 = animalObjs(i).meanAlternativePalatabilities{j}{1}(1);
            meanDur2 = mean([animalObjs(i).bouts{j}{2}.duration]);
            curpal2 = animalObjs(i).Palatabilities{j}{2}(1);
            altpal2 = animalObjs(i).meanAlternativePalatabilities{j}{2}(1);
            allMeanDurs = [allMeanDurs meanDur1 meanDur2];
            allCurPals = [allCurPals curpal1 curpal2];
            allAltPals = [allAltPals altpal1 altpal2];
        end
    end
end
uniqueCurPals = sort(unique(allCurPals));
uniqueAltPals = sort(unique(allAltPals));
for i=1:length(allCurPals)
    allCurPalIDs(i) = find(uniqueCurPals == allCurPals(i));
    allAltPalIDs(i) = find(uniqueAltPals == allAltPals(i));
end
[P,T,STATS,TERMS]=anovan(allMeanDurs,{allCurPalIDs allAltPalIDs},'model','interaction','varnames',{'current palatability','alternative palatability'})
X = zeros(length(uniqueCurPals),length(uniqueAltPals));
outercount = 0;
for i=1:length(uniqueCurPals)
    for j=1:length(uniqueAltPals)
        curdurs = intersect(find(allCurPals == uniqueCurPals(i)),find(allAltPals == uniqueAltPals(j)));
        for k = 1:length(curdurs)
            X(outercount+k,j) = allMeanDurs(curdurs(k));
        end
    end
    outercount = outercount + length(curdurs);
end
%}

%% Figure 27
% Bout duration following switch vs alternative palatability
count = count+1;
durs = [];
altPals = [];
for i=setdiff(1:length(animalObjs),5) % exclude bb12
    for j=1:length(animalObjs(i).dates)
        for k=2:length(animalObjs(i).linBoutsByDay{j})
            if (animalObjs(i).linBoutsByDay{j}(k).channel ~= animalObjs(i).linBoutsByDay{j}(k-1).channel)
                durs = [durs animalObjs(i).linBoutsByDay{j}(k).duration];
                altPals = [altPals animalObjs(i).relativePalatabilitiesLicks(animalObjs(i).solnConverter(animalObjs(i).linBoutsByDay{j}(k-1).solution))];
            end
        end
    end
end
[rho,pval] = corr(durs',altPals');
lf = polyfit(altPals,durs,1);
figure;
scatter(altPals,durs,100,'k.')
x=0:.01:3; hold on; plot(x,lf(1)*x + lf(2))
xlabel('Alternative palatability')
ylabel('Bout duration')
text(1.2,380,['\rho = ' num2str(rho)],'fontweight','bold','fontsize',20)
text(1.2,350,['p = ' num2str(pval)],'fontweight','bold','fontsize',20)
text(1.2,320,['y = ' num2str(lf(1)) 'x + ' num2str(lf(2))],'fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    saveas(gcf,[figFolder 'bout_duration_after_switch_vs_alternative_palatability'],'fig')
    saveas(gcf,[figFolder 'bout_duration_after_switch_vs_alternative_palatability'],'eps')
    saveas(gcf,[figFolder 'bout_duration_after_switch_vs_alternative_palatability'],'png')
end

%% Figure 28
% # of switches by solution pair 1st vs. 2nd session
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
xlabels = {'H2O/H2O','A/A','A/B','A/C','B/B','B/C','C/C'};
figure;
shadedErrorBar(1:7,meanFirstVisitSwitches,stdFirstVisitSwitches,'lineprops','-r')
hold on;
shadedErrorBar(1:7,meanSecondVisitSwitches,stdSecondVisitSwitches,'lineprops','-b')
set(gca,'xtick',1:7,'xticklabels',xlabels)
ylabel('# of switches')
if (do_save)
    figName = 'nSwitches_by_solutionPair_1st_vs_2nd_exposure';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    saveas(gcf,[figFolder figName],'png')
end

%% Figure 29
Y = cell(1,length(animalObjs));
X = cell(1,length(animalObjs));
for i=1:length(animalObjs)
    count = 0;
    for j=1:length(animalObjs(i).dates)
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
    multiLinearBoutDurFits(i).bint = bint;
    multiLinearBoutDurFits(i).r = r;
    multiLinearBoutDurFits(i).rint = rint;
    multiLinearBoutDurFits(i).stats = stats;
    figure;
    scatter3(X{i}(:,1),X{i}(:,2),Y{i},50,'.')
    hold on;
end
z = [multiLinearBoutDurFits.b];
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
figure; ...
for i=1:length(animalObjs)
    Zs{i} = multiLinearBoutDurFits(i).b(2)*(0:.01:3)' + multiLinearBoutDurFits(i).b(1)*(0:.01:3) + multiLinearBoutDurFits(i).b(3);
    surf(0:.01:3,0:.01:3,Zs{i},'FaceColor',cmap(i*floor(size(cmap,1)/length(animalObjs)),:),'FaceAlpha',.3,'EdgeColor',cmap(i*floor(size(cmap,1)/length(animalObjs)),:));
    hold on;
end
xlabel('Current Palatability'); ylabel('Alternative palatability'); zlabel('Predicted Bout duration')
if (do_save)
    figName = 'palatability_alternativePalatability_boutDuration_planes';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    saveas(gcf,[figFolder figName],'png')
end
%% Figure 30
doubleExponential = @(x,lambda1,lambda2) expdf(x,lambda1) + exppdf(x,lambda2);
doubleNormal = @(x,mu1,mu2,std1,std2) normpdf(x,mu1,std1) + normpdf(x,mu2,std2);
distributionNames = {'exponential','poisson','normal','lognormal','gamma','doubleExponential','doubleNormal'};
distributionFuncs = {@exppdf,@poisspdf,@normpdf,@lognpdf,@gampdf,@doubleExponential,@doubleNormal};
bestFitParams = cell(length(animalObjs),4,length(distributionNames));
BIC_fits = zeros(length(animalObjs),4,length(distributionNames));
AIC_fits = zeros(length(animalObjs),4,length(distributionNames));
neglogliks = nan(length(animalObjs),4,length(distributionNames));
bestBICInds = nan(length(animalObjs),4);
bestAICInds = nan(length(animalObjs),4);
x = 0:1:600;
figure;
for i=1:length(animalObjs)
    for j=1:4
        subplot(2,2,j)
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
                startVec = [lambdaStart 10*lambdaStart];
                startVec = sort(startVec,'descend');
                try
                    [phat,pic] = mle(curDurations,'pdf',doubleExponential,'start',startVec);
                    negloglik = -sum(log(doubleExponential(curDurations,phat(1),phat(2))));
                    bestFitParams{i,j,k} = phat;
                    neglogliks(i,j,k) = negloglik;
                    bic = 2*log(curDurationsLength) + 2*negloglik;
                    aic = 4 + 2*negloglik;
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
                stdStartVec = [stdStart stdStart] + randn(1,2)*(stdStart/50);
                startVec = [muStartVec 5*stdStartVec];
                try
                    [phat,pic] = mle(curDurations,'pdf',doubleNormal,'start',startVec);
                    negloglik = -sum(log(doubleNormal(curDurations,phat(1),phat(2),phat(3),phat(4))));
                    bestFitParams{i,j,k} = phat;
                    neglogliks(i,j,k) = negloglik;
                    bic = 4*log(curDurationsLength) + 2*negloglik;
                    aic = 8 + 2*negloglik;
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
doubleNormal = @(x,mu1,mu2,std1,std2) (1/2)*(normpdf(x,mu1,std1) + normpdf(x,mu2,std2));
doubleNormalPhat = cell(1,length(animalObjs));
expPhat = cell(4,length(animalObjs));
iliNegLogLiks = nan(length(animalObjs),4,2);
for i=1:length(animalObjs)
    figure;
    for j=1:4
        if (any(animalObjs(i).ilisBySolution{j} < 0))
            disp([num2str(i) ' ' num2str(j) ' has negative ilis'])
            badinds = find(animalObjs(i).ilisBySolution{j} < 0);
            goodILIs = animalObjs(i).ilisBySolution{j}(setdiff(1:length(animalObjs(i).ilisBySolution),badinds));
        else
            goodILIs = animalObjs(i).ilisBySolution{j};
        end
        curMean = mean(goodILIs);
        curStd = std(goodILIs);
        pd = fitdist(goodILIs','exponential');
        expPhat{j,i} = pd.ParameterValues;
        iliNegLogLiks(i,j,1) = pd.negloglik;
        isDone = 0;
        tryNumber = 0;
        while (~isDone)
            try
                [phat,pic] = mle(goodILIs,'pdf',doubleNormal,'start',[.5*rand 5*rand .5*rand 5*rand]);
                doubleNormalPhat{j,i} = phat;
                isDone = 1;
            catch
                tryNumber = tryNumber + 1;
                disp(['doubleNorm try# ' num2str(tryNumber)])
            end
        end
        negloglik = -sum(log(doubleNormal(goodILIs,phat(1),phat(2),phat(3),phat(4))));
        iliNegLogLiks(i,j,2) = negloglik;
        subplot(2,2,j)
        histogram(goodILIs,'normalization','pdf','binwidth',.005);
        hold on;
        xvec = 0:.001:max(goodILIs);
        plot(xvec,doubleNormal(xvec,phat(1),phat(2),phat(3),phat(4)))
        xlabel('ILI (s)')
        ylabel('Probability Density')
    end
end
