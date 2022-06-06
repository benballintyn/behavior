% pizza talk figs
clear all;
figFolder = '/home/ben/phd/behavior/figures/10_11_20/';
if (~exist(figFolder,'dir'))
    mkdir(figFolder)
end
do_save = 1; % Boolean for whether to save figures or not
animals = {'bb8','bb9','bb10','bb11','bb12','bb13','bb14','bb15','bb16','bb17','bb18','bb19'};
dates1 = {'190629','190630','190701','190702','190703','190704','190705','190808','190809','190810','190811','190812','190813','190814'};
dates2 = {'190816','190817','190818','190819','190820','190821','190822','190823','190824','190825','190826','190827','190828'};
dates3 = {'191027','191028','191029','191030','191031','191101','191102','191103','191104','191105','191106','191107','191108','191109'};
dates4 = {'200124','200125','200126','200127','200128','200129','200130','200131','200201','200202','200203','200204','200205','200206'};
dates5 = {'200124','200125','200126','200127','200128','200129','200130','200131','200202','200203','200204','200205','200206','200207'}; % bb18 does not have 200201
dates6 = {'200209','200210','200211','200212','200213','200214','200215','200216','200217','200218','200219','200220','200221','200222'};
dates{1} = dates1; dates{2}=dates1; 
dates{3} = dates2; dates{4} = dates2; dates{5} = dates2; dates{6} = dates2; 
dates{7} = dates3; dates{8} = dates3;
dates{9} = dates4; dates{10} = dates4;
dates{11} = dates5;
dates{12} = dates4;
%dates{13} = dates6;
omitDays = [5 12];
if (~exist('animalObjs','var'))
    for i=1:length(animals)
        animalObjs(i) = animalData(animals{i},dates{i},'analyzed_data');
        disp(['done with ' animals{i}])
    end
end
B = load(['~/phd/talks/Pizza_Talks/nov_2019/B2.mat']); B = B.B;
f=@(B,x) B(1).*exp(B(2).*x) + B(3);

% get all bouts
allBouts = [];
for i=1:length(animals)
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


% FIGURE 1
% cdf plot of lick times
figure;
for i=1:length(animals)
    cdfplot([animalObjs(i).linLicks.onset])
    hold on;
end
xlabel('Time (s)'); ylabel('P(t_{lick} < x)')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    saveas(gcf,[figFolder 'lick_times_cdf'],'fig')
    saveas(gcf,[figFolder 'lick_times_cdf'],'eps')
end

% Figure 2
% cdf plot of bout onset times
figure;
for i=1:length(animals)
    cdfplot([animalObjs(i).linBouts.onset])
    hold on;
end
xlabel('Time (s)'); ylabel('P(t_{bout} < x)')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    saveas(gcf,[figFolder 'bout_onset_times_cdf'],'fig')
    saveas(gcf,[figFolder 'bout_onset_times_cdf'],'eps')
end

% Figure 3
% bout lengths vs time
figure;
for i=1:length(animals)
    scatter([animalObjs(i).linBouts.onset],[animalObjs(i).linBouts.duration],100,'k.')
    hold on;
end
xlabel('Time(s)'); ylabel('Bout duration (s)')
x = 300:1:3900;
plot(x,f(B,x-300))
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    saveas(gcf,[figFolder 'bout_lengths_vs_time'],'fig')
    saveas(gcf,[figFolder 'bout_lengths_vs_time'],'eps')
end

% Figure 4
% bout duration as a function of previous bout duration
figure;
xs = [];
ys = [];
for i=1:length(animals)
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
end

% Figure 5
% total consumed for all animals
xlabels = {'H2O','.01M','.1M','1M'};
figure;
for i=1:length(animals)
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
end

% Figure 6
% total licks for all animals
figure;
for i=1:length(animals)
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
end

% Figure 7
% relative palatability by amount consumed or total licks
figure;
for i=1:length(animals)
    rpalLicks(:,i) = animalObjs(i).relativePalatabilitiesLicks;
    rpalConsumed(:,i) = animalObjs(i).relativePalatabilitiesConsumed;
end
subplot(1,2,1); plot(1:4,rpalConsumed); set(gca,'xticklabels',xlabels); ylabel('Relative Palatabilitiy'); title('Relative Palatability based on consumption')
subplot(1,2,2); plot(1:4,rpalLicks); set(gca,'xticklabels',xlabels); ylabel('RelativePalatability'); title('Relative Palatability based on lick totals')
set(gcf,'Position',[10 10 1800 1000])
if (do_save)
    saveas(gcf,[figFolder 'relative_palatability_by_amount_consumed_or_licks'],'fig')
    saveas(gcf,[figFolder 'relative_palatability_by_amount_consumed_or_licks'],'eps')
end

% Figure 8
% Bout duration distributions for each solution
figure;
allBoutsBySolution = cell(1,4);
maxDur = 0;
for i=1:length(animals)
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
end

% Figure 9
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
end

% Figure 10
% Bout durations as a function of solution palatability
figure;
bd = [];
relPals = [];
for i=1:length(animals)
    for j=1:4
        b = [animalObjs(i).boutsBySolution{j}.duration];
        correctedB = b - f(B,[animalObjs(i).boutsBySolution{j}.onset]);
        bd = [bd b];
        relPals = [relPals ones(1,length(b))*animalObjs(i).relativePalatabilitiesConsumed(j)];
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
end

% Figure 11
% Bout durations as a function of alternative palatability
figure;
meanAltPals = [];
bd = [];
for i=1:length(animals)
    for j=1:length(animalObjs(i).dates)
        for k=1:length(animalObjs(i).bouts{j})
            %disp([num2str(i) ' ' num2str(j) ' ' num2str(k)])
            if (~isempty(animalObjs(i).bouts{j}{k}))
                meanAltPals = [meanAltPals animalObjs(i).meanAlternativePalatabilities{j}{k}];
                correctedB = [animalObjs(i).bouts{j}{k}.duration] - f(B,[animalObjs(i).bouts{j}{k}.onset]);
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
end

% Figure 12
% Bout durations as a function of relative palatability difference
figure;
allPalDifs = [];
allBoutDurs = [];
correctedAllBoutDurs = [];
for i=1:length(animals)
    for j=1:length(animalObjs(i).dates)
        for k=1:length(animalObjs(i).bouts{j})
            if (~isempty(animalObjs(i).bouts{j}{k}))
                palDifs = [animalObjs(i).Palatabilities{j}{k}] - [animalObjs(i).meanAlternativePalatabilities{j}{k}];
                boutDurs = [animalObjs(i).bouts{j}{k}.duration];
                correctedBoutDurs = boutDurs - f(B,[animalObjs(i).bouts{j}{k}.onset]);
                allPalDifs = [allPalDifs palDifs];
                allBoutDurs = [allBoutDurs boutDurs];
                correctedAllBoutDurs = [correctedAllBoutDurs correctedBoutDurs];
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
end

% Figure 13
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
end

% Figure 14
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
end

% Figure 15
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
for k=setdiff(1:length(animals),5) % exclude bb12
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
end

slns = [1 1; 2 2; 3 2; 4 2; 2 3; 3 3; 4 3; 2 4; 3 4; 4 4];
relPalDifs = cell(1,length(transitionProbVals));
altPals = cell(1,size(slns,1));
curPals = cell(1,size(slns,1));
for i=1:size(slns,1)
    count = 0;
    for k=setdiff(1:length(animals),5) % exclude bb12
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
end

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

% Figure 16
% # of bouts for each solution pair
figure;
xlabels = {'H2O/H2O','A/A','A/B','A/C','B/B','B/C','C/C'};
nBoutsByPair = cell(1,7);
for i=setdiff(1:length(animals),5) % exclude bb12
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

% Figure 17
% # of switches for each solution pair
figure;
xlabels = {'H2O/H2O','A/A','A/B','A/C','B/B','B/C','C/C'};
nSwitchVals = cell(1,7);
for i=setdiff(1:length(animals),5) % exclude bb12
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

% Figure 18
% probability of choosing preferred side from day before first by animal
nPriorPreferred = zeros(1,length(setdiff(1:length(animals),5)));
nPriorPreferredFrac = zeros(1,length(setdiff(1:length(animals),5)));
count = 0;
for i=setdiff(1:length(animals),5) % exclude bb12
    count = count + 1;
    for j=2:length(animalObjs(i).dates)
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
figure;
histogram(nPriorPreferredFrac,'binwidth',.05)
xlabel('Probability of choosing prior preferred side')
ylabel('# of animals')
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    saveas(gcf,[figFolder 'probability_of_choosing_prior_day_preferred_side'],'fig')
    saveas(gcf,[figFolder 'probability_of_choosing_prior_day_preferred_side'],'eps')
end

% Figure 19
% Time to first lick
xlabels = {'bb8','bb9','bb10','bb11','bb13','bb14','bb15','bb16','bb17','bb18','bb19'};
count = 0;
time2firstLick = cell(1,length(setdiff(1:length(animals),5)));
allFirstLickTimes = [];
for i=setdiff(1:length(animals),5) % exclude bb12
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
set(gca,'xtick',1:length(setdiff(1:length(animals),5)),'xticklabels',xlabels,'xticklabelrotation',20); xlim([0 length(xlabels)+1])
ylabel('Time to first lick (s)')
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    saveas(gcf,[figFolder 'time_to_first_lick_by_animal'],'fig')
    saveas(gcf,[figFolder 'time_to_first_lick_by_animal'],'eps')
end

% Figure 20
% Distribution of all times to first lick
figure;
histogram(allFirstLickTimes,'binwidth',2)
xlabel('Time to first lick (s)')
ylabel('# of occurrences')
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    saveas(gcf,[figFolder 'time_to_first_lick_distribution'],'fig')
    saveas(gcf,[figFolder 'time_to_first_lick_distribution'],'eps')
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

% Figure 21
% Bout duration following switch vs alternative palatability
count = count+1;
durs = [];
altPals = [];
for i=setdiff(1:length(animals),5) % exclude bb12
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
end