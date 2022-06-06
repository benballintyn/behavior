%% make_figures2
clear all; close all;
figFolder = '/home/ben/phd/behavior/figures/alteredPalatability/08_27_21/';
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
% Alter palatability of .01M to be 30% higher than .1M
for i=1:length(animalObjs)
    animalObjs(i).relativePalatabilitiesLicks(2) = 1.3*animalObjs(i).relativePalatabilitiesLicks(3);
end
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
    figName = 'bout_duration_vs_current_palatability_alteredPalatability';
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
    figName = 'bout_duration_vs_current_palatability_difSolDays_alteredPalatability';
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
    figName = 'bout_duration_vs_alternative_palatability_alteredPalatability';
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
    muhats(i) = expfit(durs);
    medianBoutDurs(i) = median(durs);
    palDifs(i) = allPalDifs(inds(1));
    scatter(allPalDifs(inds),allBoutDurs(inds),'k.'); hold on;
    plot([palDifs(i)-.1 palDifs(i)+.1],[muhats(i) muhats(i)],'r')
end
xlabel('Palatability difference','fontsize',15,'fontweight','bold');
ylabel('Bout duration (s)','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1400 1000])
if (do_save)
    figName = 'bout_duration_vs_palatability_difference_alteredPalatability';
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
    figName = 'bout_duration_vs_palatability_difference_linfit_corr_alteredPalatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 18: Bout duration Exponential means and medians as a function of palatability difference
figure;
subplot(1,2,1)
linfit = polyfit(palDifs,muhats,1);
linfit2 = polyfit(allPalDifs,allBoutDurs,1);
[rho,pval] = corr(palDifs',muhats');
x = min(palDifs):.01:max(palDifs);
scatter(palDifs,muhats,'k.'); hold on; h1=plot(x,x*linfit(1) + linfit(2)); 
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
    figName = 'bout_duration_vs_palatability_difference_exp_fit_lin_fit_alteredPalatability';
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
if (do_save)
    figName = 'palatability_alternativePalatability_boutDuration_planes_alteredPalatability';
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
    figName = 'all_bouts_palatability_coefficients_alteredPalatability';
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
    figName = 'all_bouts_normalized_palatability_coefficients_alteredPalatability';
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
    figName = 'bout_duration_after_stay_vs_alternative_palatability_alteredPalatability';
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
    figName = 'bout_duration_after_stay_vs_palatability_alteredPalatability';
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
    figName = 'bout_duration_after_switch_vs_alternative_palatability_alteredPalatability';
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
    figName = 'bout_duration_after_switch_vs_palatability_alteredPalatability';
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
    figName = 'Stay_vs_switch_palatability_coefficients_alteredPalatability';
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
    figName = 'Stay_vs_switch_alternative_palatability_coefficients_alteredPalatability';
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
    figName = 'Stay_vs_switch_normalized_palatability_coefficients_alteredPalatability';
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
    figName = 'Stay_vs_switch_normalized_alternative_palatability_coefficients_alteredPalatability';
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
    figName = 'Early_vs_late_palatability_coefficients_alteredPalatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

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
    figName = 'Early_vs_late_alternative_palatability_coefficients_alteredPalatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
[p1,h] = signrank(abs(earlyBs(:,1)),abs(lateBs(:,1)),'tail','right');
[p2,h] = signrank(abs(lateBs(:,1)),abs(earlyBs(:,1)),'tail','right');
if (p1 < .05)
    disp(['Early palatability coefficient magnitudes are significantly higher than late coefficient magnitudes. p = ' num2str(p1)])
elseif (p2 < .05)
    disp(['Late palatability coefficient magnitudes are significantly higher than early coefficient magnitudes. p = ' num2str(p2)])
else
    disp(['Neither early or late palatability coefficient magnitudes are significantly higher than each other'])
end
notBoxPlot([abs(earlyBs(:,1)); abs(lateBs(:,1))],[ones(size(earlyBs,1),1); ones(size(lateBs,1),1)*2])
ylabel('Palatablity coefficient magnitude','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Early','Late'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'early_vs_late_palatability_coefficient_magnitudes_alteredPalatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

figure;
[p1,h] = signrank(abs(earlyBs(:,2)),abs(lateBs(:,2)),'tail','right');
[p2,h] = signrank(abs(lateBs(:,2)),abs(earlyBs(:,2)),'tail','right');
if (p1 < .05)
    disp(['Early alternative palatability coefficient magnitudes are significantly higher than late coefficient magnitudes. p = ' num2str(p1)])
elseif (p2 < .05)
    disp(['Late alternative palatability coefficient magnitudes are significantly higher than early coefficient magnitudes. p = ' num2str(p2)])
else
    disp(['Neither early or late alternative palatability coefficient magnitudes are significantly higher than each other'])
end
notBoxPlot([abs(earlyBs(:,2)); abs(lateBs(:,2))],[ones(size(earlyBs,1),1); ones(size(lateBs,1),1)*2])
ylabel('Alternative palatablity coefficient magnitude','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Early','Late'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'early_vs_late_alternative_palatability_coefficient_magnitudes_alteredPalatability';
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
    figName = 'Early_vs_late_normalized_palatability_coefficients_alteredPalatability';
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
    figName = 'Early_vs_late_normalized alternative_palatability_coefficients_alteredPalatability';
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
xlabels = {'H2O->H2O','.01M->.01M','.1M->.01M','1M->.01M','.01M->.1M','.1M->.1M','1M->.1M','.01M->1M','.1M->1M','1M->1M'};
set(gca,'xtick',1:10,'xticklabels',xlabels,'xticklabelrotation',20,'fontsize',15,'fontweight','bold'); xlim([0 11]); ylim([0 1])
ylabel('Transition probability','fontsize',15,'fontweight','bold')
hold off;
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'transition_probability_by_solution_pair_alteredPalatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
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
xlabels = {'H2O->H2O','.01M->.01M','.1M->.01M','1M->.01M','.01M->.1M','.1M->.1M','1M->.1M','.01M->1M','.1M->1M','1M->1M'};
set(gca,'xtick',1:10,'xticklabels',xlabels,'xticklabelrotation',20,'fontsize',15,'fontweight','bold'); xlim([0 11]); ylim([0 1])
ylabel('Transition probability','fontsize',15,'fontweight','bold')
title('Early bouts only','fontsize',15,'fontweight','bold')
hold off;
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'early_transition_probability_by_solution_pair_alteredPalatability';
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
xlabels = {'H2O->H2O','.01M->.01M','.1M->.01M','1M->.01M','.01M->.1M','.1M->.1M','1M->.1M','.01M->1M','.1M->1M','1M->1M'};
set(gca,'xtick',1:10,'xticklabels',xlabels,'xticklabelrotation',20,'fontsize',15,'fontweight','bold'); xlim([0 11]); ylim([0 1])
ylabel('Transition probability','fontsize',15,'fontweight','bold')
title('Late bouts only','fontsize',15,'fontweight','bold')
hold off;
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'late_transition_probability_by_solution_pair_alteredPalatability';
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
    figName = 'transition_probability_vs_palatability_difference_linfit_alteredPalatability';
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
    figName = 'transition_probability_vs_current_palatability_alteredPalatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
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
    figName = 'transition_probability_vs_alternative_palatability_alteredPalatability';
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

glm = fitglm([pastPalatability' alternativePalatability'],wasSwitch','interactions','Distribution','binomial','Link','logit');
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
xlabel('Current palatability'); ylabel('Alternative palatability'); zlabel('Switch probability')

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
[stayGLMCurPalCoeffsNoOutliers,stayGLMCurPalCoeffsOutlierInds] = rmoutliers(stayGLMCurPalCoeffs);
[switchGLMCurPalCoeffsNoOutliers,switchGLMCurPalCoeffsOutlierInds] = rmoutliers(switchGLMCurPalCoeffs);
notBoxPlot([stayGLMCurPalCoeffsNoOutliers switchGLMCurPalCoeffsNoOutliers],[ones(1,length(stayGLMCurPalCoeffsNoOutliers)) ones(1,length(switchGLMCurPalCoeffsNoOutliers))*2])
set(gca,'xtick',[1 2],'xticklabels',{'Stay','Switch'},'fontsize',15,'fontweight','bold')
ylabel('Palatability coefficient','fontsize',15,'fontweight','bold')
subplot(1,2,2)
[stayGLMAltPalCoeffsNoOutliers,stayGLMAltPalCoeffsOutlierInds] = rmoutliers(stayGLMAltPalCoeffs);
[switchGLMAltPalCoeffsNoOutliers,switchGLMAltPalCoeffsOutlierInds] = rmoutliers(switchGLMAltPalCoeffs);
notBoxPlot([stayGLMAltPalCoeffsNoOutliers switchGLMAltPalCoeffsNoOutliers],[ones(1,length(stayGLMAltPalCoeffsNoOutliers)) ones(1,length(switchGLMAltPalCoeffsNoOutliers))*2])
set(gca,'xtick',[1 2],'xticklabels',{'Stay','Switch'},'fontsize',15,'fontweight','bold')
ylabel('Alternative palatability coefficient','fontsize',15,'fontweight','bold')
suptitle('Logistic regression coefficients: Outliers removed')
set(gcf,'Position',[10 10 1400 1000])
if (do_save)
    figName = 'logistic_regression_coefficients_palatability_alternative_palatability_alteredPalatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

%% Plot early and late switch probability surfaces
surf(palatabilityVals,alternative_palatabilityVals,earlyPredictedVals','FaceColor','b','EdgeColor','b','FaceAlpha',.3)
hold on;
surf(palatabilityVals,alternative_palatabilityVals,latePredictedVals','FaceColor','r','EdgeColor','r','FaceAlpha',.3)
xlabel('Current palatability'); ylabel('Alternative palatability'); zlabel('Switch probability')
legend({'Early','Late'})

