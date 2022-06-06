%% Load in data for figures using our bout criterion
clear all; close all;
publicationFolder = '/home/ben/phd/behavior/publication/';
dataFolder = [publicationFolder 'data/'];
figFolder = [publicationFolder 'figures/'];
if (~exist(figFolder,'dir'))
    mkdir(figFolder)
end
if (~exist(dataFolder,'dir'))
    mkdir(dataFolder)
end
do_save = 1; % Boolean for whether to save figures or not
animalInfo = load('reduced_data/animalInfo.mat'); animalInfo=animalInfo.animalInfo;

% load data into animal objects
if (~exist('animalObjs','var'))
    for i=1:length(animalInfo)
        animalObjs(i) = animalData(animalInfo(i).animal,animalInfo(i).dates,animalInfo(i).sex,'reduced_data');
        disp(['done with ' animalInfo(i).animal])
    end
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
save([dataFolder 'boutTable.mat'],'boutTable','-mat')

%% Creat indices variables for relevant subsets of bouts
difSolDayInds = find(boutTable.solutionNum ~= boutTable.alternative_solutionNum);
Ainds = find(boutTable.solutionNum == 2);
Binds = find(boutTable.solutionNum == 3);
Cinds = find(boutTable.solutionNum == 4);
isEarly = find(boutTable.early);
isLate = find(boutTable.late);

%% Figure 2A: Relative palatability by solution and animal sex
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
    figName = 'fig2A_relativePalatabilities_by_solution_and_sex';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
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

for i=2:4
    [h,p,ci,stats] = ttest2(maleRelativePalatabilities{i},femaleRelativePalatabilities{i});
    if (p < .05)
        disp(['Male and female palatabilities for solution ' num2str(i) 'are significantly different. t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p)])
    else
        disp(['No sex specific differences for solution ' num2str(i) '. t(' num2str(stats.df) ') = ' num2str(stats.tstat) ', p = ' num2str(p)])
    end
end

%% FIGURE 2B: cdf plot of lick times
cdfs = cell(1,length(animalObjs));
figure;
for i=1:length(animalObjs)
    cdfplot([animalObjs(i).linLicks.onset])
    hold on;
end
xlabel('Time (s)','fontsize',15,'fontweight','bold'); ylabel('P(t_{lick} < x)','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1000])
if (do_save)
    figName = 'fig2B_lick_times_cdf';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

%% Figure 2C: exponential fits to bout durations, log plot
%Adurations = boutTable.duration(intersect(Ainds,difSolDayInds));
%Bdurations = boutTable.duration(intersect(Binds,difSolDayInds));
%Cdurations = boutTable.duration(intersect(Cinds,difSolDayInds));
minDurCutoff = 2;
Adurations = boutTable.duration(Ainds); Adurations = Adurations(Adurations > minDurCutoff);
Bdurations = boutTable.duration(Binds); Bdurations = Bdurations(Bdurations > minDurCutoff);
Cdurations = boutTable.duration(Cinds); Cdurations = Cdurations(Cdurations > minDurCutoff);


BIN_WIDTH_A = 1;
BIN_WIDTH_B = 1;
BIN_WIDTH_C = .75;

[nA,edgesA] = histcounts(Adurations,'binwidth',BIN_WIDTH_A);
[nB,edgesB] = histcounts(Bdurations,'binwidth',BIN_WIDTH_B);
[nC,edgesC] = histcounts(Cdurations,'binwidth',BIN_WIDTH_C);

Azeros = find(nA == 0);
Bzeros = find(nB == 0);
Czeros = find(nC == 0);
firstAzero = Azeros(1);
firstBzero = Bzeros(1);
firstCzero = Czeros(1);

Ax_vals = (minDurCutoff + (BIN_WIDTH_A/2)):BIN_WIDTH_A:(minDurCutoff + (firstAzero - 1)*BIN_WIDTH_A);
Bx_vals = (minDurCutoff + (BIN_WIDTH_B/2)):BIN_WIDTH_B:(minDurCutoff + (firstBzero - 1)*BIN_WIDTH_B);
Cx_vals = (minDurCutoff + (BIN_WIDTH_C/2)):BIN_WIDTH_C:(minDurCutoff + (firstCzero - 1)*BIN_WIDTH_C);

linFitA = polyfit(Ax_vals,log(nA(1:(firstAzero - 1))),1);
linFitB = polyfit(Bx_vals,log(nB(1:(firstBzero - 1))),1);
linFitC = polyfit(Cx_vals,log(nC(1:(firstCzero - 1))),1);

figure;
plot(Ax_vals,nA(1:(firstAzero - 1))); set(gca,'Yscale','log')
hold on;
plot(Ax_vals,exp(linFitA(1)*Ax_vals + linFitA(2)))
ylabel('No. of occurences (log-scale)','fontweight','bold','fontsize',20)
xlabel('Bout duration (sec)','fontweight','bold','fontsize',20)
legend({'Data','Linear model'})
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'fig2C_exponential_fit_solnA';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

figure;
plot(Bx_vals,nB(1:(firstBzero - 1))); set(gca,'Yscale','log')
hold on;
plot(Bx_vals,exp(linFitB(1)*Bx_vals + linFitB(2)))
ylabel('No. of occurences (log-scale)','fontweight','bold','fontsize',20)
xlabel('Bout duration (sec)','fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'fig2C_exponential_fit_solnB';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

figure;
plot(Cx_vals,nC(1:(firstCzero - 1))); set(gca,'Yscale','log')
hold on;
plot(Cx_vals,exp(linFitC(1)*Cx_vals + linFitC(2)))
ylabel('No. of occurences (log-scale)','fontweight','bold','fontsize',20)
xlabel('Bout duration (sec)','fontweight','bold','fontsize',20)
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'fig2C_exponential_fit_solnC';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

%% Run multilinear regressions
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

%% Get R-squared for all bout regressions
allX = [];
allY = [];
for i=1:length(animalObjs)
    ratInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    curInds = intersect(ratInds,difSolDayInds);
    [uniqueRows,ia,ic] = unique([boutTable.palatability(curInds) boutTable.alternative_palatability(curInds)],'rows');
    curX = [];
    curY = [];
    for j = 1:length(ia)
        curTableInds = curInds(find(ic == j));
        curX = [curX; uniqueRows(j,:) 1];
        curY = [curY; mean(boutTable.duration(curTableInds))];
    end
    [b,bint,r,rint,stats] = regress(curY,curX);
    Rsquared(i) = stats(1);
    pvalues(i) = stats(3);
    coeffs(i,:) = b;
    normalizedCoeffs(i,:) = b./mean(curY);
    
    allX = [allX; curX];
    allY = [allY; curY];
end

palatabilityRange = 0:.001:2;
[b,bint,r,rint,stats] = regress(allY,allX);
plane_allRats = palatabilityRange*b(1) + palatabilityRange'*b(2) + b(3);
figure;
scatter3(allX(:,1),allX(:,2),allY,'.')
hold on;
surf(palatabilityRange,palatabilityRange,plane_allRats,'FaceColor','g','EdgeColor','g')
xlabel('Current'); ylabel('Alternative'); zlabel('Duration (sec)')

figure;
notBoxPlot([normalizedCoeffs(:,1) normalizedCoeffs(:,2)])
[p,h,stats] = signrank(normalizedCoeffs(:,1),0,'tail','right');
if (p < .05)
    disp(['(mean duration) palatability coefficients are significantly greater than 0. z = ' num2str(stats.zval) ', p = ' num2str(p)])
end
[p,h,stats] = signrank(normalizedCoeffs(:,2),0,'tail','left');
if (p < .05)
    disp(['(mean duration) alternative palatability coefficients are significantly less than 0. z = ' num2str(stats.zval) ', p = ' num2str(p)])
end

%% Fig 3A: Example of multilinear regression (need to run previous section)
animal2use = 1;
vals = 0:.01:2;
figure;
subplot(1,2,1)
plot(vals,z(1,animal2use)*vals + z(2,animal2use)*animalObjs(animal2use).relativePalatabilitiesLicks(2) + z(3,animal2use))
hold on;
plot(vals,z(1,animal2use)*vals + z(2,animal2use)*animalObjs(animal2use).relativePalatabilitiesLicks(3) + z(3,animal2use))
plot(vals,z(1,animal2use)*vals + z(2,animal2use)*animalObjs(animal2use).relativePalatabilitiesLicks(4) + z(3,animal2use))
legend({'Alternative: .01M','Alternative: .1M','Alternative: 1M'},'fontsize',15,'fontweight','bold')
xlabel('Current palatability','fontsize',15,'fontweight','bold')
ylabel('Bout duration (s)','fontsize',15,'fontweight','bold')
ylim([0 25])

subplot(1,2,2)
plot(vals,z(1,animal2use)*animalObjs(animal2use).relativePalatabilitiesLicks(2) + z(2,animal2use)*vals + z(3,animal2use))
hold on;
plot(vals,z(1,animal2use)*animalObjs(animal2use).relativePalatabilitiesLicks(3) + z(2,animal2use)*vals + z(3,animal2use))
plot(vals,z(1,animal2use)*animalObjs(animal2use).relativePalatabilitiesLicks(4) + z(2,animal2use)*vals + z(3,animal2use))
legend({'Current: .01M','Current: .1M','Current: 1M'},'fontsize',15,'fontweight','bold')
xlabel('Alternative palatability','fontsize',15,'fontweight','bold')
ylabel('Bout duration (s)','fontsize',15,'fontweight','bold')
%ylim([0 25])
set(gcf,'Position',[10 10 1400 600])

if (do_save)
    figName = 'fig3A_Multilinear_regression_example';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

%% Fig 3B: Normalized regression coefficients boxplot
figure;
notBoxPlot([normalizedZ(1,:) normalizedZ(2,:)],[ones(1,length(normalizedZ(1,:))) ones(1,length(normalizedZ(2,:)))*2])
ylabel('Normalized regression coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Palatability (\alpha)','Alternative palatability (\beta)'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'fig3B_all_bouts_normalized_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

%% Run multilinear regressions for early vs. late bouts
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
    
    [uniqueRows,ia,ic] = unique(earlyX,'rows');
    meanEarlyX = [];
    meanEarlyY = [];
    for j=1:length(ia)
        curXinds = find(ic == j);
        meanEarlyX = [meanEarlyX; uniqueRows(j,:) 1];
        meanEarlyY = [meanEarlyY; mean(earlyY(curXinds))];
    end
    [b,bint,r,rint,stats] = regress(meanEarlyY,meanEarlyX);
    earlyRsquared(i) = stats(1);
    earlyPvalues(i) = stats(3);
    earlyMeanCoeffs(i,:) = b;
    earlyMeanNormalizedCoeffs(i,:) = b./mean(meanEarlyY);
    
    [uniqueRows,ia,ic] = unique(lateX,'rows');
    meanLateX = [];
    meanLateY = [];
    for j=1:length(ia)
        curXinds = find(ic == j);
        meanLateX = [meanLateX; uniqueRows(j,:) 1];
        meanLateY = [meanLateY; mean(lateY(curXinds))];
    end
    [b,bint,r,rint,stats] = regress(meanLateY,meanLateX);
    lateRsquared(i) = stats(1);
    latePvalues(i) = stats(3);
    lateMeanCoeffs(i,:) = b;
    lateMeanNormalizedCoeffs(i,:) = b./mean(meanLateY);
end

figure;
notBoxPlot([earlyMeanNormalizedCoeffs(:,1) lateMeanNormalizedCoeffs(:,1)])
[p,h,stats] = signrank(earlyMeanNormalizedCoeffs(:,1),lateMeanNormalizedCoeffs(:,1));
if (p < .05)
    disp(['(mean duration) early and late palatability coefficients are significantly different. z = ' num2str(stats.zval) ', p = ' num2str(p)])
else
    disp(['(mean duration) early and late palatability coefficients are not significantly different. z = ' num2str(stats.zval) ', p = ' num2str(p)])
end

figure;
notBoxPlot([earlyMeanNormalizedCoeffs(:,2) lateMeanNormalizedCoeffs(:,2)])
[p,h,stats] = signrank(earlyMeanNormalizedCoeffs(:,2),lateMeanNormalizedCoeffs(:,2));
if (p < .05)
    disp(['(mean duration) early and late alternative palatability coefficients are significantly different. z = ' num2str(stats.zval) ', p = ' num2str(p)])
else
    disp(['(mean duration) early and late alternative palatability coefficients are not significantly different. z = ' num2str(stats.zval) ', p = ' num2str(p)])
end

%% Fig 4: early vs. late normalized regression coefficients (need to run prior section)
figure;
notBoxPlot([earlyNormalizedBs(:,1); lateNormalizedBs(:,1)],[ones(size(earlyNormalizedBs,1),1); ones(size(lateNormalizedBs,1),1)*2])
ylabel('Normalized palatability coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'early','late'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])

if (do_save)
    figName = 'fig4_early_vs_late_normalized_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

figure;
notBoxPlot([earlyNormalizedBs(:,2); lateNormalizedBs(:,2)],[ones(size(earlyNormalizedBs,1),1); ones(size(lateNormalizedBs,1),1)*2])
ylabel({'Normalized Alternative','palatability coefficient'},'fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'early','late'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])

if (do_save)
    figName = 'fig4_early_vs_late_normalized_alternative_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

%% Run multilinear regressions for stay vs. switch bouts
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
    
    [uniqueRows,ia,ic] = unique(switchX,'rows');
    meanSwitchX = [];
    meanSwitchY = [];
    for j=1:length(ia)
        curXinds = find(ic == j);
        meanSwitchX = [meanSwitchX; uniqueRows(j,:) 1];
        meanSwitchY = [meanSwitchY; mean(switchY(curXinds))];
    end
    [b,bint,r,rint,stats] = regress(meanSwitchY,meanSwitchX);
    switchRsquared(i) = stats(1);
    switchPvalues(i) = stats(3);
    switchMeanCoeffs(i,:) = b;
    switchMeanNormalizedCoeffs(i,:) = b./mean(meanSwitchY);
    
    stayX = [stayX ones(size(stayX,1),1)];
    [b,bint,r,rint,stats] = regress(stayY',stayX);
    stayBs(count,:) = b;
    stayNormalizedBs(count,:) = b./mean(stayY);
    stayBINTs{count} = bint;
    stayRs{count} = r;
    stayRINTs{count} = rint;
    staySTATs{count} = stats;
    
    [uniqueRows,ia,ic] = unique(stayX,'rows');
    meanStayX = [];
    meanStayY = [];
    for j=1:length(ia)
        curXinds = find(ic == j);
        meanStayX = [meanStayX; uniqueRows(j,:) 1];
        meanStayY = [meanStayY; mean(stayY(curXinds))];
    end
    [b,bint,r,rint,stats] = regress(meanStayY,meanStayX);
    stayRsquared(i) = stats(1);
    stayPvalues(i) = stats(3);
    stayMeanCoeffs(i,:) = b;
    stayMeanNormalizedCoeffs(i,:) = b./mean(meanStayY);
    
    firstStayX = [firstStayX ones(size(firstStayX,1),1)];
    [b,bint,r,rint,stats] = regress(firstStayY',firstStayX);
    firstStayBs(count,:) = b;
    firstStayNormalizedBs(count,:) = b./mean(stayY);
    firstStayBINTs{count} = bint;
    firstStayRs{count} = r;
    firstStayRINTs{count} = rint;
    firstStaySTATs{count} = stats;
end

figure;
notBoxPlot([stayMeanNormalizedCoeffs(:,1) switchMeanNormalizedCoeffs(:,1)])
[p,h,stats] = signrank(stayMeanNormalizedCoeffs(:,1),switchMeanNormalizedCoeffs(:,1));
if (p < .05)
    disp(['(mean duration) stay and switch palatability coefficients are significantly different. z = ' num2str(stats.zval) ', p = ' num2str(p)])
else
    disp(['(mean duration) stay and switch palatability coefficients are not significantly different. z = ' num2str(stats.zval) ', p = ' num2str(p)])
end

figure;
notBoxPlot([stayMeanNormalizedCoeffs(:,2) switchMeanNormalizedCoeffs(:,2)])
[p,h,stats] = signrank(stayMeanNormalizedCoeffs(:,2),0);
if (p < .05)
    disp(['(mean duration) stay alternative palatability coefficients are significantly different from 0. z = ' num2str(stats.zval) ', p = ' num2str(p)])
else
    disp(['(mean duration) stay alternative palatability coefficients are not significantly different from 0. z = ' num2str(stats.zval) ', p = ' num2str(p)])
end
[p,h,stats] = signrank(switchMeanNormalizedCoeffs(:,2),0,'tail','left');
if (p < .05)
    disp(['(mean duration) switch alternative palatability coefficients are significantly less than 0. z = ' num2str(stats.zval) ', p = ' num2str(p)])
else
    disp(['(mean duration) switch alternative palatability coefficients are not significantly less than 0. z = ' num2str(stats.zval) ', p = ' num2str(p)])
end
[p,h,stats] = signrank(stayMeanNormalizedCoeffs(:,2),switchMeanNormalizedCoeffs(:,2),'tail','right');
if (p < .05)
    disp(['(mean duration) stay alternative palatability coefficients are significantly greater than switch. z = ' num2str(stats.zval) ', p = ' num2str(p)]);
else
    disp(['(mean duration) stay and switch alternative palatability coefficients are not significantly different. z = ' num2str(stats.zval) ', p = ' num2str(p)]);
end

%% Fig 5: stay vs. switch normalized regression coefficients (need to run prior section)
figure;
notBoxPlot([stayNormalizedBs(:,1); switchNormalizedBs(:,1)],[ones(size(stayNormalizedBs,1),1); ones(size(switchBs,1),1)*2])
ylabel('Normalized palatability coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Stay','Switch'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'fig5_stay_vs_switch_normalized_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end
5784
figure;
notBoxPlot([stayNormalizedBs(:,2); switchNormalizedBs(:,2)],[ones(size(stayNormalizedBs,1),1); ones(size(switchNormalizedBs,1),1)*2])
ylabel({'Normalized alternative','palatability coefficient'},'fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Stay','Switch'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'fig5_stay_vs_switch_normalized_alternative_palatability_coefficients';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

%% Compute transition probabilities
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

%% Get common source/target transition probabilities
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
%sourcePairTitles = {'P_{a -> b} vs. P_{a -> c}','P_{b -> a} vs. P_{b -> c}','P_{c -> a} vs. P_{c -> b}'};
%targetPairTitles = {'P_{b -> a} vs. P_{c -> a}','P_{a -> b} vs. P_{c -> b}','P_{a -> c} vs. P_{b -> c}'};
%sourcePairStrings = {{'P_{a -> b}','P_{a -> c}'},{'P_{b -> a}','P_{b -> c}'},{'P_{c -> a}','P_{c -> b}'}};
%targetPairStrings = {{'P_{b -> a}','P_{c -> a}'},{'P_{a -> b}','P_{c -> b}'},{'P_{a -> c}','P_{b -> c}'}};
sourcePairTitles = {'P_{.01M -> .1M} vs. P_{.01M -> 1M}','P_{.1M -> .01M} vs. P_{.1M -> 1M}','P_{1M -> .01M} vs. P_{1M -> .1M}'};
targetPairTitles = {'P_{.1M -> .01M} vs. P_{1M -> .01M}','P_{.01M -> .1M} vs. P_{1M -> .1M}','P_{.01M -> 1M} vs. P_{.1M -> 1M}'};
sourcePairStrings = {{'P_{.01M -> .1M}','P_{.01M -> 1M}'},{'P_{.1M -> .01M}','P_{.1M -> 1M}'},{'P_{1M -> .01M}','P_{1M -> .1M}'}};
targetPairStrings = {{'P_{.1M -> .01M}','P_{1M -> .01M}'},{'P_{.01M -> .1M}','P_{1M -> .1M}'},{'P_{.01M -> 1M}','P_{.1M -> 1M}'}};
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
    end
end
commonTargetTransVals = cell(length(commonTargetPairInds),2);
for i=1:length(commonTargetPairInds)
    for j=1:length(animalObjs)
        % All bouts
        curTransMat1 = boutTransitionProbs{j}{commonTargetPairInds{i}(1,1),commonTargetPairInds{i}(1,2)};
        curTransMat2 = boutTransitionProbs{j}{commonTargetPairInds{i}(2,1),commonTargetPairInds{i}(2,2)}; 
        commonTargetTransVals{i,1} = [commonTargetTransVals{i,1} curTransMat1(commonTargetPairInds{i}(1,3),commonTargetPairInds{i}(1,4))];
        commonTargetTransVals{i,2} = [commonTargetTransVals{i,2} curTransMat2(commonTargetPairInds{i}(2,3),commonTargetPairInds{i}(2,4))];
    end
end

%% Fig 9A: Common source transition probabilities
figure;
for i=1:length(commonSourcePairInds)
    subplot(1,3,i)
    %{
    [p1,h1] = signrank(commonSourceTransVals{i,1},commonSourceTransVals{i,2},'tail','right');
    [p2,h2] = signrank(commonSourceTransVals{i,1},commonSourceTransVals{i,2},'tail','left');
    [p3,h3] = signrank(commonSourceTransVals{i,1},commonSourceTransVals{i,2});
    p1 = p1*6; % Bonferroni correction
    p2 = p2*6; % Bonferroni correction
    if (p1 < .05)
        curText = [sourcePairStrings{i}{1} ' > ' sourcePairStrings{i}{2}];
    elseif (p2 < .05)
        curText = [sourcePairStrings{i}{1} ' < ' sourcePairStrings{i}{2}];
    else
        curText = [sourcePairStrings{i}{1} ' = ' sourcePairStrings{i}{2}];
    end
    %}
    [pMiddle,~,statsMiddle] = signrank(commonSourceTransVals{i,1},commonSourceTransVals{i,2});
    [pHigher,~,statsHigher] = signrank(commonSourceTransVals{i,1},commonSourceTransVals{i,2},'tail','right');
    [pLower,~,statsLower] = signrank(commonSourceTransVals{i,1},commonSourceTransVals{i,2},'tail','left');
    pMiddle = pMiddle*6;
    pHigher = pHigher*6;
    pLower = pLower*6;
    if (pMiddle < .05)
        if (pHigher < .05)
            disp(['(source) ' sourcePairStrings{i}{1} ' > ' sourcePairStrings{i}{2} ' p = ' num2str(pHigher) ', z = ', num2str(statsHigher.zval)])
        elseif (pLower < .05)
            disp(['(source) ' sourcePairStrings{i}{1} ' < ' sourcePairStrings{i}{2} ' p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
        else
            disp(['(source) ' sourcePairStrings{i}{1} ' is significantly different than ' sourcePairStrings{i}{2} ' p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
        end
    else
        disp(['(source) ' sourcePairStrings{i}{1} ' is not different than ' sourcePairStrings{i}{2} ' p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
    end
    notBoxPlot([commonSourceTransVals{i,1} commonSourceTransVals{i,2}],[ones(1,length(commonSourceTransVals{i,1})) ones(1,length(commonSourceTransVals{i,2}))*2.5]);
    ylim([0 1])
    ylabel('Transition probability','fontsize',15,'fontweight','bold')
    set(gca,'xtick',[1 2.5],'xticklabels',{sourcePairStrings{i}{1},sourcePairStrings{i}{2}},'fontsize',15,'fontweight','bold')
    %title(curText,'fontsize',15,'fontweight','bold')
end
set(gcf,'Position',[10 10 1800 500])
if (do_save)
    figName = 'fig9A_common_source';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

%% Fig 9B: Common target transition probabilities
figure;
for i=1:length(commonTargetPairInds)
    subplot(1,3,i)
    %{
    [p1,h1] = signrank(commonTargetTransVals{i,1},commonTargetTransVals{i,2},'tail','right');
    [p2,h2] = signrank(commonTargetTransVals{i,1},commonTargetTransVals{i,2},'tail','left');
    [p3,h3] = signrank(commonTargetTransVals{i,1},commonTargetTransVals{i,2});
    p1 = p1*6; % Bonferroni correction
    p2 = p2*6; % Bonferroni correction
    if (p1 < .05)
        curText = [targetPairStrings{i}{1} ' > ' targetPairStrings{i}{2}];
    elseif (p2 < .05)
        curText = [targetPairStrings{i}{1} ' < ' targetPairStrings{i}{2}];
    else
        curText = [targetPairStrings{i}{1} ' = ' targetPairStrings{i}{2}];
    end
    %}
    [pMiddle,~,statsMiddle] = signrank(commonTargetTransVals{i,1},commonTargetTransVals{i,2});
    [pHigher,~,statsHigher] = signrank(commonTargetTransVals{i,1},commonTargetTransVals{i,2},'tail','right');
    [pLower,~,statsLower] = signrank(commonTargetTransVals{i,1},commonTargetTransVals{i,2},'tail','left');
    pMiddle = pMiddle*6;
    pHigher = pHigher*6;
    pLower = pLower*6;
    if (pMiddle < .05)
        if (pHigher < .05)
            disp(['(target) ' targetPairStrings{i}{1} ' > ' targetPairStrings{i}{2} ' p = ' num2str(pHigher) ', z = ', num2str(statsHigher.zval)])
        elseif (pLower < .05)
            disp(['(target) ' targetPairStrings{i}{1} ' < ' targetPairStrings{i}{2} ' p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
        else
            disp(['(target) ' targetPairStrings{i}{1} ' is significantly different than ' targetPairStrings{i}{2} ' p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
        end
    else
        disp(['(target) ' targetPairStrings{i}{1} ' is not different than ' targetPairStrings{i}{2} ' p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
    end
    notBoxPlot([commonTargetTransVals{i,1} commonTargetTransVals{i,2}],[ones(1,length(commonTargetTransVals{i,1})) ones(1,length(commonTargetTransVals{i,2}))*2.5]);
    ylim([0 1])
    ylabel('Transition probability','fontsize',15,'fontweight','bold')
    set(gca,'xtick',[1 2.5],'xticklabels',{targetPairStrings{i}{1},targetPairStrings{i}{2}},'fontsize',15,'fontweight','bold')
    %title(curText,'fontsize',15,'fontweight','bold')
end
set(gcf,'Position',[10 10 1800 500])
if (do_save)
    figName = 'fig9B_common_target';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

%% Load data using 200ms bout criteria
clear all; close all;
publicationFolder = '/home/ben/phd/behavior/publication/';
dataFolder = [publicationFolder '/data/'];
figFolder = [publicationFolder '/figures/'];
do_save = 1; % Boolean for whether to save figures or not
animalInfo = load('analyzed_data/animalInfo.mat'); animalInfo=animalInfo.animalInfo;

if (~exist('animalObjs','var'))
    for i=1:length(animalInfo)
        animalObjs(i) = animalData(animalInfo(i).animal,animalInfo(i).dates,animalInfo(i).sex,'analyzed_data');
        disp(['done with ' animalInfo(i).animal])
    end
end

% Create table with all bouts and applicable factors
[boutTable] = createBoutDataTable(animalObjs,'boutType','200');

% get all bouts
allBouts = [];
for i=1:length(animalObjs)
    for j=1:length(animalObjs(i).dates)
        for k=1:length(animalObjs(i).bouts200{j})
            if (~isempty(animalObjs(i).bouts200{j}{k}))
                if (any([animalObjs(i).bouts200{j}{k}.duration] < 0))
                    disp([animalObjs(i).name ' ' animalObjs(i).dates{j} ' has negative duration bout'])
                end
                allBouts = [allBouts animalObjs(i).bouts200{j}{k}];
            end
        end
    end
end

%% Creat indices variables for relevant subsets of bouts
difSolDayInds = find(boutTable.solutionNum ~= boutTable.alternative_solutionNum);
Ainds = find(boutTable.solutionNum == 2);
Binds = find(boutTable.solutionNum == 3);
Cinds = find(boutTable.solutionNum == 4);
isEarly = find(boutTable.early);
isLate = find(boutTable.late);

%% Supp Fig 1: Exponential fits to bout durations using 200ms bout criterion
minDurCutoff = 2;
Adurations = boutTable.duration(Ainds); Adurations = Adurations(Adurations > minDurCutoff);
Bdurations = boutTable.duration(Binds); Bdurations = Bdurations(Bdurations > minDurCutoff);
Cdurations = boutTable.duration(Cinds); Cdurations = Cdurations(Cdurations > minDurCutoff);


BIN_WIDTH_A = 2;
BIN_WIDTH_B = 2;
BIN_WIDTH_C = 1;

[nA,edgesA] = histcounts(Adurations,'binwidth',BIN_WIDTH_A);
[nB,edgesB] = histcounts(Bdurations,'binwidth',BIN_WIDTH_B);
[nC,edgesC] = histcounts(Cdurations,'binwidth',BIN_WIDTH_C);

Azeros = find(nA == 0);
Bzeros = find(nB == 0);
Czeros = find(nC == 0);
firstAzero = Azeros(1);
firstBzero = Bzeros(1);
firstCzero = Czeros(1);


Ax_vals = (minDurCutoff + (BIN_WIDTH_A/2)):BIN_WIDTH_A:(minDurCutoff + (firstAzero - 1)*BIN_WIDTH_A);
Bx_vals = (minDurCutoff + (BIN_WIDTH_B/2)):BIN_WIDTH_B:(minDurCutoff + (firstBzero - 1)*BIN_WIDTH_B);
Cx_vals = (minDurCutoff + (BIN_WIDTH_C/2)):BIN_WIDTH_C:(minDurCutoff + (firstCzero - 1)*BIN_WIDTH_C);

linFitA = polyfit(Ax_vals,log(nA(1:(firstAzero - 1))),1);
linFitB = polyfit(Bx_vals,log(nB(1:(firstBzero - 1))),1);
linFitC = polyfit(Cx_vals,log(nC(1:(firstCzero - 1))),1);

figure;
plot(Ax_vals,nA(1:(firstAzero - 1))); set(gca,'Yscale','log')
hold on;
plot(0:max(Ax_vals),exp(linFitA(1)*(0:max(Ax_vals)) + linFitA(2)))
ylabel('No. of occurences (log-scale)','fontweight','bold','fontsize',20)
xlabel('Bout duration (sec)','fontweight','bold','fontsize',20)
legend({'Data','Linear model'},'fontweight','bold','fontsize',20)
xlim([0 firstAzero*BIN_WIDTH_A])
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'suppFig1_exponential_fit_solnA_200ms';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

figure;
plot(Bx_vals,nB(1:(firstBzero - 1))); set(gca,'Yscale','log')
hold on;
plot(0:max(Bx_vals),exp(linFitB(1)*(0:max(Bx_vals)) + linFitB(2)))
ylabel('No. of occurences (log-scale)','fontweight','bold','fontsize',20)
xlabel('Bout duration (sec)','fontweight','bold','fontsize',20)
xlim([0 firstBzero*BIN_WIDTH_B])
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'suppFig1_exponential_fit_solnB_200ms';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

figure;
plot(Cx_vals,nC(1:(firstCzero - 1))); set(gca,'Yscale','log')
hold on;
plot(0:max(Cx_vals),exp(linFitC(1)*(0:max(Cx_vals)) + linFitC(2)))
ylabel('No. of occurences (log-scale)','fontweight','bold','fontsize',20)
xlabel('Bout duration (sec)','fontweight','bold','fontsize',20)
xlim([0 firstCzero*BIN_WIDTH_C])
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'suppFig1_exponential_fit_solnC_200ms';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

%% Do multilinear regressions for bout duration with 200ms bout criterion
difSolInds = find(boutTable.solutionNum ~= boutTable.alternative_solutionNum);
for i=1:length(animalObjs)
    ratInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    curInds = intersect(difSolInds,ratInds);
    X = [boutTable.palatability(curInds) boutTable.alternative_palatability(curInds) ones(length(curInds),1)];
    Y = boutTable.duration(curInds);
    
    [b,bint,r,rint,stats] = regress(Y,X);
    multiLinearBoutDurFits(i).b = b;
    multiLinearBoutDurFits(i).normalizedB = b./mean(Y);
    multiLinearBoutDurFits(i).bint = bint;
    multiLinearBoutDurFits(i).r = r;
    multiLinearBoutDurFits(i).rint = rint;
    multiLinearBoutDurFits(i).stats = stats;
end
z = [multiLinearBoutDurFits.b];
normalizedZ = [multiLinearBoutDurFits.normalizedB];

%% Get R-squared for all bout regressions (200ms criterion)
clear allX allY Rsquared pvalues coeffs normalizedCoeffs
allX = [];
allY = [];
for i=1:length(animalObjs)
    ratInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    curInds = intersect(ratInds,difSolDayInds);
    [uniqueRows,ia,ic] = unique([boutTable.palatability(curInds) boutTable.alternative_palatability(curInds)],'rows');
    curX = [];
    curY = [];
    for j = 1:length(ia)
        curTableInds = curInds(find(ic == j));
        curX = [curX; uniqueRows(j,:) 1];
        curY = [curY; mean(boutTable.duration(curTableInds))];
    end
    [b,bint,r,rint,stats] = regress(curY,curX);
    Rsquared(i) = stats(1);
    pvalues(i) = stats(3);
    coeffs(i,:) = b;
    normalizedCoeffs(i,:) = b./mean(curY);
    
    allX = [allX; curX];
    allY = [allY; curY];
end

palatabilityRange = 0:.001:2;
[b,bint,r,rint,stats] = regress(allY,allX);
plane_allRats = palatabilityRange*b(1) + palatabilityRange'*b(2) + b(3);
figure;
scatter3(allX(:,1),allX(:,2),allY,'.')
hold on;
surf(palatabilityRange,palatabilityRange,plane_allRats,'FaceColor','g','EdgeColor','g')
xlabel('Current'); ylabel('Alternative'); zlabel('Duration (sec)')

figure;
notBoxPlot([normalizedCoeffs(:,1) normalizedCoeffs(:,2)])
[p,h,stats] = signrank(normalizedCoeffs(:,1),0,'tail','right');
if (p < .05)
    disp(['(mean duration) palatability coefficients are significantly greater than 0. z = ' num2str(stats.zval) ', p = ' num2str(p)])
end
[p,h,stats] = signrank(normalizedCoeffs(:,2),0,'tail','left');
if (p < .05)
    disp(['(mean duration) alternative palatability coefficients are significantly less than 0. z = ' num2str(stats.zval) ', p = ' num2str(p)])
end

%% Fig 6A: normalized regression coefficients for bouts with 200ms criterion
figure;
notBoxPlot([normalizedZ(1,:) normalizedZ(2,:)],[ones(1,length(normalizedZ(1,:))) ones(1,length(normalizedZ(2,:)))*2])
ylabel('Normalized regression coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Palatability','Alternative palatability'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'fig6A_all_bouts_normalized_palatability_coefficients_200ms';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'eps')
    print([figFolder figName],'-dpng','-r600')
end

[pHigher,~,statsHigher] = signrank(normalizedZ(1,:),0,'tail','right');
[pLower,~,statsLower] = signrank(normalizedZ(1,:),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(normalizedZ(1,:),0);
if (pMiddle < .05)
    if (pHigher < .05)
        disp(['Normalized palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
    elseif (pLower < .05)
        disp(['Normalized palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
    else
        disp(['Normalized palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
    end
else
    disp(['Normalized palatability coefficients are not significantly different from 0. p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
end

[pHigher,~,statsHigher] = signrank(normalizedZ(2,:),0,'tail','right');
[pLower,~,statsLower] = signrank(normalizedZ(2,:),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(normalizedZ(2,:),0);
if (pMiddle < .05)
    if (pHigher < .05)
        disp(['Normalized alternative palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
    elseif (pLower < .05)
        disp(['Normalized alternative palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
    else
        disp(['Normalized alternative palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
    end
else
    disp(['Normalized alternative palatability coefficients are not significanlty different from 0. p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
end

%% Do multilinear regressions for stay vs. switch bouts with 200ms criterion
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
difSolInds = find(boutTable.solutionNum ~= boutTable.alternative_solutionNum);
for i=1:length(animalObjs)
    count = count+1;
    switchX = [];
    switchY = [];
    stayX = [];
    stayY = [];
    ratInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    for j=animalObjs(i).difSolutionDays
        dayInds = find(strcmp(boutTable.date,animalObjs(i).dates{j}));
        curInds = intersect(ratInds,dayInds);
        curBoutNums = boutTable.boutNumber(curInds);
        [~,sortedOrder] = sort(curBoutNums);
        sortedInds = curInds(sortedOrder);
        for k=2:length(sortedInds)
            curPal = boutTable.palatability(sortedInds(k));
            curDur = boutTable.duration(sortedInds(k));
            if (boutTable.channel(sortedInds(k)) ~= boutTable.channel(sortedInds(k-1)))
                switchCount = switchCount+1;
                switchRats{switchCount} = animalObjs(i).name;
                curAltPal = boutTable.alternative_palatability(sortedInds(k));
                switchDurs = [switchDurs curDur];
                switchAltPals = [switchAltPals curAltPal];
                switchPals = [switchPals curPal];
                switchX = [switchX; [curPal curAltPal]];
                switchY = [switchY curDur];
            else
                stayCount = stayCount + 1;
                stayRats{stayCount} = animalObjs(i).name;
                curAltPal = boutTable.alternative_palatability(sortedInds(k));
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
    
    [uniqueRows,ia,ic] = unique(switchX,'rows');
    meanSwitchX = [];
    meanSwitchY = [];
    for j=1:length(ia)
        curXinds = find(ic == j);
        meanSwitchX = [meanSwitchX; uniqueRows(j,:) 1];
        meanSwitchY = [meanSwitchY; mean(switchY(curXinds))];
    end
    [b,bint,r,rint,stats] = regress(meanSwitchY,meanSwitchX);
    switchRsquared(i) = stats(1);
    switchPvalues(i) = stats(3);
    switchMeanCoeffs(i,:) = b;
    switchMeanNormalizedCoeffs(i,:) = b./mean(meanSwitchY);
    
    stayX = [stayX ones(size(stayX,1),1)];
    [b,bint,r,rint,stats] = regress(stayY',stayX);
    stayBs(count,:) = b;
    stayNormalizedBs(count,:) = b./mean(stayY);
    stayBINTs{count} = bint;
    stayRs{count} = r;
    stayRINTs{count} = rint;
    staySTATs{count} = stats;
    
    [uniqueRows,ia,ic] = unique(stayX,'rows');
    meanStayX = [];
    meanStayY = [];
    for j=1:length(ia)
        curXinds = find(ic == j);
        meanStayX = [meanStayX; uniqueRows(j,:) 1];
        meanStayY = [meanStayY; mean(stayY(curXinds))];
    end
    [b,bint,r,rint,stats] = regress(meanStayY,meanStayX);
    stayRsquared(i) = stats(1);
    stayPvalues(i) = stats(3);
    stayMeanCoeffs(i,:) = b;
    stayMeanNormalizedCoeffs(i,:) = b./mean(meanStayY);
end

figure;
notBoxPlot([stayMeanNormalizedCoeffs(:,1) switchMeanNormalizedCoeffs(:,1)])
[p,h,stats] = signrank(stayMeanNormalizedCoeffs(:,1),switchMeanNormalizedCoeffs(:,1));
if (p < .05)
    disp(['(mean duration) stay and switch palatability coefficients are significantly different. z = ' num2str(stats.zval) ', p = ' num2str(p)])
else
    disp(['(mean duration) stay and switch palatability coefficients are not significantly different. z = ' num2str(stats.zval) ', p = ' num2str(p)])
end

figure;
notBoxPlot([stayMeanNormalizedCoeffs(:,2) switchMeanNormalizedCoeffs(:,2)])
[p,h,stats] = signrank(stayMeanNormalizedCoeffs(:,2),0);
if (p < .05)
    disp(['(mean duration) stay alternative palatability coefficients are significantly different from 0. z = ' num2str(stats.zval) ', p = ' num2str(p)])
else
    disp(['(mean duration) stay alternative palatability coefficients are not significantly different from 0. z = ' num2str(stats.zval) ', p = ' num2str(p)])
end
[p,h,stats] = signrank(switchMeanNormalizedCoeffs(:,2),0,'tail','left');
if (p < .05)
    disp(['(mean duration) switch alternative palatability coefficients are significantly less than 0. z = ' num2str(stats.zval) ', p = ' num2str(p)])
else
    disp(['(mean duration) switch alternative palatability coefficients are not significantly less than 0. z = ' num2str(stats.zval) ', p = ' num2str(p)])
end
[p,h,stats] = signrank(stayMeanNormalizedCoeffs(:,2),switchMeanNormalizedCoeffs(:,2),'tail','right');
if (p < .05)
    disp(['(mean duration) stay alternative palatability coefficients are significantly greater than switch. z = ' num2str(stats.zval) ', p = ' num2str(p)]);
else
    disp(['(mean duration) stay and switch alternative palatability coefficients are not significantly different. z = ' num2str(stats.zval) ', p = ' num2str(p)]);
end

%% Fig 6C: normalized regression coefficients for stay vs. switch bouts with 200ms criterion
figure;
notBoxPlot([stayNormalizedBs(:,1); switchNormalizedBs(:,1)],[ones(size(stayNormalizedBs,1),1); ones(size(switchBs,1),1)*2])
ylabel('Normalized palatability coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Stay','Switch'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'fig6C_stay_vs_switch_normalized_palatability_coefficients_200ms';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end
[pHigher,~,statsHigher] = signrank(stayNormalizedBs(:,1),0,'tail','right');
[pLower,~,statsLower] = signrank(stayNormalizedBs(:,1),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(stayNormalizedBs(:,1),0);
if (pHigher < .05)
    disp(['Stay normalized palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
elseif (pLower < .05)
    disp(['Stay normalized palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
elseif (pMiddle < .05)
    disp(['Stay normalized palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
else
    disp(['Stay normalized paltability coefficients are not significantly different from 0. p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
end

[pHigher,~,statsHigher] = signrank(switchNormalizedBs(:,1),0,'tail','right');
[pLower,~,statsLower] = signrank(switchNormalizedBs(:,1),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(switchNormalizedBs(:,1),0);
if (pHigher < .05)
    disp(['Switch normalized palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
elseif (pLower < .05)
    disp(['Switch normalized palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
elseif (pMiddle < .05)
    disp(['Switch normalized palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
else
    disp(['Switch normalized paltability coefficients are not significantly different from 0. p = ' num2str(p) ', z = ' num2str(statsMiddle.zval)])
end

[pMiddle,~,statsMiddle] = signrank(stayNormalizedBs(:,1),switchNormalizedBs(:,1));
if (pMiddle < .05)
    [pHigher,~,statsHigher] = signrank(stayNormalizedBs(:,1),switchNormalizedBs(:,1),'tail','right');
    [pLower,~,statsLower] = signrank(stayNormalizedBs(:,1),switchNormalizedBs(:,1),'tail','left');
    if (pHigher < .05)
        disp(['Stay normalized palatability coefficients are greater than those for switches. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
    elseif (pLower < .05)
        disp(['Stay normalized palatability coefficients are significantly less than those for switches. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
    else
        disp(['Stay normalized palatability coefficients are significantly different from those for switches. (two-tailed) p = ' num2str(pMiddle) ', z =  ' num2str(statsMiddle.zval)])
    end
else
    disp(['Stay normalized palatability coefficients are not significantly different from those for switches. p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
end

figure;
notBoxPlot([stayNormalizedBs(:,2); switchNormalizedBs(:,2)],[ones(size(stayNormalizedBs,1),1); ones(size(switchNormalizedBs,1),1)*2])
ylabel({'Normalized alternative','palatability coefficient'},'fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Stay','Switch'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'fig6C_stay_vs_switch_normalized_alternative_palatability_coefficients_200ms';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end
[pHigher,~,statsHigher] = signrank(stayNormalizedBs(:,2),0,'tail','right');
[pLower,~,statsLower] = signrank(stayNormalizedBs(:,2),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(stayNormalizedBs(:,2),0);
if (pHigher < .05)
    disp(['Stay normalized alternative palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
elseif (pLower < .05)
    disp(['Stay normalized alternative palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
elseif (pMiddle < .05)
    disp(['Stay normalized alternative palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
else
    disp(['Stay normalized alternative paltability coefficients are not significantly different from 0. p = ' num2str(p) ', z = ' num2str(statsMiddle.zval)])
end

[pHigher,~,statsHigher] = signrank(switchNormalizedBs(:,2),0,'tail','right');
[pLower,~,statsLower] = signrank(switchNormalizedBs(:,2),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(switchNormalizedBs(:,2),0);
if (pHigher < .05)
    disp(['Switch normalized alternative palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
elseif (pLower < .05)
    disp(['Switch normalized alternative palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
elseif (pMiddle < .05)
    disp(['Switch normalized alternative palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
else
    disp(['Switch normalized alternative paltability coefficients are not significantly different from 0. p = ' num2str(p)])
end

[pMiddle,~,statsMiddle] = signrank(stayNormalizedBs(:,2),switchNormalizedBs(:,2));
if (pMiddle < .05)
    [pHigher,~,statsHigher] = signrank(stayNormalizedBs(:,2),switchNormalizedBs(:,2),'tail','right');
    [pLower,~,statsLower] = signrank(stayNormalizedBs(:,2),switchNormalizedBs(:,2),'tail','left');
    if (pHigher < .05)
        disp(['Stay normalized alternative palatability coefficients are greater than those for switches. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
    elseif (pLower < .05)
        disp(['Stay normalized alternative palatability coefficients are significantly less than those for switches. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
    else
        disp(['Stay normalized alternative palatability coefficients are significantly different from those for switches. (two-tailed) p = ' num2str(p) ', z = ' num2str(statsMiddle.zval)])
    end
end

%% Do multilinear regressions for early vs. late bouts with 200ms criterion
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
    
    [uniqueRows,ia,ic] = unique(earlyX,'rows');
    meanEarlyX = [];
    meanEarlyY = [];
    for j=1:length(ia)
        curXinds = find(ic == j);
        meanEarlyX = [meanEarlyX; uniqueRows(j,:) 1];
        meanEarlyY = [meanEarlyY; mean(earlyY(curXinds))];
    end
    [b,bint,r,rint,stats] = regress(meanEarlyY,meanEarlyX);
    earlyRsquared(i) = stats(1);
    earlyPvalues(i) = stats(3);
    earlyMeanCoeffs(i,:) = b;
    earlyMeanNormalizedCoeffs(i,:) = b./mean(meanEarlyY);
    
    [uniqueRows,ia,ic] = unique(lateX,'rows');
    meanLateX = [];
    meanLateY = [];
    for j=1:length(ia)
        curXinds = find(ic == j);
        meanLateX = [meanLateX; uniqueRows(j,:) 1];
        meanLateY = [meanLateY; mean(lateY(curXinds))];
    end
    [b,bint,r,rint,stats] = regress(meanLateY,meanLateX);
    lateRsquared(i) = stats(1);
    latePvalues(i) = stats(3);
    lateMeanCoeffs(i,:) = b;
    lateMeanNormalizedCoeffs(i,:) = b./mean(meanLateY);
end

%% Fig 7: Normalized regression coefficients for early vs. late bouts with 200ms bout criterion
figure;
notBoxPlot([earlyNormalizedBs(:,1); lateNormalizedBs(:,1)],[ones(size(earlyNormalizedBs,1),1); ones(size(lateNormalizedBs,1),1)*2])
ylabel('Normalized palatability coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'early','late'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'fig7_early_vs_late_normalized_palatability_coefficients_bouts200';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

[pHigher,~,statsHigher] = signrank(earlyNormalizedBs(:,1),0,'tail','right');
[pLower,~,statsLower] = signrank(earlyNormalizedBs(:,1),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(earlyNormalizedBs(:,1),0);
if (pHigher < .05)
    disp(['Early normalized palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
elseif (pLower < .05)
    disp(['Early normalized palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
elseif (pMiddle < .05)
    disp(['Early normalized palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
else
    disp(['Early normalized paltability coefficients are not significantly different from 0. p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
end

[pHigher,~,statsHigher] = signrank(lateNormalizedBs(:,1),0,'tail','right');
[pLower,~,statsLower] = signrank(lateNormalizedBs(:,1),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(lateNormalizedBs(:,1),0);
if (pHigher < .05)
    disp(['Late normalized palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
elseif (pLower < .05)
    disp(['Late normalized palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
elseif (pMiddle < .05)
    disp(['Late normalized palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
else
    disp(['Late normalized paltability coefficients are not significantly different from 0. p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
end

[pHigher,~,statsHigher] = signrank(earlyNormalizedBs(:,1),lateNormalizedBs(:,1),'tail','right');
[pLower,~,statsLower] = signrank(earlyNormalizedBs(:,1),lateNormalizedBs(:,1),'tail','left');
[pMiddle,~,statsMiddle] = signrank(earlyNormalizedBs(:,1),lateNormalizedBs(:,1));
if (pHigher < .05)
    disp(['Early normalized palatability coefficients are greater than late. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
elseif (pLower < .05)
    disp(['Early normalized palatability coefficients are less than late. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
elseif (pMiddle < .05)
    disp(['Early normalized palatability coefficients are significantly different from late. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
else
    disp(['Early normalized paltability coefficients are not significantly different from late. p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
end

figure;
notBoxPlot([earlyNormalizedBs(:,2); lateNormalizedBs(:,2)],[ones(size(earlyNormalizedBs,1),1); ones(size(lateNormalizedBs,1),1)*2])
ylabel({'Normalized Alternative','palatability coefficient'},'fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'early','late'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'fig7_early_vs_late_normalized alternative_palatability_coefficients_bouts200';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

[pHigher,~,statsHigher] = signrank(earlyNormalizedBs(:,2),0,'tail','right');
[pLower,~,statsLower] = signrank(earlyNormalizedBs(:,2),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(earlyNormalizedBs(:,2),0);
if (pHigher < .05)
    disp(['Early normalized palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
elseif (pLower < .05)
    disp(['Early normalized palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
elseif (pMiddle < .05)
    disp(['Early normalized palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
else
    disp(['Early normalized paltability coefficients are not significantly different from 0. p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
end

[pHigher,~,statsHigher] = signrank(lateNormalizedBs(:,2),0,'tail','right');
[pLower,~,statsLower] = signrank(lateNormalizedBs(:,2),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(lateNormalizedBs(:,2),0);
if (pHigher < .05)
    disp(['Late normalized palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
elseif (pLower < .05)
    disp(['Late normalized palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
elseif (pMiddle < .05)
    disp(['Late normalized palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
else
    disp(['Late normalized paltability coefficients are not significantly different from 0. p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
end

[pHigher,~,statsHigher] = signrank(earlyNormalizedBs(:,2),lateNormalizedBs(:,2),'tail','right');
[pLower,~,statsLower] = signrank(earlyNormalizedBs(:,2),lateNormalizedBs(:,2),'tail','left');
[pMiddle,~,statsMiddle] = signrank(earlyNormalizedBs(:,2),lateNormalizedBs(:,2));
if (pHigher < .05)
    disp(['Early normalized palatability coefficients are greater than late. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
elseif (pLower < .05)
    disp(['Early normalized palatability coefficients are less than late. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
elseif (pMiddle < .05)
    disp(['Early normalized palatability coefficients are significantly different from late. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
else
    disp(['Early normalized paltability coefficients are not significantly different from late. p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
end

%% Load data with altered palatability (Pal(.01M) = 1.3x Pal(.1M))
clear all; close all;
publicationFolder = '/home/ben/phd/behavior/publication/';
dataFolder = [publicationFolder '/data/'];
figFolder = [publicationFolder '/figures/'];

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

%% Creat indices variables for relevant subsets of bouts
difSolDayInds = find(boutTable.solutionNum ~= boutTable.alternative_solutionNum);
Ainds = find(boutTable.solutionNum == 2);
Binds = find(boutTable.solutionNum == 3);
Cinds = find(boutTable.solutionNum == 4);
isEarly = find(boutTable.early);
isLate = find(boutTable.late);

%% Do multilinear regressions for bout duration using altered palatablity
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
end
z = [multiLinearBoutDurFits.b];
normalizedZ = [multiLinearBoutDurFits.normalizedB];

%% Get R-squared for all bout regressions (altered palatability)
clear allX allY Rsquared pvalues coeffs normalizedCoeffs
allX = [];
allY = [];
for i=1:length(animalObjs)
    ratInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    curInds = intersect(ratInds,difSolDayInds);
    [uniqueRows,ia,ic] = unique([boutTable.palatability(curInds) boutTable.alternative_palatability(curInds)],'rows');
    curX = [];
    curY = [];
    for j = 1:length(ia)
        curTableInds = curInds(find(ic == j));
        curX = [curX; uniqueRows(j,:) 1];
        curY = [curY; mean(boutTable.duration(curTableInds))];
    end
    [b,bint,r,rint,stats] = regress(curY,curX);
    Rsquared(i) = stats(1);
    pvalues(i) = stats(3);
    coeffs(i,:) = b;
    normalizedCoeffs(i,:) = b./mean(curY);
    
    allX = [allX; curX];
    allY = [allY; curY];
end

palatabilityRange = 0:.001:2;
[b,bint,r,rint,stats] = regress(allY,allX);
plane_allRats = palatabilityRange*b(1) + palatabilityRange'*b(2) + b(3);
figure;
scatter3(allX(:,1),allX(:,2),allY,'.')
hold on;
surf(palatabilityRange,palatabilityRange,plane_allRats,'FaceColor','g','EdgeColor','g')
xlabel('Current'); ylabel('Alternative'); zlabel('Duration (sec)')

figure;
notBoxPlot([normalizedCoeffs(:,1) normalizedCoeffs(:,2)])
[p,h,stats] = signrank(normalizedCoeffs(:,1),0,'tail','right');
if (p < .05)
    disp(['(mean duration) palatability coefficients are significantly greater than 0. z = ' num2str(stats.zval) ', p = ' num2str(p)])
end
[p,h,stats] = signrank(normalizedCoeffs(:,2),0,'tail','left');
if (p < .05)
    disp(['(mean duration) alternative palatability coefficients are significantly less than 0. z = ' num2str(stats.zval) ', p = ' num2str(p)])
end

%% Fig 8A: Normalized regression coefficients for altered palatability
figure;
notBoxPlot([normalizedZ(1,:) normalizedZ(2,:)],[ones(1,length(normalizedZ(1,:))) ones(1,length(normalizedZ(2,:)))*2])
ylabel('Normalized regression coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Palatability','Alternative palatability'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'fig7A_all_bouts_normalized_palatability_coefficients_alteredPalatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

[pHigher,~,statsHigher] = signrank(normalizedZ(1,:),0,'tail','right');
[pLower,~,statsLower] = signrank(normalizedZ(1,:),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(normalizedZ(1,:),0);
if (pMiddle < .05)
    if (pHigher < .05)
        disp(['Normalized palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
    elseif (pLower < .05)
        disp(['Normalized palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
    else
        disp(['Normalized palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
    end
else
    disp(['Normalized palatability coefficients are not significantly different from 0. p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
end

[pHigher,~,statsHigher] = signrank(normalizedZ(2,:),0,'tail','right');
[pLower,~,statsLower] = signrank(normalizedZ(2,:),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(normalizedZ(2,:),0);
if (pMiddle < .05)
    if (pHigher < .05)
        disp(['Normalized alternative palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
    elseif (pLower < .05)
        disp(['Normalized alternative palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
    else
        disp(['Normalized alternative palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
    end
else
    disp(['Normalized alternative palatability coefficients are not significanlty different from 0. p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
end

%% Do multilinear regressions for stay vs. switch bouts with alterered palatability
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
    
    [uniqueRows,ia,ic] = unique(switchX,'rows');
    meanSwitchX = [];
    meanSwitchY = [];
    for j=1:length(ia)
        curXinds = find(ic == j);
        meanSwitchX = [meanSwitchX; uniqueRows(j,:) 1];
        meanSwitchY = [meanSwitchY; mean(switchY(curXinds))];
    end
    [b,bint,r,rint,stats] = regress(meanSwitchY,meanSwitchX);
    switchRsquared(i) = stats(1);
    switchPvalues(i) = stats(3);
    switchMeanCoeffs(i,:) = b;
    switchMeanNormalizedCoeffs(i,:) = b./mean(meanSwitchY);
    
    stayX = [stayX ones(size(stayX,1),1)];
    [b,bint,r,rint,stats] = regress(stayY',stayX);
    stayBs(count,:) = b;
    stayNormalizedBs(count,:) = b./mean(stayY);
    stayBINTs{count} = bint;
    stayRs{count} = r;
    stayRINTs{count} = rint;
    staySTATs{count} = stats;
   
    [uniqueRows,ia,ic] = unique(stayX,'rows');
    meanStayX = [];
    meanStayY = [];
    for j=1:length(ia)
        curXinds = find(ic == j);
        meanStayX = [meanStayX; uniqueRows(j,:) 1];
        meanStayY = [meanStayY; mean(stayY(curXinds))];
    end
    [b,bint,r,rint,stats] = regress(meanStayY,meanStayX);
    stayRsquared(i) = stats(1);
    stayPvalues(i) = stats(3);
    stayMeanCoeffs(i,:) = b;
    stayMeanNormalizedCoeffs(i,:) = b./mean(meanStayY);
end


%% Fig 8B: Normalized regression coefficients for stay vs. switch bouts with altered palatability
figure;
notBoxPlot([stayNormalizedBs(:,1); switchNormalizedBs(:,1)],[ones(size(stayNormalizedBs,1),1); ones(size(switchBs,1),1)*2])
ylabel('Normalized palatability coefficient','fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Stay','Switch'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'fig8B_stay_vs_switch_normalized_palatability_coefficients_alteredPalatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

[pHigher,~,statsHigher] = signrank(stayNormalizedBs(:,1),0,'tail','right');
[pLower,~,statsLower] = signrank(stayNormalizedBs(:,1),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(stayNormalizedBs(:,1),0);
if (pHigher < .05)
    disp(['Stay normalized palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
elseif (pLower < .05)
    disp(['Stay normalized palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
elseif (pMiddle < .05)
    disp(['Stay normalized palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
else
    disp(['Stay normalized paltability coefficients are not significantly different from 0. p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
end

[pHigher,~,statsHigher] = signrank(switchNormalizedBs(:,1),0,'tail','right');
[pLower,~,statsLower] = signrank(switchNormalizedBs(:,1),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(switchNormalizedBs(:,1),0);
if (pHigher < .05)
    disp(['Switch normalized palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
elseif (pLower < .05)
    disp(['Switch normalized palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
elseif (pMiddle < .05)
    disp(['Switch normalized palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
else
    disp(['Switch normalized paltability coefficients are not significantly different from 0. p = ' num2str(p) ', z = ' num2str(statsMiddle.zval)])
end

[pMiddle,~,statsMiddle] = signrank(stayNormalizedBs(:,1),switchNormalizedBs(:,1));
if (pMiddle < .05)
    [pHigher,~,statsHigher] = signrank(stayNormalizedBs(:,1),switchNormalizedBs(:,1),'tail','right');
    [pLower,~,statsLower] = signrank(stayNormalizedBs(:,1),switchNormalizedBs(:,1),'tail','left');
    if (pHigher < .05)
        disp(['Stay normalized palatability coefficients are greater than those for switches. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
    elseif (pLower < .05)
        disp(['Stay normalized palatability coefficients are significantly less than those for switches. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
    else
        disp(['Stay normalized palatability coefficients are significantly different from those for switches. (two-tailed) p = ' num2str(pMiddle) ', z =  ' num2str(statsMiddle.zval)])
    end
else
    disp(['Stay normalized palatability coefficients are not significantly different from those for switches. p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
end

figure;
notBoxPlot([stayNormalizedBs(:,2); switchNormalizedBs(:,2)],[ones(size(stayNormalizedBs,1),1); ones(size(switchNormalizedBs,1),1)*2])
ylabel({'Normalized alternative','palatability coefficient'},'fontsize',15,'fontweight','bold')
set(gca,'xtick',[1 2],'xticklabels',{'Stay','Switch'},'fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1200 1200])
if (do_save)
    figName = 'fig8B_stay_vs_switch_normalized_alternative_palatability_coefficients_alteredPalatability';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end
[pHigher,~,statsHigher] = signrank(stayNormalizedBs(:,2),0,'tail','right');
[pLower,~,statsLower] = signrank(stayNormalizedBs(:,2),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(stayNormalizedBs(:,2),0);
if (pHigher < .05)
    disp(['Stay normalized alternative palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
elseif (pLower < .05)
    disp(['Stay normalized alternative palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
elseif (pMiddle < .05)
    disp(['Stay normalized alternative palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
else
    disp(['Stay normalized alternative paltability coefficients are not significantly different from 0. p = ' num2str(p) ', z = ' num2str(statsMiddle.zval)])
end

[pHigher,~,statsHigher] = signrank(switchNormalizedBs(:,2),0,'tail','right');
[pLower,~,statsLower] = signrank(switchNormalizedBs(:,2),0,'tail','left');
[pMiddle,~,statsMiddle] = signrank(switchNormalizedBs(:,2),0);
if (pHigher < .05)
    disp(['Switch normalized alternative palatability coefficients are greater than 0. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
elseif (pLower < .05)
    disp(['Switch normalized alternative palatability coefficients are less than 0. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
elseif (pMiddle < .05)
    disp(['Switch normalized alternative palatability coefficients are significantly different from 0. (two-tailed) p = ' num2str(pMiddle) ', z = ' num2str(statsMiddle.zval)])
else
    disp(['Switch normalized alternative paltability coefficients are not significantly different from 0. p = ' num2str(p)])
end

[pMiddle,~,statsMiddle] = signrank(stayNormalizedBs(:,2),switchNormalizedBs(:,2));
if (pMiddle < .05)
    [pHigher,~,statsHigher] = signrank(stayNormalizedBs(:,2),switchNormalizedBs(:,2),'tail','right');
    [pLower,~,statsLower] = signrank(stayNormalizedBs(:,2),switchNormalizedBs(:,2),'tail','left');
    if (pHigher < .05)
        disp(['Stay normalized alternative palatability coefficients are greater than those for switches. (right-tailed) p = ' num2str(pHigher) ', z = ' num2str(statsHigher.zval)])
    elseif (pLower < .05)
        disp(['Stay normalized alternative palatability coefficients are significantly less than those for switches. (left-tailed) p = ' num2str(pLower) ', z = ' num2str(statsLower.zval)])
    else
        disp(['Stay normalized alternative palatability coefficients are significantly different from those for switches. (two-tailed) p = ' num2str(p) ', z = ' num2str(statsMiddle.zval)])
    end
end

%% Extract data about prior day preferred side choices
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

%% Fig 10A: Histogram of # of days rats first chose the prior day's preferred side
figure;
histogram(nPriorPreferredFrac,'binwidth',.05)
xlabel('Probability of choosing prior preferred side','fontsize',15,'fontweight','bold')
ylabel('# of animals','fontsize',15,'fontweight','bold')
set(gcf,'Position',[10 10 1600 1200])
if (do_save)
    figName = 'fig10A_probability_of_choosing_prior_day_preferred_side';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

%% Fig 10B: Comparison of mouse prior day preferred selection to bionomial distribution
figure;
plot(0:totalBinomTestCount,binomPDF); hold on; plot([sum(nPriorPreferred) sum(nPriorPreferred)],[0 .1],'r')
xlabel('# of prior preferred chosen','fontsize',15,'fontweight','bold');
ylabel('Probability density','fontsize',15,'fontweight','bold')
if (do_save)
    figName = 'fig10B_prior_day_preferred_binomial_test';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

%% Fig ?: Current solution palatability and alternative solution palatability coefficients as a function of time since alternative visit
difSolDayInds = find(boutTable.solutionNum ~= boutTable.alternative_solutionNum);
allX = cell(1,1);
allY = cell(1,1);
X = cell(1,1);
Y = cell(1,1);
for i=1:length(animalObjs)
    allX{i} = [];
    allY{i} = [];
    ratInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    for j=1:length(animalObjs(i).difSolutionDays)
        dayInds = find(strcmp(boutTable.date,animalObjs(i).dates{animalObjs(i).difSolutionDays(j)}));
        curBouts = boutTable(intersect(ratInds,dayInds),:);
        [~,sortedInds] = sort(curBouts.onset);
        curBouts = curBouts(sortedInds,:);
        if (height(curBouts) < 2)
            continue;
        else
            for k=2:height(curBouts)
                pastBouts = curBouts(1:(k-1),:);
                altChannelBoutInds = find(pastBouts.channel ~= curBouts.channel(k));
                if (isempty(altChannelBoutInds))
                    continue;
                else
                    lastAltVisitBout = pastBouts(max(altChannelBoutInds),:);
                    timeSinceLastAltVisit = curBouts.onset(k) - lastAltVisitBout.offset;
                    allX{i} = [allX{i}; curBouts.palatability(k) lastAltVisitBout.palatability timeSinceLastAltVisit];
                    allY{i} = [allY{i}; curBouts.duration(k)];
                end
            end
        end
    end
end


%stdErr = std(altSolCoeffs,1,'omitnan')./sqrt(sum(~isnan(altSolCoeffs),1));
%figure;
%shadedErrorBar(bins,mean(altSolCoeffs,1,'omitnan'),stdErr)
figure;
notBoxPlot(altSolCoeffs)
ylabel({'Normalized alternative','palatability coefficients'})
set(gca,'xticklabels',{['\Delta t_{last visit} < ' num2str(bins(1))],['\Delta t_{last visit} < ' num2str(bins(2))]})
set(gcf,'Position',[10 10 1200 600])
[p,h] = ranksum(altSolCoeffs(:,1),altSolCoeffs(:,2))
if (do_save)
    figName = 'fig11A_normalized_alt_coefficients_by_time_since_last_visit';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

minbins = 10:1:60;
pvals = nan(1,length(minbins));
for t=1:length(minbins)
    bins = [minbins(t) 6000];
    normalizedCoeffs = nan(length(bins),3);
    altSolCoeffs = nan(length(animalObjs),length(bins));
    for i=1:length(animalObjs)
        for j = 1:length(bins)
            if (j == 1)
                inds = find(allX{i}(:,3) < bins(j));
            else
                inds = find(allX{i}(:,3) < bins(j) & allX{i}(:,3) > bins(j-1));
            end
            X{j} = [allX{i}(inds,1:2) ones(length(inds),1)];
            Y{j} = allY{i}(inds);
            if (size(X{j},1) < 5)
                continue;
            else
                [b,bint,r,rint,stats] = regress(Y{j},X{j});
                normalizedCoeffs(j,:) = b./mean(Y{j});
                altSolCoeffs(i,j) = normalizedCoeffs(j,2);
            end
        end
    end
    [p,h] = ranksum(altSolCoeffs(:,1),altSolCoeffs(:,2));
    pvals(t) = p;
end

figure;
plot(minbins,pvals); xlabel(['\Delta t_{last visit} threshold'])
ylabel('Wilcoxon Rank-Sum test p-value')
set(gcf,'Position',[10 10 1200 800])
if (do_save)
    figName = 'fig11B_normalized_alt_coefficients_by_time_since_last_visit_pvalues';
    saveas(gcf,[figFolder figName],'fig')
    saveas(gcf,[figFolder figName],'epsc')
    saveas(gcf,[figFolder figName],'svg')
    print([figFolder figName],'-dpng','-r600')
end

%% Look at effect of time with no distractor. Plot bout duration vs. time since last alt visit for bouts w/w.o. intervening distractor
difSolDayInds = find(boutTable.solutionNum ~= boutTable.alternative_solutionNum);
allXWithDistractor = cell(1,1);
allYWithDistractor = cell(1,1);
allXWithoutDistractor = cell(1,1);
allYWithoutDistractor = cell(1,1);
X = cell(1,1);
Y = cell(1,1);
for i=1:length(animalObjs)
    allXWithDistractor{i} = [];
    allYWithDistractor{i} = [];
    allXWithoutDistractor{i} = [];
    allYWithoutDistractor{i} = [];
    ratInds = find(strcmp(boutTable.rat,animalObjs(i).name));
    for j=1:length(animalObjs(i).difSolutionDays)
        dayInds = find(strcmp(boutTable.date,animalObjs(i).dates{animalObjs(i).difSolutionDays(j)}));
        curBouts = boutTable(intersect(ratInds,dayInds),:);
        [~,sortedInds] = sort(curBouts.onset);
        curBouts = curBouts(sortedInds,:);
        if (height(curBouts) < 2)
            continue;
        else
            for k=2:height(curBouts)
                pastBouts = curBouts(1:(k-1),:);
                altChannelBoutInds = find(pastBouts.channel ~= curBouts.channel(k));
                if (isempty(altChannelBoutInds))
                    continue;
                else
                    lastAltVisitBout = pastBouts(max(altChannelBoutInds),:);
                    timeSinceLastAltVisit = curBouts.onset(k) - lastAltVisitBout.offset;
                    wasDistractor = curBouts.channel(k) == pastBouts.channel(end);
                    if (wasDistractor)
                        allXWithDistractor{i} = [allXWithDistractor{i}; curBouts.palatability(k) lastAltVisitBout.palatability timeSinceLastAltVisit];
                        allYWithDistractor{i} = [allYWithDistractor{i}; curBouts.duration(k)];
                    else
                        allXWithoutDistractor{i} = [allXWithoutDistractor{i}; curBouts.palatability(k) lastAltVisitBout.palatability timeSinceLastAltVisit];
                        allYWithoutDistractor{i} = [allYWithoutDistractor{i}; curBouts.duration(k)];
                    end
                end
            end
        end
    end
end

bins = [60 6000];
normalizedCoeffsWithDistractor = nan(length(bins),3);
normalizedCoeffsWithoutDistractor = nan(length(bins),3);
altSolCoeffsWithDistractor = nan(length(animalObjs),length(bins));
altSolCoeffsWithoutDistractor = nan(length(animalObjs),length(bins));
allXWD = cell(1,2); allXWD{1} = []; allXWD{2} = [];
allXWOD = cell(1,2); allXWOD{1} = []; allXWOD{2} = [];
allYWD = cell(1,2); allYWD{1} = []; allYWD{2} = [];
allYWOD = cell(1,2); allYWOD{1} = []; allYWOD{2} = [];
for i=1:length(animalObjs)
    for j = 1:length(bins)
        if (j == 1)
            indsWD = find(allXWithDistractor{i}(:,3) < bins(j));
            indsWOD = find(allXWithoutDistractor{i}(:,3) < bins(j));
        else
            indsWD = find(allXWithDistractor{i}(:,3) < bins(j) & allXWithDistractor{i}(:,3) > bins(j-1));
            indsWOD = find(allXWithoutDistractor{i}(:,3) < bins(j) & allXWithoutDistractor{i}(:,3) > bins(j-1));
        end
        
        XWD{j} = [allXWithDistractor{i}(indsWD,1:2) ones(length(indsWD),1)];
        XWOD{j} = [allXWithoutDistractor{i}(indsWOD,1:2) ones(length(indsWOD),1)];
        YWD{j} = allYWithDistractor{i}(indsWD);
        YWOD{j} = allYWithoutDistractor{i}(indsWOD);
        allXWD{j} = [allXWD{j}; [allXWithDistractor{i}(indsWD,:) ones(length(indsWD),1)]];
        allXWOD{j} = [allXWOD{j}; [allXWithoutDistractor{i}(indsWOD,:) ones(length(indsWOD),1)]];
        allYWD{j} = [allYWD{j}; allYWithDistractor{i}(indsWD)];
        allYWOD{j} = [allYWOD{j}; allYWithoutDistractor{i}(indsWOD)];
        skipFlag = false;
        if (size(XWD{j},1) < 2)
            disp(['animal ' num2str(i) ' bin ' num2str(j) ' thrown out (N = ' num2str(size(X{j},1)) ') : WD'])
            skipFlag = true;
        end
        if (size(XWOD{j},1) < 2)
            disp(['animal ' num2str(i) ' bin ' num2str(j) ' thrown out (N = ' num2str(size(X{j},1)) ') : WOD'])
            skipFlag = true;
        end
        if (skipFlag)
            continue;
        else
            [b,bint,r,rint,stats] = regress(YWD{j},XWD{j});
            normalizedCoeffsWithDistractor(j,:) = b./mean(YWD{j});
            altSolCoeffsWithDistractor(i,j) = normalizedCoeffsWithDistractor(j,2);
            
            [b,bint,r,rint,stats] = regress(YWOD{j},XWOD{j});
            normalizedCoeffsWithoutDistractor(j,:) = b./mean(YWOD{j});
            altSolCoeffsWithoutDistractor(i,j) = normalizedCoeffsWithoutDistractor(j,2);
        end
    end
end

figure;
subplot(1,2,1)
notBoxPlot([altSolCoeffsWithoutDistractor(:,1) altSolCoeffsWithDistractor(:,1)])
set(gca,'xticklabel',{'No distractor', 'Distractor'})
ylabel(['\beta (\Delta t < ' num2str(bins(1)) 's)'])

subplot(1,2,2)
notBoxPlot([altSolCoeffsWithoutDistractor(:,2) altSolCoeffsWithDistractor(:,1)])
set(gca,'xticklabels',{'No distractor', 'Distractor'})
ylabel(['\beta (\Delta t < ' num2str(bins(2)) 's)'])

set(gcf,'Position',[10 10 1800 800])

figure;
histogram(allXWD{1}(:,3),'normalization','pdf','BinWidth',1)
hold on;
histogram(allXWOD{1}(:,3),'normalization','pdf','BinWidth',1)
legend({'Distractor','No Distractor'})
xlabel(['\Delta t'])
ylabel('Probability density')

deltaTWD = allXWD{1}(:,3);
deltaTWOD = allXWOD{1}(:,3);
[sortedDeltaTWD, sortIndsTWD] = sort(deltaTWD);
[sortedDeltaTWOD, sortIndsTWOD] = sort(deltaTWOD);

diffMatrix = abs(sortedDeltaTWD - sortedDeltaTWOD');
[M,uR,uC] = matchpairs(diffMatrix,.1);
origIndsWD = sortIndsTWD(M(:,1));
origIndsWOD = sortIndsTWOD(M(:,2));

figure;
histogram(allYWD{1}(origIndsWD),'BinWidth',2)
hold on;
histogram(allYWOD{1}(origIndsWOD),'BinWidth',2)

[bWD,bintWD,rWD,rintWD,statsWD] = regress(allYWD{1}(origIndsWD),allXWD{1}(origIndsWD,[1:2 4]));
[bWOD,bintWOD,rWOD,rintWOD,statsWOD] = regress(allYWOD{1}(origIndsWOD),allXWOD{1}(origIndsWOD,[1:2 4]));
regSTATSWD = regstats(allYWD{1}(origIndsWD),allXWD{1}(origIndsWD,[1 2]));
regSTATSWOD = regstats(allYWOD{1}(origIndsWOD),allXWOD{1}(origIndsWOD,[1 2]));


%{
allXNoDistractor = [];
allYNoDistractor = [];
allXDistractor = [];
allYDistractor = [];
for i=1:length(animalObjs)
    allXNoDistractor = [allXNoDistractor; allXWithoutDistractor{i}];
    allYNoDistractor = [allYNoDistractor; allYWithoutDistractor{i}];
    allXDistractor = [allXDistractor; allXWithDistractor{i}];
    allYDistractor = [allYDistractor; allYWithDistractor{i}];
end

figure;
scatter(allXNoDistractor(:,3),allYNoDistractor,20,'filled','MarkerEdgeColor','red','MarkerFaceColor','red')
hold on;
scatter(allXDistractor(:,3),allYDistractor,20,'filled','MarkerEdgeColor','blue','MarkerFaceColor','blue')
legend({'No distractor','Distractor'})
xlabel('Time since last visit to alternative (s)')
ylabel('Bout duration')
%}