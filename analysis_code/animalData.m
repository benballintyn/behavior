classdef animalData < matlab.mixin.Copyable
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name
        dates
        sex
        datadir
        amountsConsumed
        totalConsumed
        totalBottles
        totalLicks
        relativePalatabilitiesConsumed
        relativePalatabilitiesLicks
        solutions
        Palatabilities
        Palatabilities200
        meanAlternativePalatabilities
        meanAlternativePalatabilities200
        maxAlternativePalatabilities
        maxAlternativePalatabilities200
        licks
        nLicksByBoxSide
        bouts
        bouts200
        bout_split_time
        earlyBouts
        earlyBouts200
        lateBouts
        lateBouts200
        nonFirstBouts
        nonEarlyBouts
        boutsBySolution
        boutsBySolutionDifSolDays
        boutsGMcriteria
        linLicks
        linLicksSameSolDays
        linLicksDifSolDays
        linBouts
        linBoutsByDay
        linSameSolBouts
        linDifSolBouts
        ilis
        ilisBySolution
        timeToFirstLick
        sameSolutionDays
        difSolutionDays
    end
    
    methods
        function obj = animalData(name,dates,sex,datadir)
            obj.name = name;
            obj.dates = dates;
            obj.sex = sex;
            obj.datadir = datadir;
            obj.totalConsumed = zeros(1,4);
            obj.totalBottles = zeros(1,4);
            obj.totalLicks = zeros(1,4);
            obj.loadMetadata();
            obj.loadLicks();
            obj.getRelativePalatabilities();
            obj.loadBouts();
            obj.sortBoutsByDay();
            obj.getBoutsWithCriteria();
            obj.getILIs();
        end
        
        function loadLicks(obj)
            alllicks = [];
            linLicksSameSolDays = [];
            linLicksDifSolDays = [];
            obj.timeToFirstLick = zeros(1,length(obj.dates));
            for i=1:length(obj.dates)
                t1 = 3900;
                tmp = load([obj.datadir '/' obj.dates{i} '/' obj.name '/licks.mat']); tmp = tmp.licks;
                for j=1:length(tmp)
                    if (isfield(tmp{j},'certainty'))
                        tmp{j} = rmfield(tmp{j},'certainty');
                    end
                    alllicks = [alllicks tmp{j}];
                    if (ismember(i,obj.sameSolutionDays))
                        linLicksSameSolDays = [linLicksSameSolDays tmp{j}];
                    else
                        linLicksDifSolDays = [linLicksDifSolDays tmp{j}];
                    end
                    if (~isempty(tmp{j}))
                        if (tmp{j}(1).onset < t1)
                            t1 = tmp{j}(1).onset;
                        end
                        solnID = obj.solnConverter(tmp{j}(1).solution);
                        obj.totalLicks(solnID) = obj.totalLicks(solnID) + length(tmp{j});
                    end
                end
                obj.licks{i} = tmp;
                obj.timeToFirstLick(i) = t1;
            end
            obj.linLicks = alllicks;
            obj.linLicksSameSolDays = linLicksSameSolDays;
            obj.linLicksDifSolDays = linLicksDifSolDays;
        end
        
        function loadBouts(obj)
            obj.boutsBySolution = cell(1,4);
            obj.boutsBySolutionDifSolDays = cell(1,4);
            allSameSolBouts = [];
            allDifSolBouts = [];
            for i=1:length(obj.dates)
                tmp = load([obj.datadir '/' obj.dates{i} '/' obj.name '/bouts.mat']); tmp=tmp.bouts;
                for j=1:length(tmp)
                    if (isempty(tmp{j}))
                        continue;
                    else
                        singleLickBouts = ([tmp{j}.nlicks] == 1);
                        tmp{j} = tmp{j}(~singleLickBouts);
                        if (isempty(tmp{j}))
                            continue;
                        end
                        if (ismember(i,obj.sameSolutionDays))
                            allSameSolBouts = [allSameSolBouts tmp{j}];
                        else
                            allDifSolBouts = [allDifSolBouts tmp{j}];
                        end
                        solnID = obj.solnConverter(tmp{j}(1).solution);
                        box_side = tmp{j}(1).box_side;
                        if (strcmp(box_side,'left'))
                            obj.nLicksByBoxSide(i,1) = sum([tmp{j}.nlicks]);
                        elseif (strcmp(box_side,'right'))
                            obj.nLicksByBoxSide(i,2) = sum([tmp{j}.nlicks]);
                        else
                            error('Box side not recognized')
                        end
                        obj.Palatabilities{i}{j} = ones(1,length(tmp{j}))*obj.relativePalatabilitiesLicks(solnID);
                        alternatives = setdiff(obj.solutions(i,:),solnID);
                        if (isempty(alternatives))
                            alternatives = solnID;
                        end
                        meanAltPal = mean(obj.relativePalatabilitiesLicks(alternatives));
                        maxAltPal = max(obj.relativePalatabilitiesLicks(alternatives));
                        obj.meanAlternativePalatabilities{i}{j} = ones(1,length(tmp{j}))*meanAltPal;
                        obj.maxAlternativePalatabilities{i}{j} = ones(1,length(tmp{j}))*maxAltPal;
                        %obj.totalLicks(solnID) = obj.totalLicks(solnID) + sum([tmp{j}.nlicks]);
                        obj.boutsBySolution{solnID} = [obj.boutsBySolution{solnID} tmp{j}];
                        if (ismember(i,obj.difSolutionDays))
                            obj.boutsBySolutionDifSolDays{solnID} = [obj.boutsBySolutionDifSolDays{solnID} tmp{j}];
                        end
                    end
                end
                obj.bouts{i} = tmp;
            end
            obj.linSameSolBouts = allSameSolBouts;
            obj.linDifSolBouts = allDifSolBouts;
            obj.linBouts = [allSameSolBouts allDifSolBouts];
            obj.getNonFirstBouts();
            obj.getNonEarlyBouts();
            obj.getBouts200();
        end
        
        function getNonFirstBouts(obj)
            for i=1:length(obj.dates)
                firstBoutStart = inf;
                earliestBoutChan = nan;
                for j=1:length(obj.bouts{i})
                    if (~isempty(obj.bouts{i}{j}))
                        if (obj.bouts{i}{j}(1).onset < firstBoutStart)
                            earliestBoutChan = j;
                        end
                    end
                end
                tmp = obj.bouts{i};
                if (length(tmp{earliestBoutChan}) > 1)
                    tmp{earliestBoutChan} = tmp{earliestBoutChan}(2:end);
                else
                    tmp{earliestBoutChan} = [];
                end
                obj.nonFirstBouts{i} = tmp;
            end
        end
        
        function getNonEarlyBouts(obj)
            for i=1:length(obj.dates)
                for j=1:length(obj.bouts{i})
                    if (~isempty(obj.bouts{i}{j}))
                        onsets = [obj.bouts{i}{j}.onset];
                        nonEarlyBouts = find(onsets > 600);
                        if (~isempty(nonEarlyBouts))
                            obj.nonEarlyBouts{i}{j} = obj.bouts{i}{j}(nonEarlyBouts);
                        end
                    else
                        obj.nonEarlyBouts{i}{j} = [];
                    end
                end
            end
        end
        
        function getBouts200(obj)
            for i=1:length(obj.licks)
                for j=1:length(obj.licks{i})
                    curLickArr = obj.licks{i}{j};
                    if (~isempty(curLickArr))
                        box_side = curLickArr(1).box_side;
                        channel = curLickArr(1).channel;
                    end
                    curBouts200 = getBouts2(curLickArr,.2);
                    newCurBouts200 = [];
                    if (~isempty(curBouts200))
                        solnID = obj.solnConverter(curBouts200(1).solution);
                        alternatives = setdiff(obj.solutions(i,:),solnID);
                        if (isempty(alternatives))
                            alternatives = solnID;
                        end
                        for l=1:length(curBouts200)
                            if (curBouts200(l).nlicks > 1)
                                curBouts200(l).box_side = box_side;
                                curBouts200(l).channel = channel;
                                newCurBouts200 = [newCurBouts200 curBouts200(l)];
                            end
                        end
                    end
                    obj.bouts200{i}{j} = newCurBouts200;
                    obj.Palatabilities200{i}{j} = ones(1,length(newCurBouts200))*obj.relativePalatabilitiesLicks(solnID);
                    meanAltPal = mean(obj.relativePalatabilitiesLicks(alternatives));
                    maxAltPal = max(obj.relativePalatabilitiesLicks(alternatives));
                    obj.meanAlternativePalatabilities200{i}{j} = ones(1,length(newCurBouts200))*meanAltPal;
                    obj.maxAlternativePalatabilities200{i}{j} = ones(1,length(newCurBouts200))*maxAltPal;
                end
            end
        end
        
        function sortBoutsByDay(obj)
            for i=1:length(obj.dates)
                obj.linBoutsByDay{i} = [];
                for j=1:length(obj.bouts{i})
                    obj.linBoutsByDay{i} = [obj.linBoutsByDay{i} obj.bouts{i}{j}];
                end
                [~,inds] = sort([obj.linBoutsByDay{i}.onset]);
                obj.linBoutsByDay{i} = obj.linBoutsByDay{i}(inds);
            end
        end
        
        function getILIs(obj)
            obj.ilis = [];
            obj.ilisBySolution = cell(1,4);
            for i=1:4
                for j=1:length(obj.boutsBySolution{i})
                    onsets = [obj.boutsBySolution{i}(j).licks.onset];
                    offsets = [obj.boutsBySolution{i}(j).licks.offset];
                    curILIs = onsets(2:end)-offsets(1:end-1);
                    curmax = max(curILIs);
                    if (curmax > 20)
                        disp([obj.name ' ' num2str(i) ' ' num2str(j)])
                    end
                    obj.ilisBySolution{i} = [obj.ilisBySolution{i} curILIs];
                    obj.ilis = [obj.ilis curILIs];
                end
            end
        end
        
        function getBoutsWithCriteria(obj)
            x = 0:1:3900;
            F = ksdensity([obj.linLicks.onset],x,'Bandwidth',200);
            Fcum = cumsum(F);
            diff2 = diff(Fcum,2);
            [~,min_diff2_ind] = min(diff2);
            obj.bout_split_time = x(min_diff2_ind);
            for i=1:length(obj.dates)
                % 2s bouts
                for j=1:length(obj.bouts{i})
                    if (~isempty(obj.bouts{i}{j}))
                        onsets = [obj.bouts{i}{j}.onset];
                        earlyBoutInds = find(onsets <= obj.bout_split_time);
                        lateBoutInds = find(onsets > obj.bout_split_time);
                        obj.earlyBouts{i}{j} = obj.bouts{i}{j}(earlyBoutInds);
                        obj.lateBouts{i}{j} = obj.bouts{i}{j}(lateBoutInds);
                    else
                        obj.earlyBouts{i}{j} = [];
                        obj.lateBouts{i}{j} = [];
                    end
                end
                % 200ms bouts
                for j=1:length(obj.bouts200{i})
                    if (~isempty(obj.bouts200{i}{j}))
                        onsets = [obj.bouts200{i}{j}.onset];
                        earlyBoutInds = find(onsets <= obj.bout_split_time);
                        lateBoutInds = find(onsets > obj.bout_split_time);
                        obj.earlyBouts200{i}{j} = obj.bouts200{i}{j}(earlyBoutInds);
                        obj.lateBouts200{i}{j} = obj.bouts200{i}{j}(lateBoutInds);
                    else
                        obj.earlyBouts200{i}{j} = [];
                        obj.lateBouts200{i}{j} = [];
                    end
                end
            end
        end
        
        function loadMetadata(obj)
            for i=1:length(obj.dates)
                metadata = load([obj.datadir '/' obj.dates{i} '/' obj.name '/metadata.mat']); metadata=metadata.metadata;
                count=0;
                if (~isempty(metadata.leftSolution))
                    count=count+1;
                    solnID = obj.solnConverter(metadata.leftSolution);
                    obj.solutions(i,count) = solnID;
                    obj.amountsConsumed(i,count) = metadata.leftConsumed;
                    if (~isnan(metadata.leftConsumed))
                        obj.totalConsumed(solnID) = obj.totalConsumed(solnID) + metadata.leftConsumed;
                    end
                    obj.totalBottles(solnID) = obj.totalBottles(solnID) + 1;
                end
                if (~isempty(metadata.middleSolution))
                    count=count+1;
                    solnID = obj.solnConverter(metadata.middleSolution);
                    obj.solutions(i,count) = solnID;
                    obj.amountsConsumed(i,count) = metadata.middleConsumed;
                    if (~isnan(metadata.middleConsumed))
                        obj.totalConsumed(solnID) = obj.totalConsumed(solnID) + metadata.middleConsumed;
                    end
                    obj.totalBottles(solnID) = obj.totalBottles(solnID) + 1;
                end
                if (~isempty(metadata.rightSolution))
                    count=count+1;
                    solnID = obj.solnConverter(metadata.rightSolution);
                    obj.solutions(i,count) = solnID;
                    obj.amountsConsumed(i,count) = metadata.rightConsumed;
                    if (~isnan(metadata.rightConsumed))
                        obj.totalConsumed(solnID) = obj.totalConsumed(solnID) + metadata.rightConsumed;
                    end
                    obj.totalBottles(solnID) = obj.totalBottles(solnID) + 1;
                end
            end
        end
        
        function getRelativePalatabilities(obj)
            %{
            totalConsumed = obj.totalConsumed;
            totalLicks = obj.totalLicks;
            scales = obj.totalBottles./obj.totalBottles(1);
            scaledConsumed = totalConsumed./scales;
            scaledLicks = totalLicks./scales;
            obj.relativePalatabilitiesConsumed = scaledConsumed./scaledConsumed(1);
            obj.relativePalatabilitiesLicks = scaledLicks./scaledLicks(1);
            %}
            dayCounts = zeros(1,4);
            consumed_sum = zeros(1,4);
            licks_sum = zeros(1,4);
            obj.sameSolutionDays = [];
            for i=1:4
                for j=1:size(obj.solutions,1)
                    if (all(obj.solutions(j,:) == [i i]))
                        dayCounts(i) = dayCounts(i) + 1;
                        consumed_sum(i) = consumed_sum(i) + sum(obj.amountsConsumed(j,:));
                        for k=1:length(obj.licks{j})
                            licks_sum(i) = licks_sum(i) + length(obj.licks{j}{k});
                        end
                        obj.sameSolutionDays = [obj.sameSolutionDays j];
                    end
                end
            end
            scales = dayCounts./dayCounts(1);
            scaledConsumed = consumed_sum./scales;
            scaledLicks = licks_sum./scales;
            obj.relativePalatabilitiesConsumed = scaledConsumed./scaledConsumed(1);
            obj.relativePalatabilitiesLicks = scaledLicks./scaledLicks(1);
            obj.difSolutionDays = setdiff(1:length(obj.dates),obj.sameSolutionDays);
        end
        
        function plotLicks(obj,day)
            plotLicks(obj.datadir,day,obj.name)
        end
    end

    methods(Static)
        function convSoln=solnConverter(soln)
            if (ischar(soln))
                if (contains(soln,'H2O'))
                    convSoln = 1;
                elseif (contains(soln,'.01M') || contains(soln,'A'))
                    convSoln = 2;
                elseif (contains(soln,'.1M') || contains(soln,'B'))
                    convSoln = 3;
                elseif (contains(soln,'1M') || contains(soln,'C'))
                    convSoln = 4;
                else
                    error('Solution is a character array but not recognized')
                end
            elseif (isnumeric(soln))
                if (soln == 1)
                    convSoln = 'H2O';
                elseif (soln == 2)
                    convSoln = '.01M';
                elseif (soln == 3)
                    convSoln = '.1M';
                elseif (soln == 4)
                    convSoln = '1M';
                else
                    error('Solution is numeric but not recognized')
                end
            else
                error('Type of soln not recognized (not char or numeric)')
            end
        end
    end
end

