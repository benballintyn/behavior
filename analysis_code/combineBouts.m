function [allBouts,meanAlts,maxAlts,altDifs,totalSampleDur,totalSampleAltDifs] = combineBouts(basedir,dates,animal,varargin)
allBouts = [];
meanAlts = [];
maxAlts = [];
altDifs = [];
totalSampleDur = [];
totalSampleAltDifs = [];
for i=1:length(dates)
    curdir = [basedir '/' dates{i} '/' animal];
    bouts = load([curdir '/bouts.mat']);
    bouts = bouts.bouts;
    if (contains(bouts{1}(1).solution,'H2O'))
        continue;
    end
    if (nargin > 3)
        customPalatability = varargin{1};
        if (length(customPalatability)~=4)
            error('custom palatability input is wrong size')
        else
            for j=1:length(bouts)
                allRewards(j) = lookupCustomPalatability(bouts{j}(1).solution,customPalatability);
            end
        end
    else
        for j=1:length(bouts)
            allRewards(j) = getPalatability(getLogNaClConcentration(bouts{j}(1).solution));
        end
    end
    for j=1:length(bouts)
        if (nargin > 3)
            curReward = lookupCustomPalatability(bouts{j}(1).solution,customPalatability);
        else
            curReward = getPalatability(getLogNaClConcentration(bouts{j}(1).solution));
        end
        allBouts = [allBouts bouts{j}];
        otherRewards = setdiff(allRewards,curReward);
        if (isempty(otherRewards))
            otherRewards = curReward;
        end
        totalSampleDur = [totalSampleDur sum([bouts{j}.duration])];
        totalSampleAltDifs = [totalSampleAltDifs (curReward - mean(otherRewards))];
        meanAlts = [meanAlts ones(1,length(bouts{j}))*mean(otherRewards)];
        maxAlts = [maxAlts ones(1,length(bouts{j}))*max(otherRewards)];
        altDifs = [altDifs (curReward - ones(1,length(bouts{j}))*mean(otherRewards))];
    end
end

function pal=lookupCustomPalatability(solution,palatabilities)
    if (contains(solution,'H2O'))
        pal = palatabilities(1);
    elseif (contains(solution,'A'))
        pal = palatabilities(2);
    elseif (contains(solution,'B'))
        pal = palatabilities(3);
    elseif (contains(solution,'C'))
        pal = palatabilities(4);
    else
        error('solution not recognized')
    end
end
end

