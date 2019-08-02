function [allBouts,meanAlts,maxAlts,altDifs] = combineBouts(basedir,dates,animal)
allBouts = [];
meanAlts = [];
maxAlts = [];
altDifs = [];
for i=1:length(dates)
    curdir = [basedir '/' dates{i} '/' animal];
    bouts = load([curdir '/bouts.mat']);
    bouts = bouts.bouts;
    for j=1:length(bouts)
        allRewards(j) = getPalatability(getLogNaClConcentration(bouts{j}(1).solution));
    end
    for j=1:length(bouts)
        curReward = getPalatability(getLogNaClConcentration(bouts{j}(1).solution));
        allBouts = [allBouts bouts{j}];
        otherRewards = setdiff(allRewards,curReward);
        if (isempty(otherRewards))
            otherRewards = curReward;
        end
        meanAlts = [meanAlts ones(1,length(bouts{j}))*mean(otherRewards)];
        maxAlts = [maxAlts ones(1,length(bouts{j}))*max(otherRewards)];
        altDifs = [altDifs (curReward - ones(1,length(bouts{j}))*mean(otherRewards))];
    end
end
end

