function [meanCorr,meanPval,maxCorr,maxPval,values] = getBoutLengthCorrelations(allBouts)
boutLengths = zeros(length(allBouts),1);
meanValues = zeros(length(allBouts),1);
maxValues = zeros(length(allBouts),1);
solutions = [];
solnValues = [];
for i=1:length(allBouts)
    if (ismember(allBouts(i).solution,solutions))
        continue;
    else
        solutions = [solutions allBouts(i).solution];
        solnValues = [solnValues getPalatability(getLogNaClConcentration(allBouts(i).solution))];
    end
end
for i=1:length(allBouts)
    boutLengths(i) = allBouts(i).duration;
    meanValues(i) = mean(setdiff(solnValues,getPalatability(getLogNaClConcentration(allBouts(i).solution))));
    maxValues(i) = max(setdiff(solnValues,getPalatability(getLogNaClConcentration(allBouts(i).solution))));
end
[meanCorr,meanPval] = corr(boutLengths,meanValues);
[maxCorr, maxPval] = corr(boutLengths,maxValues);
values(:,1) = boutLengths;
values(:,2) = meanValues;
values(:,3) = maxValues;
end

