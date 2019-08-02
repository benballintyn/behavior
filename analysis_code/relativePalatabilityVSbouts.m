function [allBouts,palDifs] = relativePalatabilityVSbouts(basedir,dates,animal)
allBouts = [];
palDifs = [];
for i=1:length(dates)
    curdir = [basedir '/' dates{i} '/' animal];
    bouts = load([curdir '/bouts.mat']); bouts = bouts.bouts;
    for j=1:length(bouts)
        licksBySoln(j) = sum([bouts{j}.nlicks]);
    end
    relPal = licksBySoln./sum(licksBySoln);
    for j=1:length(bouts)
        onesVec = ones(1,length(bouts));
        onesVec(j) = 0;
        curRelPal = relPal(j) - mean(relPal(logical(onesVec)));
        allBouts = [allBouts bouts{j}];
        palDifs = [palDifs ones(1,length(bouts{j}))*curRelPal];
    end
end
end

