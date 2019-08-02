function [stateVec] = makeStateVec(bouts,nT,dt,endBaseline,shuffleBouts,shuffleTime)
stateVec = ones(1,nT);
allBouts = [];
for i=1:length(bouts)
    for j=1:length(bouts{i})
        allBouts = [allBouts bouts{i}(j)];
    end
end
if (strcmp(shuffleBouts,'yes'))
    randOrder = randperm(length(allBouts));
    newBouts = allBouts;
    for i=1:length(allBouts)
        newBouts(i).solution = allBouts(randOrder(i)).solution;
    end
else
    newBouts = allBouts;
end

for i=1:length(newBouts)
    onsetInd = (round(newBouts(i).onset - endBaseline)*(1/dt));
    offsetInd = (round(newBouts(i).offset - endBaseline)*(1/dt));
    if (contains(newBouts(i).solution,'A'))
        stateVec(onsetInd:offsetInd) = 2;
    elseif (contains(newBouts(i).solution,'B'))
        stateVec(onsetInd:offsetInd) = 3;
    elseif (contains(newBouts(i).solution,'C'))
        stateVec(onsetInd:offsetInd) = 4;
    else
        error('bout solution is not valid')
    end
end
if (strcmp(shuffleTime,'yes'))
    stateVec = stateVec(randperm(length(stateVec)));
end

end

