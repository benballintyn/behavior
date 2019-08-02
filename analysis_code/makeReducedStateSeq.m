function [reducedStateVec] = makeReducedStateSeq(basedir,date,animal)
expStartOffset = 300;
curdir = [basedir '/' date '/' animal];
metadata = load([curdir '/metadata.mat']); metadata=metadata.metadata;
lchan = metadata.leftChannel;
mchan = metadata.middleChannel;
rchan = metadata.rightChannel;
lsoln = metadata.leftSolution;
msoln = metadata.middleSolution;
rsoln = metadata.rightSolution;
if (isempty(mchan))
    if (strcmp(lsoln,rsoln))
        [~,order] = sort([lchan rchan]);
    else
        [~,order] = sort({lsoln,rsoln});
    end
else % not robust to identical solutions in 3-bottle task
    [~,order] = sort({lsoln,msoln,rsoln});
    %[~,order] = sort([lchan mchan rchan]);
end
bouts = load([curdir '/bouts.mat']); bouts=bouts.bouts;
reducedStateVec = ones(1,36000);
for i=1:length(order)
    curBoutInd = find(order == i);
    for j=1:length(bouts{curBoutInd})
        t1 = round((bouts{curBoutInd}(j).onset-expStartOffset)/.1);
        t2 = round((bouts{curBoutInd}(j).offset-expStartOffset)/.1);
        if ((t2 - t1) < .1)
            reducedStateVec(t1) = i+1;
        else
            reducedStateVec(t1:t2) = i+1;
        end
    end
end
end

