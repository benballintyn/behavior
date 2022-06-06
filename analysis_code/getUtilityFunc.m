function [func] = getUtilityFunc(driveDir,animal,dates)
allLickOnsets = [];
for i=1:length(dates)
    curdir = [driveDir '/' dates{i} '/' animal];
    licks = load([curdir '/licks.mat']); licks=licks.licks;
    for j=1:length(licks)
        if (isempty(licks{j}))
            continue;
        else
            allLickOnsets = [allLickOnsets licks{j}.onset];
        end
    end
end

parmhat = expfit(allLickOnsets);
func = @(t) exp(-t/parmhat);
end

