function [allLicks] = getAllLicks(basedir,dates,animals)
allLicks = [];
for i=1:length(dates)
    for j=1:length(animals)
        dataDir = [basedir '/' dates{i} '/' animals{j}];
        licks = load([dataDir '/licks.mat']); licks=licks.licks;
        for k=1:length(licks)
            if (isfield(licks{k},'certainty'))
                licks{k} = rmfield(licks{k},'certainty');
            end
            allLicks = [allLicks licks{k}];
        end
    end
end
end

