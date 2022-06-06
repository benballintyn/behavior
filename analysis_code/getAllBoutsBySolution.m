function [allBoutsBySolution] = getAllBoutsBySolution(basedir,dataDates,animal)
allBoutsBySolution = cell(1,4);
for i=1:length(dataDates)
    dataDir = [basedir '/' dataDates{i} '/' animal];
    bouts = load([dataDir '/bouts.mat']); bouts = bouts.bouts;
    for j=1:length(bouts)
        if (contains(bouts{j}(1).solution,'H2O'))
            allBoutsBySolution{1} = [allBoutsBySolution{1} bouts{j}];
        elseif (contains(bouts{j}(1).solution,'A'))
            allBoutsBySolution{2} = [allBoutsBySolution{2} bouts{j}];
        elseif (contains(bouts{j}(1).solution,'B'))
            allBoutsBySolution{3} = [allBoutsBySolution{3} bouts{j}];
        elseif (contains(bouts{j}(1).solution,'C'))
            allBoutsBySolution{4} = [allBoutsBySolution{4} bouts{j}];
        else
            error(['Solution of bouts{' num2str(j) '} from date ' dataDates{i} ' cannot be recognized'])
        end
    end
end
end

