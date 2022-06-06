function [naclBouts,h2oBouts,naclBoutDurations] = analyzeSolutionPairs(basedir,dates,animal)
naclBouts = cell(3,3);
h2oBouts = [];
for i=1:length(dates)
    dataDir = [basedir '/' dates{i} '/' animal];
    bouts = load([dataDir '/bouts.mat']); bouts=bouts.bouts;
    if (contains(bouts{1}(1).solution,'A'))
        ind1=1;
    elseif (contains(bouts{1}(1).solution,'B'))
        ind1=2;
    elseif (contains(bouts{1}(1).solution,'C'))
        ind1=3;
    elseif (contains(bouts{1}(1).solution,'H2O'))
        ind1=4;
    else
        error(['solution from bouts{' num2str(1) '} on day ' dates{i} ' cannot be recognized'])
    end
    if (contains(bouts{2}(1).solution,'A'))
        ind2=1;
    elseif (contains(bouts{2}(1).solution,'B'))
        ind2=2;
    elseif (contains(bouts{2}(1).solution,'C'))
        ind2=3;
    elseif (contains(bouts{2}(1).solution,'H2O'))
        ind2=4;
    else
        error(['solution from bouts{' num2str(2) '} on day ' dates{i} ' cannot be recognized'])
    end
    if (ind1 < 4 && ind2 < 4)
        naclBouts{ind1,ind2} = [naclBouts{ind1,ind2} bouts{1}];
        naclBouts{ind2,ind1} = [naclBouts{ind2,ind1} bouts{2}];
    elseif (ind1 == 4 && ind2 == 4)
        h2oBouts = [h2oBouts bouts{1} bouts{2}];
    else
        error('indices cannot be correct')
    end
end
for i=1:3
    for j=1:3
        naclBoutDurations{i,j} = [naclBouts{i,j}.duration];
    end
end
end

