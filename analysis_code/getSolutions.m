function [solutions,nSolns] = getSolutions(basedir,animal,dates)
firstDir = [basedir '/' dates{1} '/' animal];
if (~exist([firstDir '/manual_flag.mat'],'file'))
    error([firstDir ' was not manually analyzed'])
else
    bouts1 = load([firstDir '/bouts.mat']); bouts1=bouts1.bouts;
    nSolns = length(bouts1);
    solutions = cell(length(dates),nSolns);
    for i=1:length(dates)
        curdir = [basedir '/' dates{i} '/' animal];
        metadata = load([curdir '/metadata.mat']); metadata=metadata.metadata;
        if (isempty(metadata.leftSolution))
            solutions{i,1} = 'N/A';
        else
            solutions{i,1} = metadata.leftSolution;
        end
        if (isempty(metadata.middleSolution))
            solutions{i,2} = 'N/A';
        else
            solutions{i,2} = metadata.middleSolution;
        end
        if (isempty(metadata.rightSolution))
            solutions{i,3} = 'N/A';
        else
            solutions{i,3} = metadata.rightSolution;
        end
    end   
end
end

