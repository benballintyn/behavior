function [T,transitionCount,totalTransitions,nswitches,nbouts] = getBoutTransitionMatrix(bouts,solnIDs)
nBottles = length(bouts);
T = zeros(nBottles);
transitionCount = zeros(nBottles);
totalTransitions = zeros(1,nBottles);
allBouts = [];
allBoutIDs = [];
foundIDs = [];
for i=1:length(bouts)
    if (~isempty(bouts{i}))
        foundIDs = [foundIDs solnConverter(bouts{i}(1).solution)];
    end
end
for i=1:length(bouts)
    if (~isempty(bouts{i}))
        id = solnConverter(bouts{i}(1).solution);
        if (length(unique(solnIDs)) > 1 && id == min(solnIDs))
            boutSolnIDs(i) = 1;
        elseif (length(unique(solnIDs)) > 1 && id == max(solnIDs))
            boutSolnIDs(i) = length(solnIDs);
        elseif (length(unique(solnIDs)) > 1)
            boutSolnIDs(i) = 2;
        else
            boutSolnIDs(i) = i;
        end
    else
        if (length(unique(solnIDs)) == 1)
            boutSolnIDs(i) = i;
        else
            id = setdiff(solnIDs,foundIDs);
            if (id == min(solnIDs))
                boutSolnIDs(i) = 1;
            elseif (id == max(solnIDs))
                boutSolnIDs(i) = length(solnIDs);
            else
                boutSolnIDs(i) = 2;
            end
        end
    end
end
for i=1:length(bouts)
    allBouts = [allBouts bouts{i}];
    allBoutIDs = [allBoutIDs ones(1,length(bouts{i}))*boutSolnIDs(i)];
end
[~,sortedInds] = sort([allBouts.onset]);
sortedIDs = allBoutIDs(sortedInds);
nswitches = 0;
for i=2:length(sortedIDs)
    if (sortedIDs(i) ~= sortedIDs(i-1))
        nswitches = nswitches + 1;
    end
    transitionCount(sortedIDs(i),sortedIDs(i-1)) = transitionCount(sortedIDs(i),sortedIDs(i-1)) + 1;
    totalTransitions(sortedIDs(i-1)) = totalTransitions(sortedIDs(i-1)) + 1;
end
for i=1:nBottles
    if (totalTransitions(i) ~= 0)
        T(:,i) = transitionCount(:,i)./totalTransitions(i);
    end
end
nbouts = length(allBouts);
end

function convSoln=solnConverter(soln)
    if (contains(soln,'H2O'))
        convSoln = 1;
    elseif (contains(soln,'.01M') || contains(soln,'A'))
        convSoln = 2;
    elseif (contains(soln,'.1M') || contains(soln,'B'))
        convSoln = 3;
    elseif (contains(soln,'1M') || contains(soln,'C'))
        convSoln = 4;
    elseif (soln == 1)
        convSoln = 'H2O';
    elseif (soln == 2)
        convSoln = '.01M';
    elseif (soln == 3)
        convSoln = '.1M';
    elseif (soln == 4)
        convSoln = '1M';
    else
        error('soln not recognized')
    end
end