function relabelBoutAndLickSolutions(basedir,date,animal)
curdir = [basedir '/' date '/' animal];
data = load([curdir '/data.mat']); data=data.data;
bouts = load([curdir '/bouts.mat']); bouts=bouts.bouts;
licks = load([curdir '/licks.mat']); licks=licks.licks;

for i=1:length(bouts)
    if (~isempty(bouts{i}))
        for j=1:length(bouts{i})
            if (~strcmp(data(i).solution,bouts{i}(j).solution))
                disp(['Changing bout{' num2str(i) '}(' num2str(j) ').solution to ' getSolutionName(data(i).solution)])
                bouts{i}(j).solution = getSolutionName(data(i).solution);
            end
        end
    end
end

for i=1:length(licks)
    if (~isempty(licks{i}))
        for j=1:length(licks{i})
            if (~strcmp(data(i).solution,licks{i}(j).solution))
                disp(['Changing lick{' num2str(i) '}(' num2str(j) ').solution to ' getSolutionName(data(i).solution)])
                licks{i}(j).solution = getSolutionName(data(i).solution);
            end
        end
    end
end
save([curdir '/licks.mat'],'licks','-mat')
save([curdir '/bouts.mat'],'bouts','-mat')
end

function name = getSolutionName(solnString)
if (contains(solnString,'H2O'))
    name = 'H2O';
elseif (contains(solnString,'A'))
    name = 'A';
elseif (contains(solnString,'B'))
    name = 'B';
elseif (contains(solnString,'C'))
    name = 'C';
else
    error('solnString not recognized')
end
end

