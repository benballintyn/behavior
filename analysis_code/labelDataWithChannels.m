function [data] = labelDataWithChannels(basedir,dataDate,animal)
datadir = [basedir '/' dataDate '/' animal];
data = load([datadir '/data.mat']); data=data.data;
bouts = load([datadir '/bouts.mat']); bouts=bouts.bouts;
licks = load([datadir '/licks.mat']); licks=licks.licks;
metadata = load([datadir '/metadata.mat']); metadata=metadata.metadata;
for i=1:length(data)
    solnConcentration{i} = getLogNaClConcentration(data(i).solution);
end
if (length(data) == 2)
    if (solnConcentration{1} == solnConcentration{2} || (isnan(solnConcentration{1}) && isnan(solnConcentration{2})))
        idx1 = strfind(data(1).datafile,'_');
        if (length(idx1) ~= 3)
            error('something went wrong processing data(1).datafile')
        else
            datafileNum(1) = data(1).datafile(idx1(3)-1);
        end
        idx2 = strfind(data(2).datafile,'_');
        if (length(idx2) ~= 3)
            error('something went wrong processing data(2).datafile')
        else
            datafileNum(2) = data(2).datafile(idx2(3)-1);
        end
        lchan = metadata.leftChannel;
        rchan = metadata.rightChannel;
        if (lchan < rchan)
            [~,leftChannelInd] = min(datafileNum);
            data(leftChannelInd).box_side = 'left';
            data(leftChannelInd).channel = metadata.leftChannel;
            data(leftChannelInd).amount_consumed = metadata.leftConsumed;
            for i=1:length(bouts{leftChannelInd})
                bouts{leftChannelInd}(i).box_side = 'left';
                bouts{leftChannelInd}(i).channel = metadata.leftChannel;
                if (~strcmp(data(leftChannelInd).solution,bouts{leftChannelInd}(i).solution))
                    disp(['Changing bouts{' num2str(leftChannelInd) '}(' num2str(i) ').solution to ' data(leftChannelInd).solution])
                    bouts{leftChannelInd}(i).solution = data(leftChannelInd).solution;
                end
            end
            for i=1:length(licks{leftChannelInd})
                licks{leftChannelInd}(i).box_side = 'left';
                licks{leftChannelInd}(i).channel = metadata.leftChannel;
                if (~strcmp(data(leftChannelInd).solution,licks{leftChannelInd}(i).solution))
                    disp(['Changing licks{' num2str(leftChannelInd) '}(' num2str(i) ').solution to ' data(leftChannelInd).solution])
                    licks{leftChannelInd}(i).solution = data(leftChannelInd).solution;
                end
            end
            metadata.leftChannelDatafile = data(leftChannelInd).datafile;
            [~,rightChannelInd] = max(datafileNum);
            data(rightChannelInd).box_side = 'right';
            data(rightChannelInd).channel = metadata.rightChannel;
            data(rightChannelInd).amount_consumed = metadata.rightConsumed;
            for i=1:length(bouts{rightChannelInd})
                bouts{rightChannelInd}(i).box_side = 'right';
                bouts{rightChannelInd}(i).channel = metadata.rightChannel;
                if (~strcmp(data(rightChannelInd).solution,bouts{rightChannelInd}(i).solution))
                    disp(['Changing bouts{' num2str(rightChannelInd) '}(' num2str(i) ').solution to ' data(rightChannelInd).solution])
                    bouts{rightChannelInd}(i).solution = data(rightChannelInd).solution;
                end
            end
            for i=1:length(licks{rightChannelInd})
                licks{rightChannelInd}(i).box_side = 'right';
                licks{rightChannelInd}(i).channel = metadata.rightChannel;
                if (~strcmp(data(rightChannelInd).solution,licks{rightChannelInd}(i).solution))
                    disp(['Changing licks{' num2str(rightChannelInd) '}(' num2str(i) ').solution to ' data(rightChannelInd).solution])
                    licks{rightChannelInd}(i).solution = data(rightChannelInd).solution;
                end
            end
            metadata.rightChannelDatafile = data(rightChannelInd).datafile;
        elseif (rchan < lchan)
            [~,rightChannelInd] = min(datafileNum);
            data(rightChannelInd).box_side = 'right';
            data(rightChannelInd).channel = metadata.rightChannel;
            data(rightChannelInd).amount_consumed = metadata.rightConsumed;
            for i=1:length(bouts{rightChannelInd})
                bouts{rightChannelInd}(i).box_side = 'right';
                bouts{rightChannelInd}(i).channel = metadata.rightChannel;
                if (~strcmp(data(rightChannelInd).solution,bouts{rightChannelInd}(i).solution))
                    disp(['Changing bouts{' num2str(rightChannelInd) '}(' num2str(i) ').solution to ' data(rightChannelInd).solution])
                    bouts{rightChannelInd}(i).solution = data(rightChannelInd).solution;
                end
            end
            for i=1:length(licks{rightChannelInd})
                licks{rightChannelInd}(i).box_side = 'right';
                licks{rightChannelInd}(i).channel = metadata.rightChannel;
                if (~strcmp(data(rightChannelInd).solution,licks{rightChannelInd}(i).solution))
                    disp(['Changing licks{' num2str(rightChannelInd) '}(' num2str(i) ').solution to ' data(rightChannelInd).solution])
                    licks{rightChannelInd}(i).solution = data(rightChannelInd).solution;
                end
            end
            metadata.rightChannelDataFile = data(rightChannelInd).datafile;
            [~,leftChannelInd] = max(datafileNum);
            data(leftChannelInd).box_side = 'left';
            data(leftChannelInd).channel = metadata.leftChannel;
            data(leftChannelInd).amount_consumed = metadata.leftConsumed;
            for i=1:length(bouts{leftChannelInd})
                bouts{leftChannelInd}(i).box_side = 'left';
                bouts{leftChannelInd}(i).channel = metadata.leftChannel;
                if (~strcmp(data(leftChannelInd).solution,bouts{leftChannelInd}(i).solution))
                    disp(['Changing bouts{' num2str(leftChannelInd) '}(' num2str(i) ').solution to ' data(leftChannelInd).solution])
                    bouts{leftChannelInd}(i).solution = data(leftChannelInd).solution;
                end
            end
            for i=1:length(licks{leftChannelInd})
                licks{leftChannelInd}(i).box_side = 'left';
                licks{leftChannelInd}(i).channel = metadata.leftChannel;
                if (~strcmp(data(leftChannelInd).solution,licks{leftChannelInd}(i).solution))
                    disp(['Changing licks{' num2str(leftChannelInd) '}(' num2str(i) ').solution to ' data(leftChannelInd).solution])
                    licks{leftChannelInd}(i).solution = data(leftChannelInd).solution;
                end
            end
            metadata.leftChannelDataFile = data(leftChannelInd).datafile;
        else
            error('right and left channel numbers are the same. this is bad')
        end
        save([datadir '/data.mat'],'data','-mat')
        save([datadir '/licks.mat'],'licks','-mat')
        save([datadir '/bouts.mat'],'bouts','-mat')
        save([datadir '/metadata.mat'],'metadata','-mat')
    else
        concentration{1} = getConcentrationFromSolnName(data(1).solution);
        concentration{2} = getConcentrationFromSolnName(data(2).solution);
        if (strcmp(concentration{1},metadata.leftSolution))
            data(1).box_side = 'left';
            data(1).channel = metadata.leftChannel;
            data(1).amount_consumed = metadata.leftConsumed;
            for i=1:length(bouts{1})
                bouts{1}(i).box_side = 'left';
                bouts{1}(i).channel = metadata.leftChannel;
                if (~strcmp(data(1).solution,bouts{1}(i).solution))
                    disp(['Changing bouts{1}(' num2str(i) ').solution to ' data(1).solution])
                    bouts{1}(i).solution = data(1).solution;
                end
            end
            for i=1:length(licks{1})
                licks{1}(i).box_side = 'left';
                licks{1}(i).channel = metadata.leftChannel;
                if (~strcmp(data(1).solution,licks{1}(i).solution))
                    disp(['Changing licks{1}(' num2str(i) ').solution to ' data(1).solution])
                    licks{1}(i).solution = data(1).solution;
                end
            end
            metadata.leftChannelDataFile = data(1).datafile;
            data(2).box_side = 'right';
            data(2).channel = metadata.rightChannel;
            data(2).amount_consumed = metadata.rightConsumed;
            for i=1:length(bouts{2})
                bouts{2}(i).box_side = 'right';
                bouts{2}(i).channel = metadata.rightChannel;
                if (~strcmp(data(2).solution,bouts{2}(i).solution))
                    disp(['Changing bouts{2}(' num2str(i) ').solution to ' data(2).solution])
                    bouts{2}(i).solution = data(2).solution;
                end
            end
            for i=1:length(licks{2})
                licks{2}(i).box_side = 'right';
                licks{2}(i).channel = metadata.rightChannel;
                if (~strcmp(data(2).solution,licks{2}(i).solution))
                    disp(['Changing licks{2}(' num2str(i) ').solution to ' data(2).solution])
                    licks{2}(i).solution = data(2).solution;
                end
            end
            metadata.rightChannelDataFile = data(2).datafile;
        elseif (strcmp(concentration{1},metadata.rightSolution))
            data(1).box_side = 'right';
            data(1).channel = metadata.rightChannel;
            data(1).amount_consumed = metadata.rightConsumed;
            for i=1:length(bouts{1})
                bouts{1}(i).box_side = 'right';
                bouts{1}(i).channel = metadata.rightChannel;
                if (~strcmp(data(1).solution,bouts{1}(i).solution))
                    disp(['Changing bouts{1}(' num2str(i) ').solution to ' data(1).solution])
                    bouts{1}(i).solution = data(1).solution;
                end
            end
            for i=1:length(licks{1})
                licks{1}(i).box_side = 'right';
                licks{1}(i).channel = metadata.rightChannel;
                if (~strcmp(data(1).solution,licks{1}(i).solution))
                    disp(['Changing licks{1}(' num2str(i) ').solution to ' data(1).solution])
                    licks{1}(i).solution = data(1).solution;
                end
            end
            metadata.rightChannelDataFile = data(1).datafile;
            data(2).box_side = 'left';
            data(2).channel = metadata.leftChannel;
            data(2).amount_consumed = metadata.leftConsumed;
            for i=1:length(bouts{2})
                bouts{2}(i).box_side = 'left';
                bouts{2}(i).channel = metadata.leftChannel;
                if (~strcmp(data(2).solution,bouts{2}(i).solution))
                    disp(['Changing bouts{2}(' num2str(i) ').solution to ' data(2).solution])
                    bouts{2}(i).solution = data(2).solution;
                end
            end
            for i=1:length(licks{2})
                licks{2}(i).box_side = 'left';
                licks{2}(i).channel = metadata.leftChannel;
                if (~strcmp(data(2).solution,licks{2}(i).solution))
                    disp(['Changing licks{2}(' num2str(i) ').solution to ' data(2).solution])
                    licks{2}(i).solution = data(2).solution;
                end
            end
            metadata.leftChannelDataFile = data(2).datafile;
        else
            error('data(1).solution could not be matched to either solution in metadata')
        end
        save([datadir '/data.mat'],'data','-mat')
        save([datadir '/licks.mat'],'licks','-mat')
        save([datadir '/bouts.mat'],'bouts','-mat')
        save([datadir '/metadata.mat'],'metadata','-mat')
    end
end
disp(['done labeling data for ' animal ' ' dataDate])
end

