function [goodlicks] = findLicksWithLSTM(data,manualCorrection,savedir)
if (length(data) == 2)
    lickNet = load('analysis_code/bilstmLickNet2Bottle.mat'); lickNet=lickNet.bilstmLickNet;
    licks = cell(1,2);
    validLicks = cell(1,2);
    lickCount = zeros(1,2);
elseif (length(data) == 3)
    lickNet = load('analysis_code/lickNet3Bottle.mat'); lickNet=lickNet.lickNet3Bottle;
    licks = cell(1,3);
    validLicks = cell(1,3);
    lickCount = zeros(1,3);
else
    error('Length of data is not 2 or 3. What did you do wrong?')
end
X = zeros(length(data),length(data(1).raw_voltage));
for i=1:length(data)
    X(i,:) = data(i).raw_voltage';
end
subsetLength = 100000;
nbatches = ceil(length(data(1).raw_voltage)/subsetLength);
lickProbs = [];
for i=1:nbatches
    disp(['Prediction batch ' num2str(i) ' / ' num2str(nbatches)])
    start_ind = (i-1)*subsetLength + 1;
    if (i == nbatches)
        end_ind = length(data(1).raw_voltage);
    else
        end_ind = i*subsetLength;
    end
    [lickNet,lickprobs] = predictAndUpdateState(lickNet,X(:,start_ind:end_ind),'ExecutionEnvironment','gpu');
    lickProbs = [lickProbs lickprobs];
end
[~,maxVec] = max(lickProbs,[],1);
inLick = 0;
curLick = struct();
curSolnInd = [];
for t=1:length(maxVec)
    if (maxVec(t) == 1)
        if (inLick)
            curLick.offset_ind = t-1;
            curLick.offset = data(curSolnInd).tvec(t-1);
            curLick.offsetVal = data(curSolnInd).raw_voltage(t-1);
            curLick.duration = curLick.offset - curLick.onset;
            curLick.raw_voltage = data(curSolnInd).raw_voltage(curLick.onset_ind:curLick.offset_ind);
            curLick.smoothed_voltage = data(curSolnInd).smoothed_voltage(curLick.onset_ind:curLick.offset_ind);
            [maxVal,maxInd] = max(curLick.raw_voltage);
            curLick.maxInd = curLick.onset_ind + maxInd - 1;
            curLick.maxTime = data(curSolnInd).tvec(curLick.maxInd);
            curLick.maxVal = data(curSolnInd).raw_voltage(curLick.maxInd);
            curLick.certainty = mean(lickProbs(curSolnInd+1,curLick.onset_ind:curLick.offset_ind) - lickProbs(1,curLick.onset_ind:curLick.offset_ind));
            lickCount(curSolnInd) = lickCount(curSolnInd) + 1;
            licks{curSolnInd}(lickCount(curSolnInd)) = curLick;
            curLick = struct();
            curSolnInd = [];
            inLick = 0;
        else
            continue;
        end
    else
        if (inLick)
            continue;
        else
            inLick=1;
            curSolnInd = maxVec(t)-1;
            curLick.onset_ind = t;
            curLick.onset = data(curSolnInd).tvec(t);
            curLick.onsetVal = data(curSolnInd).raw_voltage(t);
            curLick.solution = data(curSolnInd).solution;
        end
    end
end
for i=1:length(licks)
    [refinedlicks] = refineOnsetsOffsets(licks{i},data(i));
    refinedLicks{i} = refinedlicks;
end
for i=1:length(refinedLicks)
    nValid = 0;
    for j=1:length(refinedLicks{i})
        if (refinedLicks{i}(j).duration < .01 || refinedLicks{i}(j).duration > .125)
            continue;
        else
            nValid = nValid+1;
            validLicks{i}(nValid) = refinedLicks{i}(j);
        end 
    end
end
if (strcmp(manualCorrection,'yes'))
    if (isempty(validLicks))
        goodlicks = [];
    else
        for i=1:length(licks)
            %[~,sortedInds] = sort([licks{i}.certainty]);
            %[newlicks] = manualLickID(data(i),licks{i}(sortedInds));
            %[~,sortedOrder] = sort([newlicks.onset]);
            %goodlicks{i} = newlicks(sortedOrder);
            if (i == 1)
                goodlicks{i} = manualLickID(data(i),validLicks{i},savedir,validLicks,i);
            else
                goodlicks{i} = manualLickID(data(i),validLicks{i},savedir,goodlicks,i);
            end
        end
    end
else
    goodlicks = validLicks;
end
end

