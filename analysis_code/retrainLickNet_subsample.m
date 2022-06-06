function [lickNet] = retrainLickNet_subsample(netType,dates,animals,batchSize)
if (length(dates) > 10)
    dateNumMult = ceil(length(dates)/10);
end
for i=1:(dateNumMult*5)
    disp(['Sample ' num2str(i) '/' num2str(dateNumMult*5)])
    randDates = randperm(length(dates));
    for j = 1:10
        sampledDates{j} = dates{randDates(j)};
    end
    sampledDates
    disp('creating input from data ...')
    [Xtr,Ytr,Xval,Yval,Xtst,Ytst] = convertLicks2LSTMinput('data',sampledDates,animals,3000,.5);
    disp('input created ...')
    nbatches = ceil(length(Xtr)/batchSize);
    i=1;
    initialLearnRate = .005;
    while i < nbatches
        try
            disp(['Batch ' num2str(i) ' / ' num2str(nbatches)])
            if (strcmp(netType,'2bottle'))
                lickNet = load(['analysis_code/bilstmLickNet2Bottle.mat']); lickNet=lickNet.bilstmLickNet;
                save('analysis_code/oldLickNet.mat','lickNet','-mat')
            elseif (strcmp(netType,'3bottle'))
                lickNet = load(['analysis_code/bilstmLickNet3Bottle.mat']); lickNet=lickNet.bilstmLickNet;
            else
                error('Wrong netType given')
            end
            trRandInds = randperm(length(Xtr)); trRandInds=trRandInds(1:batchSize);
            valRandInds = randperm(length(Xval)); valRandInds=valRandInds(1:length(Xval));
            options = trainingOptions('adam','MaxEpochs',50,'MiniBatchSize',256,...
                    'Shuffle','every-epoch','ExecutionEnvironment','gpu',...
                    'ValidationData',{Xval(valRandInds),Yval(valRandInds)},'ValidationFrequency',10,...
                    'ValidationPatience',10,'Plots','training-progress',...
                    'LearnRateSchedule','piecewise','LearnRateDropFactor',.75,...
                    'LearnRateDropPeriod',1,'InitialLearnRate',initialLearnRate/(i/2),'Plots','none');
            lickNet = trainNetwork(Xtr(trRandInds),Ytr(trRandInds),lickNet.Layers,options);
            if (strcmp(netType,'2bottle'))
                bilstmLickNet = lickNet;
                disp('saving ...')
                save('analysis_code/bilstmLickNet2Bottle.mat','bilstmLickNet','-mat')
            elseif (strcmp(netType,'3bottle'))
                disp('saving ...')
                save('analysis_code/bilstmLickNet3Bottle.mat','bilstmLickNet','-mat')
            else
                warning('netType for saving not recognized')
                error('netType for saving not recognized')
                break;
            end
            i = i+1;
        catch err
            err
            warning('That annoying bug happened')
            continue;
        end
    end
end
end

