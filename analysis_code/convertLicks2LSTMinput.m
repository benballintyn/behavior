function [Xtr,Ytr,Xval,Yval,Xtst,Ytst] = convertLicks2LSTMinput(basedir,dates,animals,wndwLength,fracOverlap)
totalDataPoints = 0;
dates
for i=1:length(dates)
    for j=1:length(animals)
        curdir = [basedir '/' dates{i} '/' animals{j}];
        if (exist(curdir,'dir') && exist([curdir '/licks.mat'],'file'))
            disp([dates{i} '/' num2str(animals{j})])
            data=load([curdir '/data.mat']); data=data.data;
            licks=load([curdir '/licks.mat']); licks=licks.licks;
            dataLength = length(data(1).raw_voltage);
            lickVec = ones(1,dataLength);
            overlapLength = floor(fracOverlap*wndwLength);
            startOffset = ceil(fracOverlap*wndwLength);
            npts = 1 + (dataLength - wndwLength)/overlapLength;
            for k=1:length(licks)
                if (~isempty(licks{k}))
                    onsets = [licks{k}.onset_ind];
                    offsets = [licks{k}.offset_ind];
                    for l=1:length(onsets)
                        lickVec(onsets(l):offsets(l)) = k+1;
                    end
                end
            end
            for k=1:npts
                s = (k-1)*(wndwLength - startOffset) + 1;
                e = s + wndwLength - 1;
                x = zeros(length(data),wndwLength);
                y = zeros(1,wndwLength);
                for l=1:length(data)
                    x(l,:) = data(l).raw_voltage(s:e);
                    y = lickVec(s:e);
                end
                totalDataPoints = totalDataPoints+1;
                X{totalDataPoints,1} = x;
                Y{totalDataPoints,1} = categorical(y);
            end
        end
    end
end
randInds = randperm(length(X));
X = X(randInds);
Y = Y(randInds);
tstValFrac = floor(.05*length(X));
trEnd = length(X) - 2*tstValFrac;
valEnd = trEnd + tstValFrac;
tstEnd = valEnd + tstValFrac;
Xtr = X(1:trEnd);
Ytr = Y(1:trEnd);
Xval = X(trEnd+1:valEnd);
Yval = Y(trEnd+1:valEnd);
Xtst = X(valEnd+1:tstEnd);
Ytst = Y(valEnd+1:tstEnd);
end

