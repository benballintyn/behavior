function [T] = getTransitionMatrix(stateSeqs,nLocations)
nDays = size(stateSeqs,1);
ntimepoints = size(stateSeqs,2);
T = zeros(nDays,nLocations,nLocations);
nT = zeros(nDays,nLocations,nLocations);
for i=1:nDays
    for t=1:ntimepoints-1
        curState = stateSeqs(i,t);
        nextState = stateSeqs(i,t+1);
        nT(i,nextState,curState) = nT(i,nextState,curState)+1;
    end
    for s=1:nLocations
        T(i,:,s) = squeeze(nT(i,:,s))./sum(squeeze(nT(i,:,s)));
    end
end
end

