function [metadata] = makeMetadata(basedir,date,animal,lSln,lChn,lCnsm,rSln,rChn,rCnsm,mSln,mChn,mCnsm,w)
metadata.savedir = [basedir '/' date '/' animal];
metadata.leftSolution = lSln;
metadata.leftChannel = lChn;
metadata.leftConsumed = lCnsm;
metadata.rightSolution = rSln;
metadata.rightChannel = rChn;
metadata.rightConsumed = rCnsm;
metadata.middleSolution = mSln;
metadata.middleChannel = mChn;
metadata.middleConsumed = mCnsm;
metadata.weight = w;
save([basedir '/' date '/' animal '/metadata.mat'],'metadata','-mat')
end

