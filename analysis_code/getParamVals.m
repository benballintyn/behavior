function [paramVals,paramNames,scores] = getParamVals(basedir,modelDir,subdir)
datadir = [basedir '/' modelDir '/' subdir];
files=dir(datadir);
count=0;
for i=1:length(files)
    if (~(strcmp(files(i).name,'.') || strcmp(files(i).name,'..')) && files(i).isdir)
        count=count+1;
        p = load([datadir '/' num2str(count) '/params.mat']); p=p.params;
        paramVals(count,:) = p;
        score = load([datadir '/' num2str(count) '/score.mat']); score=score.score;
        scores(count) = score;
    end
end

info = load([datadir '/' modelDir '.mat']); info=info.(modelDir);
paramNames = info.agentParams;
paramNames = [paramNames info.envParams];
end

