% reduce data
oldDataDir = '/home/ben/phd/behavior/analyzed_data';
newDataDir = '/home/ben/phd/behavior/reduced_data';
mkdir(newDataDir)
animalInfo = load('analyzed_data/animalInfo.mat'); animalInfo=animalInfo.animalInfo;

for i=1:length(animalInfo)
    curRat = animalInfo(i).animal;
    for j=1:length(animalInfo(i).dates)
        curDate = animalInfo(i).dates{j};
        if (~exist([newDataDir '/' curDate],'dir'))
            mkdir([newDataDir '/' curDate])
        end
        curDir = [newDataDir '/' curDate '/' curRat];
        mkdir(curDir)
        
        src = [oldDataDir '/' curDate '/' curRat '/licks.mat'];
        tgt = [newDataDir '/' curDate '/' curRat '/licks.mat'];
        copyfile(src,tgt)
        
        src = [oldDataDir '/' curDate '/' curRat '/bouts.mat'];
        tgt = [newDataDir '/' curDate '/' curRat '/bouts.mat'];
        copyfile(src,tgt)
        
        src = [oldDataDir '/' curDate '/' curRat '/metadata.mat'];
        tgt = [newDataDir '/' curDate '/' curRat '/metadata.mat'];
        copyfile(src,tgt)
        
        src = [oldDataDir '/' curDate '/' curRat '/states.mat'];
        tgt = [newDataDir '/' curDate '/' curRat '/states.mat'];
        copyfile(src,tgt)
        
        src = [oldDataDir '/' curDate '/' curRat '/manual_flag.mat'];
        if (exist(src,'file'))
            tgt = [newDataDir '/' curDate '/' curRat '/manual_flag.mat'];
            copyfile(src,tgt)
        else
            disp([src ' does not exist'])
        end
    end
end