animalInfo = load('analyzed_data/animalInfo.mat'); animalInfo=animalInfo.animalInfo;

for i=1:length(animalInfo)
    for j=1:length(animalInfo(i).dates)
        bouts = load(['analyzed_data/' animalInfo(i).dates{j} '/' animalInfo(i).animal '/bouts.mat']);
        bouts = bouts.bouts;
        for k=1:length(bouts)
            for l=1:length(bouts{k})
                if (bouts{k}(l).nlicks ~= length(bouts{k}(l).licks))
                    bouts{k}(l).nlicks = length(bouts{k}(l).licks);
                    disp(['Corrected nlicks for ' animalInfo(i).animal ' ' animalInfo(i).dates{j} ' ' num2str(k) ' ' num2str(l)])
                end
            end
        end
        save(['analyzed_data/' animalInfo(i).dates{j} '/' animalInfo(i).animal '/bouts.mat'],'bouts','-mat')
    end
end