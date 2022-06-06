function [boutLengthsVsPastLicks,boutLengthsVsPastBout] = getBoutLengthVsPastLicks(basedir,dates,animal)
count=0;
count2=0;
for i=1:length(dates)
    dataDir = [basedir '/' dates{i} '/' animal];
    bouts = load([dataDir '/bouts.mat']); bouts = bouts.bouts;
    for j=1:length(bouts)
        totalLicks = 0;
        for k=1:length(bouts{j})
            count=count+1;
            boutLengthsVsPastLicks(count,1) = totalLicks;
            boutLengthsVsPastLicks(count,2) = bouts{j}(k).duration;
            if (k > 1)
                count2=count2+1;
                boutLengthsVsPastBout(count2,1) = bouts{j}(k-1).duration;
                boutLengthsVsPastBout(count2,2) = bouts{j}(k).duration;
            end
            totalLicks = totalLicks + bouts{j}(k).nlicks;
        end
    end
end
end

