function [] = plotLicks(dataTrace,tvec,licks)
onsets = [licks.onset];
onsetVals = [licks.onsetVal];
offsets = [licks.offset];
offsetVals = [licks.offsetVal];
maxVals = [licks.maxVal];
maxTimes = [licks.maxTime];
figure;
plot(tvec,dataTrace);
hold on;
plot(onsets+.002,onsetVals,'r*')
plot(offsets-.002,offsetVals,'b*')
plot(maxTimes,maxVals,'g*')
end

