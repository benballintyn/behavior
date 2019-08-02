function [] = plotLicks(dataTrace,tvec,licks,mini_bouts,bouts)
onsets = [licks.onset];
onsetVals = [licks.onsetVal];
offsets = [licks.offset];
offsetVals = [licks.offsetVal];
maxVals = [licks.maxVal];
maxTimes = [licks.maxTime];
maxy = max(dataTrace);
figure;
plot(tvec,dataTrace);
hold on;
plot(onsets+.0002,onsetVals,'r*')
plot(offsets-.0002,offsetVals,'b*')
plot(maxTimes,maxVals,'k*')
%{
for i=1:length(mini_bouts)
    plot([[mini_bouts.onset];[mini_bouts.onset]],[zeros(1,length(mini_bouts));ones(1,length(mini_bouts))*maxy],'m','linew',1)
    plot([[mini_bouts.offset];[mini_bouts.offset]],[zeros(1,length(mini_bouts));ones(1,length(mini_bouts))*maxy],'m','linew',1)
end
%}

for i=1:length(bouts)
    plot([[bouts.onset];[bouts.onset]],[zeros(1,length(bouts));ones(1,length(bouts))*maxy],'k')
    plot([[bouts.offset];[bouts.offset]],[zeros(1,length(bouts));ones(1,length(bouts))*maxy],'g')
end
end

