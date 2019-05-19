function [] = plot_trace_with_licks(data,tracetype,licks,bouts,savedir)
figure;
plotheight = (2/3)*max(data(1).(tracetype));
subplot(2,1,1); plot(data(1).tvec,data(1).(tracetype))
hold on; 
plot([licks{1}.onset],[licks{1}.onset_val],'*')
plot([licks{1}.offset],[licks{1}.offset_val],'*')
minData = min(data(1).(tracetype));
maxData = max(data(1).(tracetype));
for i=1:length(bouts{1})
    xval1 = bouts{1}(i).onset;
    xval2 = bouts{1}(i).offset;
    plot([xval1 xval1],[minData maxData],'k')
    plot([xval2 xval2],[minData maxData],'k')
end
xlabel('Time (s)')
ylabel('Filtered Voltage')
title(['Solution ' data(1).solution ' # of licks = ' num2str(length(licks{1}))])

plotheight = (2/3)*max(data(2).(tracetype));
subplot(2,1,2); plot(data(2).tvec,data(2).(tracetype))
hold on;
plot([licks{2}.onset],[licks{2}.onset_val],'*')
plot([licks{2}.offset],[licks{2}.offset_val],'*')
minData = min(data(2).(tracetype));
maxData = max(data(2).(tracetype));
for i=1:length(bouts{2})
    xval1 = bouts{2}(i).onset;
    xval2 = bouts{2}(i).offset;
    plot([xval1 xval1],[minData maxData],'k')
    plot([xval2 xval2],[minData maxData],'k')
end
xlabel('Time (s)')
ylabel('Filtered Voltage')
title(['Solution ' data(2).solution ' # of licks = ' num2str(length(licks{2}))])
savefig([savedir '/licks_with_bouts_plot.fig'])
saveas(gcf,[savedir '/licks_with_bouts_plot.tif'],'tiffn')
end

