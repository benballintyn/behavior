function [] = plot_bout_lengths(bouts,solutions,savedir)
[h,p] = kstest2([bouts{1}.nlicks],[bouts{2}.nlicks]);
figure;
data = [[bouts{1}.nlicks] [bouts{2}.nlicks]];
g = [zeros(1,length(bouts{1})) ones(1,length(bouts{2}))];
boxplot(data,g,'Labels',{solutions(1), solutions(2)})
ylabel('Licks / lick bout')
xlabel('Solution ID')
title(['p = ' num2str(p)])
savefig([savedir '/bout_length_boxplot.fig'])
saveas(gcf,[savedir '/bout_length_boxplot.tif'],'tiffn')

figure; hold on;
middlePts1 = ([bouts{1}.offset] + [bouts{1}.onset])/2;
middlePts2 = ([bouts{2}.offset] + [bouts{2}.onset])/2;
plot(middlePts1,[bouts{1}.nlicks],'LineWidth',3)
plot(middlePts2,[bouts{2}.nlicks],'LineWidth',3)
xlabel('Time (s)')
ylabel('Licks / lick bout')
legend({solutions(1),solutions(2)})
savefig([savedir '/bout_length_vs_time.fig'])
saveas(gcf,[savedir '/bout_length_vs_time.tif'],'tiffn')
end

