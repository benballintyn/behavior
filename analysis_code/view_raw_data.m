function [] = view_raw_data(basedir,date,animal)
datadir = [basedir '/' date '/' animal];
datafiles = dir([datadir '/' date animal '*.txt']);

for i=1:length(datafiles)
    raw_data = load([basedir '/' date '/' animal '/' datafiles(i).name],'-ascii');
    tsteps(:,i) = raw_data(:,2);
    voltages(:,i) = raw_data(:,1);
end
figure;
for i=1:length(datafiles)
    plot(tsteps(:,i)-tsteps(1,i),voltages(:,i));
    hold on;
end
hold off;
set(gcf,'Position',[100,50,1400,600])
end

