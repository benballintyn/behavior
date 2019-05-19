function [data] = read_datafiles(basedir,date,animal)
datadir = [basedir '/' date '/' animal];
datafiles = dir([datadir '/' date animal '*.txt']);
expt_start = 5;
for i=1:length(datafiles)
    data(i).datafile = datafiles(i).name;
    raw_data = load([datadir '/' datafiles(i).name],'-ascii');
    raw_voltage = raw_data(:,1);
    timestamps = raw_data(:,2);
    tvec = timestamps - timestamps(1);
    data(i).tvec = tvec;
    tvecdiff = diff(tvec);
    [b,a] = butter(4,[20]/(.5*(1/mean(tvecdiff))),'low');
    v = filtfilt(b,a,raw_voltage);
    data(i).raw_voltage = raw_voltage;
    data(i).filtered_voltage = v;
    temp = find(tvec > expt_start*60);
    expt_start_ind = temp(1);
    baseline = raw_voltage(1:expt_start_ind);
    data(i).zscore = (raw_voltage - mean(baseline))/std(baseline);
end
end

