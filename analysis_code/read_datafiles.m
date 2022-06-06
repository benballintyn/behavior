function [data] = read_datafiles(basedir,date,animal)
% read_datafiles
%   This function takes as input the relevant info about the date and
%   animal of the data to be extracted. It then reads the information in
%   the relevant data (.txt) files into a struct array of length nChannels
%   where nChannels is the number of channels (and number of data files)
%
%   Inputs:
%       basedir - top level directory containg all animal data
%
%       date - date in 'year/month/day' format of the data to be analyzed
%       e.g. '190808'
%
%       animal - name of animal to be analyzed. e.g. 'bb8'
%
%   Outputs:
%       data - struct array (with 1 entry per channel) containing raw and
%              filtered data
datadir = [basedir '/' date '/' animal];
datafiles = dir([datadir '/' date animal '*.txt']);
expt_start = 5; % baseline for 5 minutes in or 300 seconds
for i=1:length(datafiles)
    data(i).datafile = datafiles(i).name;
    data(i).solution = getSolutionType(data(i).datafile);
    raw_data = load([datadir '/' datafiles(i).name],'-ascii');
    raw_voltage = raw_data(:,1);
    timestamps = raw_data(:,2);
    tvec = timestamps - timestamps(1);
    data(i).tvec = tvec;
    tvecdiff = diff(tvec);
    [b,a] = butter(4,[30]/(.5*(1/mean(tvecdiff))),'low');
    v = filtfilt(b,a,raw_voltage);
    data(i).raw_voltage = raw_voltage;
    data(i).lowpass_voltage = v;
    [b,a] = butter(4,[40 80]/(.5*(1/mean(tvecdiff))),'stop');
    y = filtfilt(b,a,raw_voltage);
    %y = bandpass(raw_voltage,[100 200],1/mean(diff(tvec)));
    data(i).bandpass_voltage = y;
    temp = find(tvec > expt_start*60);
    expt_start_ind = temp(1);
    baseline = v(1:expt_start_ind);
    data(i).zlowpass = (v - mean(baseline))/std(baseline);
    baseline = y(1:expt_start_ind);
    data(i).zbandpass = (y - mean(baseline))/std(baseline);
    baseline = raw_voltage(1:expt_start_ind);
    data(i).zraw = (raw_voltage - mean(baseline))/std(baseline);
    data(i).ema_fast = getEMA(raw_voltage,.5);
    data(i).ema_slow = getEMA(raw_voltage,.25);
    data(i).smoothed_voltage = gaussSmooth(raw_voltage,11);
end
end

