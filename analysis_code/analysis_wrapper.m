function [data,licks,bouts,switches] = analysis_wrapper(basedir,date,animal,tracetype,thresh,removelicks,doplot)
savedir = [basedir '/' date '/' animal];
nDV = 6;
minDur = .05;
maxDur = .15;
ampThresh = 20;
volumes = load([basedir '/' date '/' animal '/volumes.mat']); volumes=volumes.volumes;
solutions = load([basedir '/' date '/' animal '/solutions.mat']); solutions=solutions.solutions;
data = read_datafiles(basedir,date,animal);
temp = find(data(1).tvec > 5*60); endbaseline = temp(1);
save([basedir '/' date '/' animal '/data.mat'],'data','-mat')
licks1 = find_licks4(data(1).(tracetype),data(1).tvec,nDV,endbaseline);
licks2 = find_licks4(data(2).(tracetype),data(2).tvec,nDV,endbaseline);
switch removelicks
    case 'yes'
        goodlicks1 = remove_licks(licks1,thresh,ampThresh,minDur,maxDur);
        goodlicks2 = remove_licks(licks2,thresh,ampThresh,minDur,maxDur);
        licks{1} = goodlicks1;
        licks{2} = goodlicks2;
        disp([num2str(length(licks{1})) ' licks detected for solution ' data(1).solution])
        disp([num2str(length(licks{2})) ' licks detected for solution ' data(2).solution])
        bouts1 = get_lick_bouts2(licks{1});
        bouts2 = get_lick_bouts2(licks{2});
        bouts{1} = bouts1;
        bouts{2} = bouts2;
        disp([num2str(length(bouts{1})) ' lick bouts detected for solution ' data(1).solution])
        disp([num2str(length(bouts{2})) ' lick bouts detected for solution ' data(2).solution])
        save([basedir '/' date '/' animal '/licks_' tracetype '_filtered.mat'],'licks','-mat')
        save([basedir '/' date '/' animal '/lick_bouts.mat'],'bouts','-mat')
        switch doplot
            case 'yes'
                figure;
                subplot(2,1,1); plot(data(1).tvec,data(1).raw_voltage); 
                title(['Solution ' solutions(1)]); xlabel('Time (s)'); ylabel('Voltage (mV)')
                subplot(2,1,2); plot(data(2).tvec,data(2).raw_voltage); 
                title(['Solution ' solutions(2)]); xlabel('Time (s)'); ylabel('Voltage (mV)')
                savefig([savedir '/raw_voltage_plot.fig'])
                saveas(gcf,[savedir '/raw_voltage_plot.tif'],'tiffn')
                plot_trace_with_licks(data,tracetype,licks,bouts,savedir);
                plot_bout_lengths(bouts,solutions,savedir)
        end
        [lin_bouts,switches] = get_switches(bouts,solutions,savedir,doplot);
        save([basedir '/' date '/' animal '/switches.mat'],'switches','-mat')
    case 'no'
        licks{1} = licks1;
        licks{2} = licks2;
        save([basedir '/' date '/' animal '/licks_' tracetype '_unfiltered.mat'],'licks','-mat')
end
end

