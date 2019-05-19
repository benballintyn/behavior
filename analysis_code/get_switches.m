function [linear_bouts,switches] = get_switches(bouts,solutions,savedir,doplot)
ind1 = 1;
ind2 = 1;
curInd = 1;
nswitches = 0;
done = 0;
while (~done)
    if (ind1 > length(bouts{1}) && ind2 <= length(bouts{2}))
        linear_bouts(curInd) = bouts{2}(ind2);
        solns(curInd) = 1;
        %linear_bouts(curInd).solution = solutions(2);
        ind2 = ind2 + 1;
    elseif (ind2 > length(bouts{2}) && ind1 <= length(bouts{1}))
        linear_bouts(curInd) = bouts{1}(ind1);
        solns(curInd) = 0;
        %linear_bouts(curInd).solution = solutions(1);
        ind1 = ind1 + 1;
    elseif (ind1 > length(bouts{1}) && ind2 > length(bouts{2}))
        done = 1;
        continue;
    else
        if (bouts{1}(ind1).onset < bouts{2}(ind2).onset)
            linear_bouts(curInd) = bouts{1}(ind1);
            solns(curInd) = 0;
            %linear_bouts(curInd).solution = solutions(1);
            ind1 = ind1 + 1;
        else
            linear_bouts(curInd) = bouts{2}(ind2);
            solns(curInd) = 1;
            %linear_bouts(curInd).solution = solutions(2);
            ind2 = ind2 + 1;
        end
    end
    curInd = curInd + 1;
end

switch doplot
    case 'yes'
        figure; hold on;
        set(gca,'YTick',[0 1],'YTickLabel',{solutions(1), solutions(2)})
        for i=1:length(linear_bouts)
            onset = linear_bouts(i).onset;
            offset = linear_bouts(i).offset;
            soln = solns(i);
            if (i == 1)
                lastSoln = soln;
            end
            plot([onset offset],[soln soln],'LineWidth',5,'Color','k')
            if (lastSoln ~= soln)
                nswitches=nswitches+1;
                switches(nswitches).start = linear_bouts(i-1).offset;
                switches(nswitches).start_soln = linear_bouts(i-1).lastSoln;
                switches(nswitches).end = linear_bouts(i).onset;
                switches(nswitches).end_soln = soln;
                lastOffset = linear_bouts(i-1).offset;
                plot([lastOffset onset],[lastSoln soln],'LineWidth',1,'Color','k')
                lastSoln = soln;
            end
        end
        xlabel('Time (s)')
        ylabel('Solution ID')
        set(gcf,'Position',[10 10 1000 400])
        savefig([savedir '/switches.fig'])
        saveas(gcf,[savedir '/switches.tif'],'tiffn')
    case 'no'
        for i=1:length(linear_bouts)
            onset = linear_bouts(i).onset;
            offset = linear_bouts(i).offset;
            soln = solns(i);
            if (i == 1)
                lastSoln = soln;
            end
            if (lastSoln ~= soln)
                nswitches=nswitches+1;
                switches(nswitches).start = linear_bouts(i-1).offset;
                switches(nswitches).start_soln = lastSoln;
                switches(nswitches).end = linear_bouts(i).onset;
                switches(nswitches).end_soln = soln;
                lastOffset = linear_bouts(i-1).offset;
                lastSoln = soln;
            end
        end
end
end

