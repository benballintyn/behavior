function [] = analzyeParamSlices(basedir,modelDir,subdir)
[paramVals,paramNames,scores] = getParamVals(basedir,modelDir,subdir);
info = load([basedir '/' modelDir '/' subdir '/' modelDir '.mat']); info=info.(modelDir);
for i=1:length(paramNames)
    for j=1:length(paramNames)
        if (i>=j)
            continue;
        else
            F = scatteredInterpolant(paramVals(:,i),paramVals(:,j),scores');
            xmin = info.lower_bounds(i); xmax = info.upper_bounds(i);
            ymin = info.lower_bounds(j); ymax = info.upper_bounds(j);
            dx = (xmax-xmin)/1000;
            dy = (ymax-ymin)/1000;
            [Xq,Yq] = meshgrid(xmin:dx:xmax,ymin:dy:ymax);
            Vq = F(Xq,Yq);
            figure;
            imagesc([xmin xmax],[ymin ymax],Vq); colormap(jet); colorbar()
            set(gca,'ydir','normal')
            xlabel(paramNames{i}); ylabel(paramNames{j})
            hold on; plot(paramVals(:,i),paramVals(:,j),'.k','markersize',10)
        end
    end
end
end

