% This script is used to remove incorrectly labeled licks. This can be
% redundant (duplicate licks) or overlapping licks
clear all; close all;
% Load info structure containing names and dates of all animals
animalInfo = load('analyzed_data/animalInfo.mat'); animalInfo=animalInfo.animalInfo;
% Keep track of the total # of licks removed
removedCount = 0;
for i=1:length(animalInfo)
    for j=1:length(animalInfo(i).dates) % For each day of each animal
        % Load current data, bout, and lick data structures
        curdir = ['analyzed_data/' animalInfo(i).dates{j} '/' animalInfo(i).animal '/'];
        data = load([curdir 'data.mat']); data=data.data;
        bouts = load([curdir 'bouts.mat']); bouts=bouts.bouts;
        licks = load([curdir 'licks.mat']); licks=licks.licks;
        % Create a copy of the old bouts and licks
        save([curdir 'oldBouts.mat'],'bouts','-mat')
        save([curdir 'oldLicks.mat'],'licks','-mat')
        % Create a new cell array for the remaining licks
        newlicks = cell(size(licks));
        % Keep track of the number of licks removed for this day
        totalRemoved = 0;
        for k=1:length(bouts)
            for l=1:length(bouts{k})
                boutONind = bouts{k}(l).onset_ind; % current bout onset time index
                boutOFFind = bouts{k}(l).offset_ind; % current bout offset time index
                curBoutLicks = bouts{k}(l).licks; % extract all of the licks for this bout
                onsets = [curBoutLicks.onset];
                onsetVals = [curBoutLicks.onsetVal];
                offsets = [curBoutLicks.offset];
                offsetVals = [curBoutLicks.offsetVal];
                ilis = onsets(2:end)-offsets(1:end-1); % comput inter lick intervals
                bad = find(ilis < 0); % see if there are any negative (incorrect) ILIs
                if (~isempty(bad)) % if there are negative ILIs, remove bad licks
                    timeIndsON = [curBoutLicks(bad).onset_ind];
                    timeIndsOFF = [curBoutLicks(bad+1).offset_ind];
                    
                    disp([animalInfo(i).animal ' ' animalInfo(i).dates{j} ' ind ' num2str(k) ' bout ' num2str(l) ' # ' num2str(length(bad)) ' bad ' num2str(bad)])
                    
                    % Look for duplicate licks (licks whose both onsets and
                    % offsets are within 10 timepoints which corresponds to
                    % less than 10ms)
                    duplicateLicks = [];
                    for m=1:length(curBoutLicks)
                        for n=m:length(curBoutLicks)
                            if (m ~= n && (abs(curBoutLicks(m).onset_ind - curBoutLicks(n).onset_ind) < 10) && (abs(curBoutLicks(m).offset_ind - curBoutLicks(n).offset_ind) < 10))
                                duplicateLicks = [duplicateLicks n];
                            end
                        end
                    end
                    % Remove duplicate licks from the current bout
                    curBoutLicks = curBoutLicks(setdiff(1:length(curBoutLicks),duplicateLicks));
                    totalRemoved = totalRemoved + length(duplicateLicks);
                    
                    % Look for licks that are contained within another
                    % lick
                    containedLicks = [];
                    for m=1:(length(curBoutLicks)-1)
                        if (curBoutLicks(m).onset_ind <= curBoutLicks(m+1).onset_ind && curBoutLicks(m+1).offset_ind <= curBoutLicks(m).offset_ind)
                            containedLicks = [containedLicks m+1];
                        elseif (curBoutLicks(m+1).onset_ind <= curBoutLicks(m).onset_ind && curBoutLicks(m).offset_ind <= curBoutLicks(m+1).offset_ind)
                            containedLicks = [containedLicks m];
                        end
                    end
                    % Removed contained licks
                    curBoutLicks = curBoutLicks(setdiff(1:length(curBoutLicks),containedLicks));
                    totalRemoved = totalRemoved + length(containedLicks);
                    
                    % Look for consecutive licks that overlap
                    overlappingLicks = [];
                    for m=1:(length(curBoutLicks)-1)
                        if (curBoutLicks(m).offset_ind >= curBoutLicks(m+1).onset_ind)
                            overlappingLicks = [overlappingLicks m+1];
                        end
                    end
                    % Remove overlapping licks
                    curBoutLicks = curBoutLicks(setdiff(1:length(curBoutLicks),overlappingLicks));
                    totalRemoved = totalRemoved + length(overlappingLicks);
                    
                    onsets = [curBoutLicks.onset];
                    onsetVals = [curBoutLicks.onsetVal];
                    offsets = [curBoutLicks.offset];
                    offsetVals = [curBoutLicks.offsetVal];
                    ilis = onsets(2:end) - offsets(1:end-1);
                    bad = find(ilis < 0); % Check for negative ILIs in the remaining licks
                    % If there are still remaining negative ILIs plot the
                    % relevant data and licks
                    if (length(bad) > 0)
                        
                        timeIndsON = [curBoutLicks(bad).offset_ind];
                        timeIndsOFF = [curBoutLicks(bad+1).onset_ind];
                    
                        disp([animalInfo(i).animal ' ' animalInfo(i).dates{j} ' ind ' num2str(k) ' bout ' num2str(l) ' # ' num2str(length(bad)) ' bad ' num2str(bad)])
                        plot(data(k).tvec(boutONind:boutOFFind),data(k).raw_voltage(boutONind:boutOFFind))
                        hold on;
                        plot(onsets(bad),onsetVals(bad),'r*')
                        plot(onsets(bad+1),onsetVals(bad+1),'r*')
                        plot(offsets(bad),offsetVals(bad),'b*')
                        plot(offsets(bad+1),offsetVals(bad+1),'b*')
                        plot(onsets,onsetVals,'rx')
                        plot(offsets,offsetVals,'bx')
                        for m=1:length(curBoutLicks)
                            text((curBoutLicks(m).onset + curBoutLicks(m).offset)/2,975,num2str(m))
                        end
                        for m=1:length(bad)
                            plot([data(k).tvec(timeIndsON(m)) data(k).tvec(timeIndsON(m))],[0 1000],'g')
                            plot([data(k).tvec(timeIndsOFF(m)) data(k).tvec(timeIndsOFF(m))],[0 1000],'k')
                        end
                        uiwait;
                        
                        disp('============================================================================')
                    end
                    disp(['corrected #bad = ' num2str(length(bad)) ' total removed = ' num2str(totalRemoved)])
                end
                % Replace licks and associated metadata in the current bout
                bouts{k}(l).licks = curBoutLicks;
                bouts{k}(l).nlicks = length(curBoutLicks);
                bouts{k}(l).onset = bouts{k}(l).licks(1).onset;
                bouts{k}(l).onset_ind = bouts{k}(l).licks(1).onset_ind;
                bouts{k}(l).offset = bouts{k}(l).licks(end).offset;
                bouts{k}(l).offset_ind = bouts{k}(l).licks(end).offset_ind;
                bouts{k}(l).duration = bouts{k}(l).offset - bouts{k}(l).onset;
                % Replace licks with only the remaining licks in the bouts
                newlicks{k} = [newlicks{k} curBoutLicks];
            end
        end
        licks = newlicks;
        % Save the new licks and bout data structures
        save([curdir 'licks.mat'],'licks','-mat')
        save([curdir 'bouts.mat'],'bouts','-mat')
        removedCount = removedCount + totalRemoved;
    end
end