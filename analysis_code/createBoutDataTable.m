function [boutTable] = createBoutDataTable(animalObjs,varargin)
p = inputParser;
addRequired(p,'animalObjs')
addParameter(p,'boutType','normal',@ischar)
parse(p,animalObjs,varargin{:})

rat = cell(0);
date = cell(0);
day = [];
sex = cell(0);
box_side = cell(0);
channel = [];
solution = cell(0);
solutionNum = [];
alternative_solution = cell(0);
alternative_solutionNum = [];
palatability = [];
alternative_palatability = [];
onset = [];
offset = [];
duration = [];
back1_solution = cell(0);
back2_solution = cell(0);
back3_solution = cell(0);
back1_duration = [];
back2_duration = [];
back3_duration = [];
nlicks = [];
totalPrevLicks = [];
totalPrevDuration = [];
boutNumber = [];
early = [];
late = [];
count = 0;
for i=1:length(animalObjs)
    for j=1:length(animalObjs(i).dates)
        dayLinBouts = [];
        if (strcmp(p.Results.boutType,'normal'))
            for k=1:length(animalObjs(i).bouts{j})
                dayLinBouts = [dayLinBouts animalObjs(i).bouts{j}{k}];
            end
        elseif (strcmp(p.Results.boutType,'200'))
            for k=1:length(animalObjs(i).bouts200{j})
                dayLinBouts = [dayLinBouts animalObjs(i).bouts200{j}{k}];
            end
        else
            error('boutType not recognized')
        end
        dayOnsets = [dayLinBouts.onset];
        [~,sortedInds] = sort(dayOnsets);
        dayLinBouts = dayLinBouts(sortedInds);
        if (strcmp(p.Results.boutType,'normal'))
            curDayBouts = animalObjs(i).bouts{j};
        elseif (strcmp(p.Results.boutType,'200'))
            curDayBouts = animalObjs(i).bouts200{j};
        end
        for k=1:length(curDayBouts)
            if (isempty(curDayBouts{k}))
                continue;
            end
            curSolnNum = animalObjs(i).solnConverter(curDayBouts{k}(1).solution);
            alt_ind = setdiff(1:length(curDayBouts),k);
            if (~isempty(curDayBouts{alt_ind}))
                alt_solution = curDayBouts{alt_ind}(1).solution;
                altSolnNum = animalObjs(i).solnConverter(alt_solution);
            else
                altSolnNum = setdiff(animalObjs(i).solutions(j,:),curSolnNum);
                if (isempty(altSolnNum))
                    alt_solution = curDayBouts{k}(1).solution;
                    altSolnNum = animalObjs(i).solnConverter(alt_solution);
                else
                    alt_solution = animalObjs(i).solnConverter(altSolnNum);
                end
            end
            for l=1:length(curDayBouts{k})
                count = count + 1;
                %rat = [rat animalObjs(i).name];
                rat{count} = animalObjs(i).name;
                %date = [date animalObjs(i).dates{j}];
                date{count} = animalObjs(i).dates{j};
                day = [day j];
                %sex = [sex animalObjs(i).sex];
                sex{count} = animalObjs(i).sex;
                %box_side = [box_side animalObjs(i).bouts{j}{k}(l).box_side];
                box_side{count} = curDayBouts{k}(l).box_side;
                channel = [channel curDayBouts{k}(l).channel];
                %solution = [solution animalObjs(i).bouts{j}{k}(l).solution];
                solution{count} = curDayBouts{k}(l).solution;
                solutionNum = [solutionNum curSolnNum];
                alternative_solutionNum = [alternative_solutionNum altSolnNum];
                %alternative_solution = [alternative_solution alt_solution];
                alternative_solution{count} = alt_solution;
                if (strcmp(p.Results.boutType,'normal'))
                    %palatability = [palatability animalObjs(i).Palatabilities{j}{k}(l)];
                    palatability = [palatability animalObjs(i).relativePalatabilitiesLicks(curSolnNum)];
                    %alternative_palatability = [alternative_palatability animalObjs(i).meanAlternativePalatabilities{j}{k}(l)];
                    alternative_palatability = [alternative_palatability animalObjs(i).relativePalatabilitiesLicks(altSolnNum)];
                elseif (strcmp(p.Results.boutType,'200'))
                    %palatability = [palatability animalObjs(i).Palatabilities200{j}{k}(l)];
                    palatability = [palatability animalObjs(i).relativePalatabilitiesLicks(curSolnNum)];
                    %alternative_palatability = [alternative_palatability animalObjs(i).meanAlternativePalatabilities200{j}{k}(l)];
                    alternative_palatability = [alternative_palatability animalObjs(i).relativePalatabilitiesLicks(altSolnNum)];
                else
                    error('boutType not recognized')
                end
                onset = [onset curDayBouts{k}(l).onset];
                offset = [offset curDayBouts{k}(l).offset];
                duration = [duration curDayBouts{k}(l).duration];
                nlicks = [nlicks curDayBouts{k}(l).nlicks];
                isEarly = curDayBouts{k}(l).onset < animalObjs(i).bout_split_time;
                early = [early isEarly];
                late = [late ~isEarly];
                curBoutLinInd = find([dayLinBouts.onset] == curDayBouts{k}(l).onset);
                if (length(curBoutLinInd) > 1)
                    disp([num2str(i) ' ' num2str(j) ' ' num2str(k) ' ' num2str(l) '        ' num2str(curBoutLinInd)])
                    error('More than one bout index found with the same onset')
                end
                if (curBoutLinInd > 3)
                    back1_solution{count} = dayLinBouts(curBoutLinInd-1).solution;
                    back1_duration = [back1_duration dayLinBouts(curBoutLinInd-1).duration];
                    back2_solution{count} = dayLinBouts(curBoutLinInd-2).solution;
                    back2_duration = [back2_duration dayLinBouts(curBoutLinInd-2).duration];
                    back3_solution{count} = dayLinBouts(curBoutLinInd-3).solution;
                    back3_duration = [back3_duration dayLinBouts(curBoutLinInd-3).duration];
                elseif (curBoutLinInd == 3)
                    back1_solution{count} = dayLinBouts(curBoutLinInd-1).solution;
                    back1_duration = [back1_duration dayLinBouts(curBoutLinInd-1).duration];
                    back2_solution{count} = dayLinBouts(curBoutLinInd-2).solution;
                    back2_duration = [back2_duration dayLinBouts(curBoutLinInd-2).duration];
                    back3_solution{count} = nan;
                    back3_duration = [back3_duration nan];
                elseif (curBoutLinInd == 2)
                    back1_solution{count} = dayLinBouts(curBoutLinInd-1).solution;
                    back1_duration = [back1_duration dayLinBouts(curBoutLinInd-1).duration];
                    back2_solution{count} = nan;
                    back2_duration = [back2_duration nan];
                    back3_solution{count} = nan;
                    back3_duration = [back3_duration nan];
                else
                    back1_solution{count} = nan;
                    back1_duration = [back1_duration nan];
                    back2_solution{count} = nan;
                    back2_duration = [back2_duration nan];
                    back3_solution{count} = nan;
                    back3_duration = [back3_duration nan];
                end
                if (curBoutLinInd > 1)
                    totalPrevLicks = [totalPrevLicks sum([dayLinBouts(1:curBoutLinInd-1).nlicks])];
                    totalPrevDuration = [totalPrevDuration sum([dayLinBouts(1:curBoutLinInd-1).duration])];
                else
                    totalPrevLicks = [totalPrevLicks 0];
                    totalPrevDuration = [totalPrevDuration 0];
                end
                boutNumber = [boutNumber curBoutLinInd];
            end
        end
    end
end
rat = rat';
date = date';
day = day';
sex = sex';
box_side = box_side';
channel = channel';
solution = solution';
solutionNum = solutionNum';
alternative_solution = alternative_solution';
alternative_solutionNum = alternative_solutionNum';
palatability = palatability';
alternative_palatability = alternative_palatability';
onset = onset';
offset = offset';
duration = duration';
nlicks = nlicks';
back1_solution = back1_solution';
back1_duration = back1_duration';
back2_solution = back2_solution';
back2_duration = back2_duration';
back3_solution = back3_solution';
back3_duration = back3_duration';
totalPrevLicks = totalPrevLicks';
totalPrevDuration = totalPrevDuration';
boutNumber = boutNumber';
early = early';
late = late';
%{
size(rat)
size(date)
size(day)
size(sex)
size(box_side)
size(channel)
size(solution)
size(solutionNum)
size(alternative_solution)
size(alternative_solutionNum)
size(palatability)
size(alternative_palatability)
size(onset)
size(offset)
size(duration)
size(nlicks)
size(back1_solution)
size(back1_duration)
size(back2_solution)
size(back2_duration)
size(back3_solution)
size(back3_duration)
size(totalPrevLicks)
size(totalPrevDuration)
size(boutNumber)
%}
boutTable = table(rat,date,day,sex,box_side,channel,solution,solutionNum,alternative_solution,...
                  alternative_solutionNum,palatability,alternative_palatability,...
                  onset,offset,duration,nlicks,back1_solution,back1_duration,...
                  back2_solution,back2_duration,back3_solution,back3_duration,...
                  totalPrevLicks,totalPrevDuration,boutNumber,early,late);
end

