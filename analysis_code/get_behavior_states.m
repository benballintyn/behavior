function [states,stateseq] = get_behavior_states(bouts,tvec)
ntastes = length(bouts);
states = zeros(ntastes+1,length(tvec));
for i=1:length(bouts)
    for b=1:length(bouts{i})
        states(i+1,bouts{i}(b).onset_ind:bouts{i}(b).offset_ind) = 1;
    end
end
for i=1:size(states,2)
    if (sum(states(:,i)) == 0)
        states(1,i) = 1;
    end
end
for i=1:size(states,2)
   state = find(states(:,i)==1);
   stateseq(i) = state;
end
end

