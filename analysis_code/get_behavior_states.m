function [states,stateseq] = get_behavior_states(bouts,tvec)
% get_behavior_states
%   This function takes as input the bouts and timestamp vector for one
%   session and returns 2 representations of the state sequence of the
%   animal. For a session with n solutions, there are n+1 possible states.
%   Sampling at the n solutions represents n of the states and the last
%   state represents the animal not sampling any of the solutions (the
%   animal is 'wandering')
%
%   Inputs:
%       bouts - cell array of struct arrays of length # of solutions
%
%       tvec - vector of timestamps for each timestep
%
%   Outpus:
%       states - (n+1)x36000 matrix where each column has one entry that is
%                1 indicating the animal is in that state
%
%       stateseq - vector of length 36000 where each entry gives the state
%                  of the animal during that 100ms window
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
   if (length(state) > 1)
       error(['Bout overlap at t = ' num2str(tvec(i))])
   end
   stateseq(i) = state;
end
end

