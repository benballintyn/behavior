function [ilis] = getILIs(licks)
ilis = zeros(1,length(licks)-1);
for i=1:(length(licks)-1)
    ilis(i) = licks(i+1).onset - licks(i).offset;
end
end

