function [fracSame] = compareSequences(seq1,seq2)
l1 = length(seq1);
l2 = length(seq2);
if (l1 ~= l2)
    error('Sequences are not the same length')
else
    sameCount=0;
    for i=1:l1
        if (seq1(i) == seq2(i))
            sameCount = sameCount+1;
        end
    end
    fracSame = sameCount/l1;
end
end

