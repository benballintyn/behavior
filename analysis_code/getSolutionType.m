function [soln] = getSolutionType(fname)
idx = strfind(fname,'_');
if (length(idx) == 2)
    soln = fname(idx(2)-1);
elseif (length(idx) == 3)
    soln = fname((idx(1)+1):(idx(2)-1));
    num = fname(idx(3)-1);
    soln = [soln num];
else
    error('file name does not fit the correct format')
end
end

