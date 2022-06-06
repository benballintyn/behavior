function [soln] = getSolutionType(fname)
% For new file naming scheme starting with mt1 and mt2
if (contains(fname,'ch'))
    idx = strfind(fname,'_');
    if (length(idx) == 2)
        soln = fname(idx(2)-1);
    elseif (length(idx) == 3)
        soln = fname((idx(1)+5):(idx(2)-1));
        num = fname(idx(3)-1);
        if (str2num(num) > 2 && mod(str2num(num),2) == 1)
            num = '1';
        elseif (str2num(num) > 2 && mod(str2num(num),2) == 2)
            num = '2';
        end
        soln = [soln num];
    else
        error('file name does not fit the correct format')
    end
else % For original file name scheme
    idx = strfind(fname,'_');
    if (length(idx) == 2)
        soln = fname(idx(2)-1);
    elseif (length(idx) == 3)
        soln = fname((idx(1)+5):(idx(2)-1));
        num = fname(idx(3)-1);
        soln = [soln num];
    else
        error('file name does not fit the correct format')
    end
end
end

