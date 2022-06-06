function [logNaClConcentration] = getLogNaClConcentration(solnName)
if (contains(solnName,'A'))
    logNaClConcentration = -2;
elseif (contains(solnName,'B'))
    logNaClConcentration = -1;
elseif (contains(solnName,'C'))
    logNaClConcentration = 0;
elseif (contains(solnName,'H2O'))
    logNaClConcentration = NaN;
else
    error('solnName not valid')
end
end

