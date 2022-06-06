function [concentration] = getConcentrationFromSolnName(solnName)
if (contains(solnName,'A'))
    concentration = '.01M';
elseif (contains(solnName,'B'))
    concentration = '.1M';
elseif (contains(solnName,'C'))
    concentration = '1M';
elseif (contains(solnName,'H2O'))
    concentration = 'H2O';
else
    error('solnName not recognized')
end

