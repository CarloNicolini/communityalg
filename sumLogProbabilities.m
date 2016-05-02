function [stop, logP] = sumLogProbabilities(nextLogP, logP)
%SUMLOGPROBABILITIES Utility function for the computation of hypergeometric
%distributions
if(nextLogP == 0)
    stop = true;
else
    stop = false;
    if(nextLogP > logP)
        common = nextLogP;
        diffExponent = logP - common;
    else
        common = logP;
        diffExponent = nextLogP - common;
    end
    logP = common + ((log(1 + 10.^diffExponent)) / log(10));
    
    if(nextLogP - logP > -4)
        stop = true;
    end
end