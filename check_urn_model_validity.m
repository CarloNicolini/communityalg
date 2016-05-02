
function [argsok] = check_urn_model_validity(p,pi,m,mi)
%CHECK_URN_MODEL_VALIDITY
%
% This function check that the arguments provided to compute_surprise
% function are ok in the sense that they represent a realistic urn model:
% An urn with p balls inside, pi white and p-pi black, and m balls are
% drawm and one looks for the probability that at least mi of them are
% white.
%
% Inputs:   p, total number of balls in the urn
%           pi, number of white balls in the urn
%           m, number of drawn balls
%           mi, number of white balls in the drawn balls
%
%   Carlo Nicolini, Istituto Italiano di Tecnologia (2016).

argsok = true;
if ( mi<0 )
    error('Integer overflow: negative mi');
end

if ( m<0 )
    
    error('Integer overflow: negative m');
end

if ( pi<0 )
    
    error('Integer overflow: negative pi');
end

if ( p<0 )
    error('Integer overflow: negative p');
end

if (pi>p)
    error('Error: impossible urn model p<pi');
end

if (mi>m)
    error('Error: impossible urn model m<mi');
end

if ((m-mi) > (p-pi) )
    
    error('Error: impossible urn model (p-pi)<(m-mi)');
end

if (pi>p)
    error('Error: impossible urn model pi>p');
end

if (mi > pi)
    error('Error: impossible urn modelmi>pi');
end

if (mi==0)
    argsok=false;
    return;
end

if (pi==0)
    argsok=false;
    return;
end

if ((m-mi) == (p-pi) )
    argsok=false;
    return;
end

argsok=true;
end