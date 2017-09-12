function phi = phicoeff(memb_a, memb_b)

% Implementation of the phi coefficient as in 
% "Network community structure alterations in adult schizophrenia:
% "Identification and localization of alterations"
% 
% Lerman-Sinkoff, Dov B. Barch, Deanna M.
% Neuroimage, Clinical (2016)
%
n = length(memb_a);
if (length(memb_b) ~= n)
    error('Membership vectors have different length, cannot compute phi.');
end

memb_a=memb_a(:)'; % force it to be a row vector
memb_b=memb_b(:)'; % force it to be a row vector

% The first matrix C is an ncommsX lenght(ci) matrix of ones and zeros. 
% Each row holds a 1's where ci is equal to one of the unique value.
CijA = double(bsxfun(@eq,memb_a, unique(memb_a)'));
CijB = double(bsxfun(@eq,memb_b, unique(memb_b)'));
deltaijA = CijA'*CijA;
deltaijB = CijB'*CijB;

[X,Y] = meshgrid(1:n,1:n);
phi = arrayfun(@(x,y)(corrcoefXY(deltaijA(:,x),deltaijB(:,y))), X, Y);

function phi = corrcoefXY(x,y)
    phi = corrcoef(x,y);
    phi = phi(end,1);