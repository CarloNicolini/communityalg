function v = membership_agreement(m1, m2)
%MEMBERSHIP_AGREEMENT Returns a vector with the node-membership agreement
%between two memberships
if length(m1)~=length(m2)
    error('Vectors must be the same size');
end

[sjacc,sp] = membership_similarity(m1,m2);
n = length(m1);
v=zeros(length(m1),1);
