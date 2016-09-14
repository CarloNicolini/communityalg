function write_brainet(filename, nodesXYZ, nodesColor, nodesSize, nodesLabel)
% write_brainet
% write a BrainetViewer compatible .node file for visualization
if nargin<5 || isempty(nodesLabel)
    nodesLabel = repmat('-',[size(nodesXYZ,1),1]);
end
n = size(nodesXYZ,1);
if length(nodesColor)~= n || length(nodesSize)~= n
    error('unconsistent number of nodes');
end

% create the table
T=table;
T.XYZ = nodesXYZ;
T.color = reshape(nodesColor(:),[n,1]);
T.size = reshape(nodesSize(:),[n,1]);

T.nodesLabel = nodesLabel;
writetable(T,filename,'Delimiter','\t','FileType','text','WriteVariableNames',false);


