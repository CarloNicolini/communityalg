function write_brainet_community(filename, nodesXYZ, nodesColor, nodesSize, nodesLabel)
% write_brainet
% write a BrainetViewer compatible .node file for visualization but with
% every color (community) as a separate file
if isempty(nodesLabel)
    nodesLabel = repmat('-',size(nodesXYZ,1),1);
end

for c=unique(nodesColor(:))'
    filenamecomm = strcat(filename,'_c',num2str(c),'.node');
    nodes = find(nodesColor(:)==c);
    write_brainet(filenamecomm,nodesXYZ(nodes,:),nodesColor(nodes),nodesSize(nodes),nodesLabel(nodes,:));
end