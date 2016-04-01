function write_brainet(filename, nodesXYZ, nodesColor, nodesSize, nodesLabel)
% write_brainet
% write a BrainetViewer compatible .node file for visualization

n = size(nodesXYZ,1);
if isempty(nodesLabel)
    nodesLabel = repmat('-',size(nodesXYZ,1),1);
end

if length(nodesColor)~= n || length(nodesSize)~= n
    error('unconsistent number of nodes');
end

outfile = fopen(filename,'w');
for i=1:n
    fprintf(outfile,'%g\t%g\t%g\t%g\t%g\t%s\n',nodesXYZ(i,1),nodesXYZ(i,2),nodesXYZ(i,3),nodesColor(i),nodesSize(i),nodesLabel{i});
end

fclose(outfile);

