function gen_bnv_pics(W, node_coords, node_labels, memb, prefix, folder, template_file, options_file )

write_brainet_community([folder '/' prefix], node_coords, memb , ones(length(W),1), node_labels);
um=unique(memb(:))';
for i=1:length(um)
    c=um(i);
    BrainNet_MapCfg(template_file, [folder '/' prefix '_c' num2str(c) '.node'], options_file, [folder '/' prefix '_c' num2str(c) '.png'] );
    close all
end