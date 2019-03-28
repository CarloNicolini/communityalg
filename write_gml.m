function  write_gml( gmlfname, A, E, conf_true )
%WRITE_GML write network into gml file.

    file = fopen(gmlfname,'w');
    fprintf(file,'graph[\n');
    for i=1:length(A)
        fprintf(file,' node [\n');
        fprintf(file,'       id %d\n',i);
        fprintf(file,'       value %d\n',conf_true(i));
        fprintf(file,'      ]\n');
    end 
    
    for a=1:length(E)
        fprintf(file,' edge [\n');
        fprintf(file,'      source %d\n',E(a,1));
        fprintf(file,'      target %d\n',E(a,2));
        fprintf(file,'      ]\n');
    end 
    fprintf(file,']\n');
    fclose(file);
end

