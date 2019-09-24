function makeOutput(filename)
    % Load the output file.
    load([filename,'.mat']);
    % Basic parameters.
    JJ = size(uarray,1); 
    assert(JJ == 2*Lambda);
    assert(MC == size(uarray,2));
    Ntout = size(uarray,3);
    % Flatten the array of u into a single list.
    ulist = reshape(uarray,1,[]);
    dulist = reshape(duarray,1,[]);
    
    % Compute the skewness of each trajectory.
    skews = skewness(uarray,1,1);
    skews = reshape(skews,MC,Ntout);
    meanskew = skewness(uarray,1,[1,2]);
    meanskew = reshape(meanskew,1,[]);
    tt = [0: Ntout-1] * tfin/(Ntout-1);    
    skewarray = [tt; meanskew; skews];    
    skewlist = reshape(skewarray,1,[]);
    
    % Write output to a text file.
    fileID = fopen([filename,'.txt'],'w');
    % Write input parameters.
    fprintf(fileID,'#Basic integers: JJ, MC, Ntout, Nw.\n');
    out1 = [JJ,MC,Ntout,Nw]; fprintf(fileID,'%d\n', out1);
    fprintf(fileID,'#Basic floats: theta, tfin, runtime (mins).\n');
    out2 = [theta, tfin, runtime]; fprintf(fileID,'%9.5f\n', out2);
    fprintf(fileID,'# All other inputs: C2,C3,Drat,gibd,fi,dt.\n');
    out3 = [C2,C3,Drat,gibd,fi,dt]; fprintf(fileID,'%9.5f\n', out3);
    fprintf(fileID,'# Pad with zeros until index 20.\n');
    len = 20 - length(out1)-length(out2)-length(out3);
    fprintf(fileID,'%d\n', zeros(len,1) );
    % Write microstates u and du.
    fprintf(fileID,'# Values of u\n'); fprintf(fileID,'%9.5f\n', ulist);
    fprintf(fileID,'# Values of du\n'); fprintf(fileID,'%9.5f\n', dulist);
    % Write skewnesses.
    fprintf(fileID,'# skews\n'); fprintf(fileID,'%9.5f\n', skewlist);
    fclose(fileID);
end