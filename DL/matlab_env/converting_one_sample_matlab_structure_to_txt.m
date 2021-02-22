folder = 'D:\Technion studies\Keren Laboratory\Data\droplet_seq\all_samples\M139';




T = dir(fullfile(folder,'*'));
C = {T(~[T.isdir]).name}; % files in subfolder.


for jj = 1:numel(C)
    file = fullfile(folder,C{jj});
    matlab_structure = load(file);
    Data = getfield(matlab_structure, char(fieldnames(matlab_structure)));  
    display("DATA processed well");
    display(Data)
    
    % GeneName starts.
    out_path = fullfile(folder,"GeneName.txt");
    fid = fopen(out_path,'wt');
    for row = 1:size(Data.GeneName,1)
       fprintf(fid,'%s ',char(Data.GeneName(row)));
       fprintf(fid,'\n');   
    end
    fclose(fid);
    % GeneName ends. 

    % ensID starts.
    out_path = fullfile(folder,"ensID.txt");
    fid = fopen(out_path,'wt');
    for row = 1:size(Data.ensID,1)
       fprintf(fid,'%s ',char(Data.ensID(row)));
       fprintf(fid,'\n');   
    end
    fclose(fid);
    % ensID ends. 

    % sample_name starts.
    out_path = fullfile(folder,"sample_name.txt");
    fid = fopen(out_path,'wt');
    for row = 1:size(Data.sample_name,1)
       fprintf(fid,'%s ',char(Data.sample_name(row)));
       fprintf(fid,'\n');   
    end
    fclose(fid);
    % sample_name ends.

    % cells starts.     
    out_path = fullfile(folder,"counts.txt");
    fid = fopen(out_path,'wt');
    for row = 1:size(Data.counts,1)
        for col = 1:size(Data.counts,2)
            fprintf(fid,'%d ',char(Data.counts(row, col)));
        end
       fprintf(fid,'\n');   
    end
    fclose(fid);
    % cells ends.
end