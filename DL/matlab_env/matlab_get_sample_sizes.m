SAMPLES_PATH = '..\..\Data\rna_seq200k\all_samples';
sample_names = dir(fullfile(SAMPLES_PATH,'*'));
samples_folders = setdiff({sample_names([sample_names.isdir]).name},{'.','..'});
%display({S.name});
% display(samples_folders);
for sample_idx = 1:numel(samples_folders)
    files_in_folder = dir(fullfile(SAMPLES_PATH,samples_folders{sample_idx},'*')); % improve by specifying the file extension.
    for file_idx = 1:length(files_in_folder)
        n_file_length = length(files_in_folder(file_idx).name);
        if n_file_length > 4 && strcmp(files_in_folder(file_idx).name(n_file_length-3:n_file_length), '.mat')
            if files_in_folder(file_idx).name == "M102.mat"
               continue 
            end
            sample_data_path = fullfile(SAMPLES_PATH, samples_folders{sample_idx}, files_in_folder(file_idx).name);
            matlab_structure = load(sample_data_path);
            Data = getfield(matlab_structure, char(fieldnames(matlab_structure)));
            display(files_in_folder(file_idx).name);
            disp('gene length');
            disp(size(Data.counts,1));
            disp('cells length');
            disp(size(Data.counts, 2));
        end
    end
   
          
       
    
end