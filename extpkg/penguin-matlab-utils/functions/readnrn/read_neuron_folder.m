function res = read_neuron_folder(folder_name)
all_files = dir(folder_name); 
all_files = arrayfun(@(x) fullfile(x.folder, x.name), all_files, 'uni', 0);
res = read_neuron_file(all_files); 
end