function [file_info, file_recordings] = matchfile_info_recordings(file_info, file_recordings)
% this is just in case there is a mismatch, like a file being deleted or
% created with the other being corrupted or something 
info_ids = {file_info.identifier}; 
recordings_ids = {file_recordings.identifier}; 

[sorted_info, sorted_ind_info] = sort(info_ids);
[sorted_recordings, sorted_ind_recordings] = sort(recordings_ids);

if ~isequal(sorted_info,sorted_recordings) || ...
        length(info_ids) ~= length(unique(info_ids)) || ...
        length(recordings_ids) ~= length(unique(recordings_ids))
    error(['The identifiers of the 2 file list structs need to agree,' ...
        ' and cannot contain any duplicates within each struct array']); 
end

file_info = file_info(sorted_ind_info); 
file_recordings = file_recordings(sorted_ind_recordings); 
end