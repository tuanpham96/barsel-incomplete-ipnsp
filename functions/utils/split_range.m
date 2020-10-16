function [num_section, sections] = split_range(low_bound, upp_bound, section_length, section_distance)
if ~exist('section_distance', 'var'), section_distance = 0; end

if strcmp(section_length, 'full')
    num_section = 1; 
    sections = [low_bound, upp_bound];
    return; 
end

num_section = ceil((upp_bound - low_bound)/section_length); 
sections = nan(num_section, 2); 

sec_stop = low_bound - section_distance; 
for i = 1:num_section 
    sec_start = sec_stop + section_distance; 
    sec_stop = min(sec_start + section_length, upp_bound);
    sections(i,:) = [sec_start, sec_stop];
end
    
end