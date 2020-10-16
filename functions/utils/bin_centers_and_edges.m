function [centr, edges] = bin_centers_and_edges(v, bin_size)
rng_v = [min(v), max(v)] + [-1, 1] * bin_size; 
nbins = ceil(abs(diff(rng_v)) / bin_size); 
edges = linspace(rng_v(1), rng_v(2), nbins + 1); 
centr = (edges(1:end-1) + edges(2:end)) / 2; 
end