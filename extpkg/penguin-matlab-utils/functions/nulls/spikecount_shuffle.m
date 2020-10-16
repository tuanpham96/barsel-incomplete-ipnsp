function X = spikecount_shuffle(X, shuf_method, varargin)

[num_neurons, num_bins] = size(X);

switch upper(shuf_method)
    case 'BIN', shuf_fun = @bin_shuffle;
    case {'COLUMN', 'COL'}, shuf_fun = @column_shuffle;
    case {'SWAP', 'ROW'}, shuf_fun = @row_shuffle;
    case 'SHIFT', shuf_fun = @shift_shuffle;
    otherwise, error('The shuffle method "%s" is not accepted', shuf_method);
end

X = shuf_fun(X, num_neurons, num_bins, varargin{:});

end


function X = bin_shuffle(X, num_neurons, num_bins, p_nrn)
if ~exist('p_nrn', 'var'), p_nrn = 1; end
for i = 1:num_neurons
    if rand > p_nrn, continue; end
    X(i,:) = X(i, randperm(num_bins));
end
end

function X = column_shuffle(X, ~, num_bins, p_bin)
if ~exist('p_bin', 'var'), p_bin = 1; end
num_shuf_bin = round(p_bin*num_bins);
shuf_bin_idx = randperm(num_bins, num_shuf_bin);
X(:,shuf_bin_idx) = X(:,shuf_bin_idx(randperm(num_shuf_bin)));
end

function X = row_shuffle(X, num_neurons, ~, p_nrn)
if ~exist('p_nrn', 'var'), p_nrn = 1; end
num_shuf_nrn = round(p_nrn*num_neurons);
shuf_nrn_idx = randperm(num_neurons, num_shuf_nrn);
X(shuf_nrn_idx,:) = X(shuf_nrn_idx(randperm(num_shuf_nrn)),:);
end

function X = shift_shuffle(X, num_neurons, num_bins, max_dbin)
if ~exist('max_dbin', 'var'), max_dbin = ceil(num_bins/2); end
dbins = randi([-1,1]*max_dbin, num_neurons, 1);
for i = 1:num_neurons
    X(i,:) = circshift(X(i,:), dbins(i));
end
end
