function s = vec2str(v, delim, brackets, fmt)

if ~exist('delim', 'var'), delim = ', '; end
if ~ischar(delim), delim = ', '; end

if ~exist('brackets', 'var'), brackets = {'[',']'}; end 
if ~iscell(brackets), brackets = {'[',']'}; end 

if ~exist('fmt', 'var'), fmt = '%.1f'; end

if isempty(v), s = brackets; return; end

if ~isnumeric(v) || ~isvector(v) || ~isreal(v) 
    error('The input must be a real numeric vector'); 
end

if nargin > 1
    if ~ischar(delim) || ~ischar(fmt)
        error('Deliminiter and precision format inputs must be of char type');
    end
    
    errrmes_brackets = 'Either have empty brackets input or a character (or cell) array length of 2';
    if isempty(brackets)
        brackets = {'',''};
    elseif length(brackets) == 2
        if ischar(brackets), brackets = {brackets(1), brackets(2)}; end
        if ~iscell(brackets), error(errrmes_brackets); end
    else, error(errrmes_brackets);
    end
end

s = sprintf([fmt, delim], v); 
s(end-length(delim)+1:end) = ''; 
s = sprintf('%s%s%s', brackets{1}, s, brackets{2}); 

end