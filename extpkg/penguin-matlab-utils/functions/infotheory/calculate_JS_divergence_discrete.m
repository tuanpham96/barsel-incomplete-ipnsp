function d = calculate_JS_divergence_discrete(p, q, a)
% Calculate alpha-Jensen–Shannon divergence aJSD for discrete distributions
% aJSD = 1/2 * (aDKL(P || Q) + aDKL(Q || P))
% Here `log` refers to `log2`
% a = 1/2 is default JSD, a = 1 is 1/2 x Jeffreys divergence
% Usage: 
%   calculate_JS_divergence_discrete(p, q)      a = 1/2 default -> JSD 
%   calculate_JS_divergence_discrete(p, q, a)   a in [0, 1] range 

if ~exist('a', 'var'), a = 1/2; end 
if a < 0 || a > 1, error('The input "a" must be in the range of [0, 1]'); end 

p = p(:); q = q(:);

if length(p) ~= length(q), error('The 2 distribution vectors need to have the same number of elements'); end
if any(p == 0 & q == 0) && a == 1, error('Cannot handle aJSD when a = 1 and exists mutual p=0 and q=0'); end

p = p / sum(p);
q = q / sum(q);

d = p .* log2( p ./ ((1-a)*p + a*q) ) + ...
    q .* log2( q ./ ((1-a)*q + a*p) );
d = 0.5 * sum(d(~isnan(d)));
end