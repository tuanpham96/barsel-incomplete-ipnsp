function d = calculate_KL_divergence_discrete(p, q, a)
% Calculate skewed alpha-Kullback–Leibler divergence of discrete distributions
% aDKL(P || Q) = DKL(P || (1-a)*P + a*Q) = E_p(log[p/((1-a)p + aq)])
% Here `log` refers to `log2`
% Usage: 
%   calculate_KL_divergence_discrete(p, q)      a = 1 default -> DKL 
%   calculate_KL_divergence_discrete(p, q, a)   a in [0, 1] range 

if ~exist('a', 'var'), a = 1; end 
if a < 0 || a > 1, error('The input "a" must be in the range of [0, 1]'); end 

p = p(:); q = q(:);

if length(p) ~= length(q), error('The 2 distribution vectors need to have the same number of elements'); end
if any(p(q == 0) ~= 0) && a == 1, error('Any where q = 0, p has to also be 0'); end

p = p / sum(p);
q = q / sum(q);

d = p .* log2( p ./ ((1-a)*p + a*q) );
d = sum(d(~isnan(d)));
end