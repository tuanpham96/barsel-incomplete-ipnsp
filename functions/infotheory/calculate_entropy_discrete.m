function ent = calculate_entropy_discrete(p)
% Calculate entropy of a discrete distribution (`log2`)
p = p(:); 
p = p / sum(p); 
ent = -p.*log2(p);
ent = sum(ent(~isnan(ent)));  
end