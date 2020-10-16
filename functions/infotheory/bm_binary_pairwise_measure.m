num_words = 1000;
bit_length = 10; 
num_repeats = 1; 
X = logical(randi(2,num_words,bit_length)-1); 

measure_fun = @calculate_Jaccard_similarity_mex; 
tic
for i = 1:num_repeats
    distrib = get_binary_pairwise_measure_distribution(X, measure_fun);
end
toc