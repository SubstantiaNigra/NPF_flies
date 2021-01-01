%% load the data
% load('/home/daniel/Dokumenty/Internships/Data Analysis - Aarhus/Fly_data/NPF_flies/NPF_2.mat')

%% initialize variables that should be created only once (as I must load data three times
variants = zeros(1,2); % variant [0 0] already added
%% get all variants; drop variants with prob = 0

fd2_length = length(FD2); % the number of flies

% to get all baiting variants
for i = 1:fd2_length
    x = FD2{1, i}.stimprob;
    y = flip(x);
    if ((not(ismember(x, variants, 'rows')) == 1)) & (not(ismember(y, variants, 'rows')) == 1) & (not(ismember(0, x)) == 1)
        variants(end+1, :) = x;
    end
end
variants_num = length(variants); % how many different variants we have

% [0 0.05 0.1 0.15 0.2 0.3 0.45 0.6 0.75 1] - all probs in the experiment