%% calculate probability of occurence in a given x position
bins = 50;
variants_num = length(variants); % how many different variants we have 
probs = zeros(variants_num, bins); 
prob = "all conditions";

for i = 1:variants_num
    probs(i, :) = x_probability(x_pos(i,:), bins);
end

%%
% for i = 0:(variants_num-1)
%     % odd rows contain one-side stimulation
%     if ~((2*i+1) > conditions)
%         probs(2*i+1,:) = x_probability(xpos_single(i+1,:), bins);
%         % even rows contain two-side stimulation occurence probs
%         probs(2*i+2,:) = x_probability(xpos_both(i+1,:), bins);
%     end
% end

%% plot data baiting
baitings = [];
baiting_variants = []; % collect baiting-related data only
non_baitings = [];
non_baiting_variants = []; % collect non-baiting0related data only
for i = 1:variants_num
    if variants(i, 1) == variants(i, 2)
        non_baitings(end+1, :) = probs(i, :);
        non_baiting_variants(end+1, :) = variants(i, :);
    else
        baitings(end+1, :) = probs(i, :);
        baiting_variants(end+1, :) = variants(i, :);
    end
end
%% plot baiting and non-baiting conditions separately
plot_occurence_distributions(baitings, bins, length(baiting_variants), 0)
% plot_occurence_distributions(non_baitings, bins, length(non_baiting_variants), 1)
