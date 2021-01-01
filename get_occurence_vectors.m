%% initialize variables that should be created only once (as I must load data three times)
variants_num = length(variants); % how many different variants we have
x_pos = cell(variants_num, 1); % to collect all stimulation cells for all probabilities
stim_probs = [0.05 0.1 0.15 0.2 0.25 0.3 0.35];

%% find indices of the conditions in a given field
indices = cell(1); % to collect indices of variants
fd2_length = length(FD2); % the number of experiments (=cells) in a given field 

for i = 1:variants_num
    condition = variants(i, :);
    temp_idx = [];
    for j = 1:fd2_length
        fd_variant = FD2{1, j}.stimprob;
        fd_variant_flipped = flip(fd_variant); % flip the probability variants
        if ismember(fd_variant, condition, 'rows') | ismember(fd_variant_flipped, condition, 'rows')
            temp_idx(end+1) = j;
        end
    end
    indices{1, end+1} = temp_idx;
end
%% remove 1st cell from indices as it is a side effect of coding :-)
indices = indices(2:end);
%% get vector of occurence probabilities for single- and both-side stimulation
% 
% concatenate all Xpos cells into single cell array for a given variant
Xpos_length = 0;
for i = 1:variants_num
%     x_pos_temp = {}; % temporary cell array for a given variant x positions
    if not(isempty(indices(i)))
        var_idcs = cell2mat(indices(i));
        for j = var_idcs
            Xpos_length = length(FD2{1, j}.Xpos.Arena(1, :));        
            for k = 1:Xpos_length
%                 x_pos_temp{1, end+1} = FD2{1, j}.Xpos.Arena(1, k);    
                x_pos{i, end+1} = cell2mat(FD2{1, j}.Xpos.Arena(1, k));
            end
        end
%         x_pos(i, end+1) = x_pos_temp;    
    end
end

%% remove as many empty cells as possible from the cell array
[~,idc] = sort(cellfun(@isempty,x_pos),2);
s = size(x_pos);
[idr,~] = ndgrid(1:s(1),1:s(2));
x_pos = x_pos(sub2ind(s,idr,idc));

first_empty_in_row = [];
for i = 1:variants_num
    non_emptis = find(~cellfun('isempty', x_pos(i, :)));
    first_empty_in_row(end+1) = max(non_emptis);
end
idx_to_delete = max(first_empty_in_row);
