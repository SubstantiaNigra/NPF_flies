%% load the data
% load('/home/daniel/Dokumenty/Internships/Data Analysis - Aarhus/Fly_data/NPF_flies/NPF_2.mat')
% [0 0.05 0.1 0.15 0.2 0.3 0.45 0.6 0.75 1] - all probs in the experiment
%% initialize variables that should be created only once (as I must load data three times)
variants_num = length(variants); % how many different variants we have
stim_probs = [0.05 0.1 0.15 0.2 0.25 0.3 0.35];

%% find indices of the conditions in a given field
indices = cell(1); % to collect indices of variants and info about baiting
fd_length = length(FD); % the number of experiments (=cells) in a given field 

for i = 1:variants_num
    condition = variants(i, :);
    temp_idx = [];
    for j = 1:fd_length
        fd_variant = FD{1, j}.stimprob;
        fd_variant_flipped = flip(fd_variant); % flip the probability variants
        if ismember(fd_variant, condition, 'rows') | ismember(fd_variant_flipped, condition, 'rows')
            temp_idx(end+1) = j;
        elseif ismember(fd_variant, condition, 'rows') | ismember(fd_variant_flipped, condition, 'rows')
            temp_idx(end+1) = j;
        end
    end
    indices{1, end+1} = temp_idx;
end
%% remove 1st cell from indices as it is a side effect of coding :-)
indices = indices(2:end);
%% take only these flies that have at least 30 trials on both sides
var_number = 4; % [0.3 0.15]
var_idx = indices{1, var_number}(2);
fd_variant = FD{1, var_idx}.stimprob;
to_be_flipped = ~ismember(fd_variant, condition, 'rows'); % 1 if vectors should be flipped

side = 2; % if = 2, take both sides under consideration
nmin = 20; % minimal number of trials on a given side to be valid !!!
baited = FD{1, var_idx}.baiting; % if it was a baited condition
for i = var_idx
    flies_nr = length(FD{1, i}.TrialN.Arena); % how many flies in a set
    valid_flies = nan(side, flies_nr); % append 1 if a fly performed > 30 trials, 0 otherwise
    for j = 1:flies_nr
        valid1 = (~isempty(FD{1, i}.TrialN.Arena{j,1}) && FD{1, i}.TrialN.Arena{j,1} > nmin); % check left side
        valid2 = (~isempty(FD{1, i}.TrialN.Arena{j,2}) && FD{1, i}.TrialN.Arena{j,2} > nmin); % check right side of variant
        valid_flies(1, j) = valid1;
        valid_flies(2, j) = valid2;
    end
    valid_flies = logical(valid_flies);
    % sum rewards for each fly and each side
    sum_of_rewards = nan(side, flies_nr);
    for k = 1:flies_nr
        sum_of_rewards(1, k) = sum(double(~isnan(FD{1, i}.TrigEvent.Arena{k,1})));
        sum_of_rewards(2, k) = sum(double(~isnan(FD{1, i}.TrigEvent.Arena{k,2})));
    end
    maxrew1 = floor(mean(sum_of_rewards(1, valid_flies(1, :)), 'omitnan'));
    maxrew2 = floor(mean(sum_of_rewards(2, valid_flies(2, :)), 'omitnan'));
    maxrew = [maxrew1 maxrew2]; % mean reward got on each side
end

%% logistic regression (choices vs. rewards) for a single fly - statistics for each side
nhist = 10; % look at 10 trials in the past

% variables for left-side statistics
br_1 = nan(flies_nr,2*nhist);
Br_1 = cell(flies_nr,100);
FitInf1_1 = cell(flies_nr, 1);
Brkernel_1 = cell(flies_nr, 1);
XX_1 = [];
Y_1 = [];
% variables for right-side statistics
br_2 = nan(flies_nr,2*nhist);
Br_2 = cell(flies_nr,100);
FitInf1_2 = cell(flies_nr, 1);
Brkernel_2 = cell(flies_nr, 1);
XX_2 = [];
Y_2 = [];

ivec = 1:flies_nr;

for i = ivec(valid_flies(1,:))
    clear stats1_1.s stats1_2.s
    % -------------------------------------------
    
    nmax_1 = FD{1, var_idx}.TrialN.Arena{i,1}; % num of trials on left side
    nmax_2 = FD{1, var_idx}.TrialN.Arena{i,2}; % num of trials on right side
    clear X y
    Choices_1 = FD{1, var_idx}.Returns.Arena{i,1}(1,1:nmax_1);
    Choices_2 = FD{1, var_idx}.Returns.Arena{i,2}(1,1:nmax_2);
    Rewards_1 = ~isnan(FD{1, var_idx}.TrigEvent.Arena{i,1}(1:nmax_1));
    Rewards_2 = ~isnan(FD{1, var_idx}.TrigEvent.Arena{i,2}(1:nmax_2));

    Xchoices_1 = zeros(length(Rewards_1), nhist-1);
    X_1 = zeros(length(Rewards_1), nhist);
    Xchoices_2 = zeros(length(Rewards_2), nhist-1);
    X_2 = zeros(length(Rewards_2), nhist);    
    
    for o = 2:length(Rewards_1)%-(nhist-1)
        if o-nhist+1 <= 0
            X_1(o,abs(o-nhist-1):end) = (Rewards_1(1:o));
            Xchoices_1(o,abs(o-nhist-1):end) = (Choices_1(1:o-1));
        else
            X_1(o,:) = (Rewards_1(o-nhist+1:o));
            Xchoices_1(o,:) = (Choices_1(o-nhist+1:o-1));
        end
    end
    
    for o = 2:length(Rewards_2)%-(nhist-1)
        if o-nhist+1 <= 0
            X_2(o,abs(o-nhist-1):end) = (Rewards_2(1:o));
            Xchoices_2(o,abs(o-nhist-1):end) = (Choices_2(1:o-1));
        else
            X_2(o,:) = (Rewards_2(o-nhist+1:o));
            Xchoices_2(o,:) = (Choices_2(o-nhist+1:o-1));
        end
    end


    y_1 = (Choices_1);
    X_1 = zscore(X_1);
    Xchoices_1 = zscore(Xchoices_1);
    Xz = [X_1 Xchoices_1];

    [br_1(i,:),~, stats1_1.s(i)] = glmfit([X_1 Xchoices_1],y_1.','binomial', 'link','logit');

    y_2 = (Choices_2);
    X_2 = zscore(X_2);
    Xchoices_2 = zscore(Xchoices_2);
    Xz = [X_2 Xchoices_2];

    [br_2(i,:),~, stats1_2.s(i)] = glmfit([X_2 Xchoices_2],y_2.','binomial', 'link','logit');
   
    
end

[idcs1, ~] = find(~isnan(br_1));
idcs1 = unique(idcs1);
idcs1 = idcs1(logical(ismember(idcs1, ivec(valid_flies(1, :))))); %ismember(idcs1, ivec(valid_flies(1, :)));
br_1 = rmmissing(br_1);
brp_1 = nan(size(br_1));

% statsvec_1 = 

for m = 1:size(br_1, 1)
    brp_1(m,:) = br_1(m,:).*(stats1_1.s(idcs1(m)).p<0.01).';
%     brp_1(m,:) = br_1(m,:).*(stats1_1.s(m).p<0.01).';
end

[idcs2, ~] = find(~isnan(br_2));
idcs2 = unique(idcs2);
idcs2 = idcs2(logical(ismember(idcs2, ivec(valid_flies(2, :)))));%ismember(idcs2, ivec(valid_flies(2, :)));
br_2 = rmmissing(br_2);
brp_2 = nan(size(br_2));
% statsvec_2 = idcs2(logical(ismember(idcs2, ivec(valid_flies(2, :)))));

for m = 1:size(br_2, 1)
%     brp_2(m,:) = br_2(m,:).*(stats1_2.s(statsvec_2(m)).p<0.01).';
    brp_2(m,:) = br_2(m,:).*(stats1_2.s(idcs1(m)).p<0.01).';
end

figure,
subplot(2,1,1),
shadedErrorBar(1:1:nhist, flip(nanmean(brp_1(:,2:nhist+1),1)), flip(nanstd(brp_1(:,2:nhist+1),1))./sqrt(size(br_1,1)), 'lineProps',{'r-o','markerfacecolor','r'}) 
hold on 
title([mat2str(variants(var_number, :)),' - left - log reg, baited: ', mat2str(baited)])
hold on  
shadedErrorBar(2:1:nhist, flip(nanmean(brp_1(:,nhist+2:end),1)), flip(nanstd(brp_1(:,nhist+2:end),1))./sqrt(size(br_1,1)), 'lineProps',{'g-o','markerfacecolor','g'})
legend('reward weights','choice weights','')
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 3;
ax.Box = 'on';
ax.FontSize = 24;
ax.XTickLabel = {'0','1','2','3','4','5','6','7','8','9'};
ylim([-0.05 0.5])


subplot(2,1,2),
shadedErrorBar(1:1:nhist, flip(nanmean(brp_2(:,2:nhist+1),1)), flip(nanstd(brp_2(:,2:nhist+1),1))./sqrt(size(br_2,1)), 'lineProps',{'r-o','markerfacecolor','r'}) 
hold on 
title([mat2str(variants(var_number, :)),' - right - log reg, baited: ', mat2str(baited)])
hold on  
shadedErrorBar(2:1:nhist, flip(nanmean(brp_2(:,nhist+2:end),1)), flip(nanstd(brp_2(:,nhist+2:end),1))./sqrt(size(br_2,1)), 'lineProps',{'g-o','markerfacecolor','g'})
legend('reward weights','choice weights','')
ax = gca;
ax.TickDir = 'out';
ax.LineWidth = 3;
ax.Box = 'on';
ax.FontSize = 24;
ax.XTickLabel = {'0','1','2','3','4','5','6','7','8','9'};
ylim([-0.05 0.5])
