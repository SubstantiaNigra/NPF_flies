%% concatenate cell arrays - ah, my RAM...


%%
fd2_length = length(FD2); % the number of flies
stim_probs = [0 0.05 0.1 0.15 0.2 0.3 0.45 0.6 0.75 1]; % probabilities used in the experiment
prob = 0.05;

%% find indices of the conditions with probability equal only to prob or 0
both_sides_p = [];
single_side_p = [];
for i = 1:fd2_length
    if ((FD2{1, i}.stimprob(1, 1) == prob) && (FD2{1, i}.stimprob(1, 2) == 0)) || ((FD2{1, i}.stimprob(1, 1) == 0) && (FD2{1, i}.stimprob(1, 2) == prob))
        single_side_p(end+1) = i;
    elseif ((FD2{1, i}.stimprob(1, 1) == prob) && (FD2{1, i}.stimprob(1, 2) == prob))
        both_sides_p(end+1) = i;
    end
end

%% get vector of occurence probabilities for single- and both-side stimulation

% for both-side stimulation concatenate all Xpos cells 
% into single cell array
xpos = {};
xpos_length = 0;
for i = both_sides_p
    xpos_length = length(FD2{1, i}.Xpos.Arena(1, :));
    for j = 1:xpos_length
        xpos(1, end+1) = FD2{1, i}.Xpos.Arena(1, j);
    end
end

% calculate probability of occurence in a given x position
bins = 50;
probs = zeros(2, bins); % first row contains both-side stimulation and
% second row contains single-side stimulation occurence probs
probs(2,:) = x_probability(xpos, bins);


% for single-side stimulation
xpos = {};
xpos_length = 0;
for i = single_side_p
    xpos_length = length(FD2{1, i}.Xpos.Arena(1, :));
    for j = 1:xpos_length
        x = FD2{1, i}.Xpos.Arena(1, j);
        % flip Xpos vector (x) such that stimulation side indices match
        if FD2{1, i}.stimprob(1, 1) ~= prob
            x = flip(x);
        end
        xpos(1, end+1) = x;
    end
end

% calculate probability of occurence in a given x position
probs(1,:) = x_probability(xpos, bins);

%% plot the data
x_ax = linspace(0,1,50);
waterfall(x_ax, [1 2], probs)

formatSpec = "Occupancy distribution of fly populations for %d";
str = sprintf(formatSpec, prob);
title(str)
xlabel("normalized x position")
ylabel("condition: 1 - one-sided, 2 - two-sided")
zlabel("probability of occurence")
% legend("normalized x position", "condition", "probability of occurence")
grid on