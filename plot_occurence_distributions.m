function [] = plot_occurence_distributions(probs, bins, variants_number, variant)
    x_ax = linspace(0,1,bins); % normalized x positions for bins = bins
    waterfall(x_ax, linspace(1, variants_number, variants_number), probs)

%     formatSpec = "Occupancy distribution of fly populations for %d";
%     str = sprintf(formatSpec, prob);
    if variant == 0
        str = "Occupancy distribution of fly populations for baiting conditions";
        title(str);
    else
        str = "Occupancy distribution of fly populations for non-baiting conditions";
        title(str);
    end
    xlabel("normalized x position")
    ylabel("conditions")
    zlabel("probability of occurence")
    % legend("normalized x position", "condition", "probability of occurence")
    grid on