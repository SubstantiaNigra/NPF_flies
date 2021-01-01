function [probs] = x_probability(x_pos, bins)
    samples = size(x_pos); % number of cells, takes empty cells for now too
    x_max = 50; % maximal x position to which a fly can walk
    empty_cells = 0; % to count empty cells
    for i = 1:samples(1,2)
        if isempty(x_pos{1, i}) == 0 % I add {1,1} at the end  
                                    % such that function works for
                                    % get_fields
            [occurences{i}, edges] = histcounts(x_pos{1, i}, bins); % count 
%                                                 occurences of given x position
            occurences{i} = occurences{i} / sum(occurences{i}); % normalize 
%                                                 with the number of time stamps
        else
            empty_cells = empty_cells + 1; % count empty cells such that
%                                           they are subtracted from number
%                                           of flies (samples)
            occurences{i} = zeros(1, bins); % if cell was empty, just fill bins with zeros
        end
    end

    probs = zeros(1, bins); % array of summed bins from each fly to obtain mean probs
    num_of_flies = samples(1,2) - empty_cells; % number of flies
    for i = 1:bins
        for j = 1:samples(1,2)
            probs(i) = probs(i) + occurences{1, j}(i);    
        end
        probs(i) = probs(i) / num_of_flies; % obtain mean for a given probability
    end
end
