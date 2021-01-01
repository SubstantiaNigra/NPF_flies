# NPF_flies
A repository for data analysis results from the experiments conducted in Duda Lab, Aarhus.

To recreate the data from 'return_codes' folder:
- occurence distributions:
1. get_all_variants.m
  
  1.1 it creates matrix [p1 p2] with different probabilities variants

2. get_occurence_vectors.m
  
  2.1 it returns x-positions for all conditions (n x m cell array for n conditions and m cells for each condition (can be empty))

3. get_occurence_probabilities
  
  3.1 it uses 'x_probability' function to get x-position probabilities (n x bins matrix)
  
  3.2 it uses 'plot_occurence_distribution' function to plot baiting (add 0 at the end) and non-baiting (add 1 at the end of the function) conditions separately
