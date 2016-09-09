# predict-tree-growth-with-a-hyperspectral-image
Code for a least squares algorithm to select narrowband indices to predict tree growth from a single hyperspectral image
This code is associated with the following paper in the journal Ecological Applications:

Title: A hyperspectral image can predict tropical tree growth rates in single-species stands

Authors: T. Trevor Caughlin, Sarah J. Graves, Gregory P. Asner, Michiel van Breugel, Jefferson S.
Hall, Roberta E. Martin, Mark S. Ashton, Stephanie A. Bohlman

T. Trevor Caughlin is the author of the R code. This code uses the PRESS statistic to select narrowband indices, the normalized difference between two narrowbands from hyperspectral data, to predict tree growth rates. The concept is not limited to tree growth rates and could be used for any continuous metric that needs to be predicted from a hyperspectral image.
 
We recommend using the R files in this order:
1. Use "Least Squares Algorithm.R" to explore whether our approach is suitable for your dataset.
2. If you plan to use the algorithm, you should run the Randomization test in "Randomization test for least squares algorithm.R" to determine how many narrowband indices should be added to your model
3. Once you have determined an appropriate number of narrowband indices to add using the randomization test, you should run the R script "LOOCV for least squares algorithm.R" to calculate out-of-sample metrics of model fit from leave-one-out cross-validation
