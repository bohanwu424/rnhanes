# rnhanes
Associated code for paper	"A comparison of functional principal component analysis methods with accelerometry applications" (Wu, B and Van Allen, B.). Available at <a href="https://arxiv.org/abs/2105.14649" > arXiv:2105.14649</a>. 

The association between a person's physical activity and various health outcomes is an area of active research. The National Health and Nutrition Examination Survey (NHANES) data provide a valuable resource for studying these associations. NHANES accelerometry data has been used by many to measure individuals' activity levels. A common approach for analyzing accelerometry data is functional principal component analysis (FPCA). The first part of the paper uses Poisson FPCA (PFPCA), Gaussian FPCA (GFPCA), and nonnegative and regularized function decomposition (NARFD) to extract features from the count-valued NHANES accelerometry data. The second part of the paper compares logistic regression, random forests, and AdaBoost models based on GFPCA, NARFD, or PFPCA scores in the context of mortality prediction. The results show that Poisson FPCA is the best FPCA model for the inference of accelerometry data, and the AdaBoost model based on Poisson FPCA scores gives the best mortality prediction results.
