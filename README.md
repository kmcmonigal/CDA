# CDA
Covariant discriminant analysis (CDA) is a method that calculates the patterns and associated timeseries which maximize the ratio variance(X)/variance(Y) of two datasets, X and Y. This can be applied to climate data to isolate mechanisms that cause differences between products, such as between two model runs. 

It is generally advised to use CDA on datasets which contain the same set of variability. To do this, the inverse EOF of a third dataset, Z, is projected onto both X and Y prior to calculating the CDA patterns. The script eofi.m will calculate the EOFs and inverse EOFs. Climate observations are often used as the Z field. 
