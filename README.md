TEMPTED Vignette
================

## Introduction of TEMPTED

<img src="Fig_tempted.png" style="width:100.0%" />

This is a vignette for the R package `tempted`, which implements the
statistical method TEMPoral TEnsor Decomposition (TEMPTED). The goal of
TEMPTED is to perform dimensionality reduction for multivariate
longitudinal data, with a special attention to longitudinal mirobiome
studies.

**Package dependencies:** R (\>= 4.2.0), np (\>= 0.60-17), ggplot2 (\>=
3.4.0), methods (\>= 4.2.1).

**Run time** for all the example codes in this demo was within one
minute on a MacBook Pro with 2.3 GHz Intel Core i7 processor. You may
expect longer run time if you have more subjects or more features.

You can **cite this paper** for using TEMPTED:

Shi P, Martino C, Han R, Janssen S, Buck G, Serrano M, Owzar K, Knight
R, Shenhav L, Zhang AR. [Time-Informed Dimensionality Reduction for
Longitudinal Microbiome Studies.
bioRxiv.](https://www.biorxiv.org/content/10.1101/2023.07.26.550749v1)

The statistical theories behind TEMPTED can be found in this paper:

Han R, Shi P, Zhang AR. [Guaranteed Functional Tensor Singular Value
Decomposition. Journal of the American Statistical Association (2023):
1-13.](https://www.tandfonline.com/doi/full/10.1080/01621459.2022.2153689)

TEMPTED is also **implemented in Python** through the Python package
`gemelli` and as a plugin of Qiime2. It can be installed through
`pip install gemelli`. [Documentation for
`gemelli`.](https://github.com/biocore/gemelli/tree/jointrpca)

## Installation

You can directly install TEMPTED from CRAN by

``` r
install.packages("tempted")
```

**Typical installation time** is within a few minutes on a typical
desktop.

You can install the development version of TEMPTED from
[GitHub](https://github.com/pixushi/tempted) with:

``` r
# install.packages("devtools")
devtools::install_github("pixushi/tempted")
```

or download the
[tarball](https://github.com/pixushi/tempted/blob/master/tempted_0.1.0.tar.gz)
to your folder of interest and install using

``` r
install.packages("folder_of_interest/tempted_0.1.0.tar.gz", repos = NULL, type="source")
```

## Load packages for this vignette

``` r
library(gridExtra)
library(tempted)
```

## Read the example data

The example dataset is originally from [Bokulich, Nicholas A., et
al. (2016)](https://pubmed.ncbi.nlm.nih.gov/27306664/). We provide three
data objects:

- `meta_table`: A data.frame with rows representing samples and matching
  with data.frame `count_table` and `processed_table` and three columns:

  - `studyid`: character denoting the subject ID of the infants.
  - `delivery`: character denoting the delivery mode of the infants.}
  - `day_of_life`: character denoting the age of infants measured in
    days when microbiome sample was taken.

- `count_tabe`: A data.frame with rows representing samples and matching
  with data.frame `meta_table` and columns representing microbial
  features (i.e. OTUs). Each entry is a read count.

- `processed_table`: A data.frame with rows representing samples and
  matching with data.frame `meta_table` and columns representing
  microbial features (i.e. ASVs, genes). Entries do not need to be
  transformed, and will be directly used by `tempted()`. This data.frame
  is used to illustrate how `tempted()` can be used for general form of
  multivariate longitudinal data already preprocessed by user.

``` r
# match rows of different data frames
# check number of samples
table(rownames(count_table)==rownames(meta_table))
#> 
#> TRUE 
#>  852
metauni <- unique(meta_table[,c('studyid', 'delivery')])
```

## Running TEMPTED for different formats of data

We provide example codes for the following scenarios:

1.  Run TEMPTED for Microbiome Count Data (Straightforward Way)
2.  Run TEMPTED for Microbiome Compositional Data (Straightforward Way)
3.  Run TEMPTED for General Form of Multivariate Longitudinal Data
    (Straightforward Way)
4.  Run TEMPTED in Customized Way
5.  Transferring TEMPTED result from training to testing data

## Run TEMPTED for Microbiome Count Data (Straightforward Way)

### Run TEMPTED

A complete description of all parameters can be found in the pdf manual.
Here we explain the key parameters:

- `feature_table`: A sample by feature matrix.
- `time_point`: The time stamp of each sample, matched with the rows of
  `feature_table`.
- `subjectID`: The subject ID of each sample, matched with the rows of
  `feature_table`.
- `threshold`: A threshold for feature filtering for microbiome data.
  Features with zero value percentage \>= threshold will be excluded.
  Default is 0.95.
- `smooth`: Smoothing parameter for RKHS norm. Larger means smoother
  temporal loading functions. Default is set to be 1e-8. Value can be
  adjusted depending on the dataset by checking the smoothness of the
  estimated temporal loading function in plot.
- `transform`: The transformation applied to the data. “logcomp” for log
  of compositions. “comp” for compositions. “ast” for arcsine squared
  transformation. “clr” for central log ratio transformation. “logit”
  for logit transformation. “lfb” for log fold over baseline. “none” for
  no transformation. Default transform=“clr” is recommended for
  microbiome data. For data that are already transformed, use
  transform=“none”.
- `r`: Number of components to decompose into, i.e. rank of the CP type
  decomposition. Default is set to 3.
- `pct_ratio`: The percent of features to sum up for taking log ratios.
  Default is 0.05, i.e. 5%.
- `absolute`: `absolute = TRUE` means features are ranked by the
  absolute value of feature loadings, and the top `pct_ratio` percent of
  features are picked. `absolute = FALSE` means features are ranked by
  the original value of feature loadings, and the top and bottom
  `pct_ratio` percent of features are picked. Then ratio is taken as the
  abundance of the features with positive loading over the abundance of
  the features with negative loading. Input for ratio_feature.
- `pct_aggregate`: The percent of features to aggregate, features ranked
  by absolute value of the feature loading of each component. Default is
  1, which means 100% of features are aggregated. Setting
  `pct_aggregate = 0.01` means top 1% of features is aggregated, where
  features are ranked in absolute value of feature loading of each
  component. Input for aggregate_feature.
- `pseudo`: A small number to add to all the counts before normalizing
  into proportions and log transformation. Default is 1/2 of the
  smallest non-zero value that is specific for each sample. This pseudo
  count is added for `transform=c("logcomp", "clr", "logit", "lfb")`.

**IMPORTANT NOTE:** In matrix singular value decomposition, the sign of
subject scores and feature loadings can be flipped together. Similarly,
you can flip the signs of any pair of subject loadings, feature
loadings, and temporal loadings.

``` r
res_count <- tempted_all(count_table,
                         meta_table$day_of_life,
                         meta_table$studyid,
                         threshold=0.95,
                         transform="clr",
                         pseudo=0.5,
                         r=2,
                         smooth=1e-5,
                         pct_ratio=0.1,
                         pct_aggregate=1)
#> Calculate the 1th Component
#> Convergence reached at dif=3.57400400849311e-05, iter=4
#> Calculate the 2th Component
#> Convergence reached at dif=4.745809269163e-05, iter=7
```

### Low-dimensional representation of subjects

Subject loadings are stored in variable `A_hat` of the `tempted_all()`
output.

``` r
A_data <- metauni
rownames(A_data) <- A_data$studyid
A_data <- cbind(res_count$A_hat[rownames(A_data),], A_data)

p_subj <- ggplot(data=A_data, aes(x=A_data[,1], y=A_data[,2], color=delivery)) + 
  geom_point() +
  labs(x='Component 1', y='Component 2', title='subject loading') 
print(p_subj)
```

<img src="man/figures/README-plot_sub_loading-1.png" width="100%" />

### Plot the temporal loadings

Temporal loadings are stored in variable `Phi_hat` of the
`tempted_all()` output. We provide an R function `plot_time_loading()`
to plot these curves. Option `r` lets you decide how many components to
plot.

``` r
p_time <- plot_time_loading(res_count, r=2) + 
  geom_line(size=1.5) + 
  labs(title='temporal loadings', x='days')
print(p_time)
```

<img src="man/figures/README-plot_time_loading-1.png" width="100%" />

### Plot the feature loadings

Feature loadings are stored in variable `B_hat` of the `tempted_all()`
output.

``` r
p_feature <- ggplot(as.data.frame(res_count$B_hat), aes(x=PC1, y=PC2)) + 
  geom_point() + 
  labs(x='Component 1', y='Component 2', title='feature loading')
print(p_feature)
```

<img src="man/figures/README-plot_feature_loading-1.png" width="100%" />

### Plot log ratio of top features

The log ratios are stored in variable `metafeature_ratio` of the
`tempted_all()` output.

``` r
group <- unique(meta_table[,c("studyid", "delivery")])
plot_metafeature(res_count$metafeature_ratio, group) + xlab("Days of Life")
#> Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   
```

<img src="man/figures/README-plot_log_ratio-1.png" width="100%" />

### Plot subject trajectories

The subject trajectores are stored in `metafeature_aggregate` of the
`tempted_all()` output.

``` r
group <- unique(meta_table[,c("studyid", "delivery")])
plot_metafeature(res_count$metafeature_aggregate, group) + xlab("Days of Life")
#> Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   
```

<img src="man/figures/README-plot_sub_traj-1.png" width="100%" />

### Low-dimensional representation of samples

This is an alternative way to visualize the subject trajectories.
Instead of plotting against the time points, you can visualize the
samples in low-dimensional spaces.

``` r
tab_feat_obs <- res_count$metafeature_aggregate
colnames(tab_feat_obs)[2] <- 'studyid'
tab_feat_obs <- merge(tab_feat_obs, metauni)

reshape_feat_obs <- reshape(tab_feat_obs, 
                            idvar=c("studyid","timepoint") , 
                            v.names=c("value"), timevar="PC",
                            direction="wide")
colnames(reshape_feat_obs) <- 
  sub(".*value[.]", "",  colnames(reshape_feat_obs))
colnames(reshape_feat_obs)
#> [1] "studyid"   "timepoint" "delivery"  "PC1"       "PC2"

p_aggfeat_scatter <- ggplot(data=reshape_feat_obs, aes(x=PC1, y=PC2, 
                             color=timepoint, shape=delivery)) +
  geom_point() + scale_color_gradient(low = "#2b83ba", high = "#d7191c") + 
  labs(x='Component 1', y='Component 2', color='Day')
p_aggfeat_scatter
```

<img src="man/figures/README-plot_aggfeat_scatter-1.png" width="100%" />

``` r

# subsetting timepoint between 3 months and 1 year
p_aggfeat_scatter2 <- ggplot(data=dplyr::filter(reshape_feat_obs, timepoint<365 & timepoint>30),
                             aes(x=PC1, y=PC2, 
                             color=delivery)) +
  geom_point() + 
  labs(x='Component 1', y='Component 2', color='Delivery Mode')
p_aggfeat_scatter2
```

<img src="man/figures/README-plot_aggfeat_scatter-2.png" width="100%" />

## Run TEMPTED for Microbiome Compositional Data (Straightforward Way)

### Run TEMPTED

**IMPORTANT NOTE**: Different form the count data, `pseudo = NULL` is
used so that 1/2 of the smallest non-zero value is added to each sample.
This pseudo count is added for
`transform=c("logcomp", "clr", "logit", "lfb")`.

``` r
proportion_table <- count_table/rowSums(count_table)
res_proportion <- tempted_all(proportion_table,
                              meta_table$day_of_life,
                              meta_table$studyid,
                              threshold=0.95,
                              transform="clr",
                              pseudo=NULL,
                              r=2,
                              smooth=1e-5)
#> Calculate the 1th Component
#> Convergence reached at dif=3.57400414642375e-05, iter=4
#> Calculate the 2th Component
#> Convergence reached at dif=4.7458092689688e-05, iter=7
```

### Low-dimensional representation of subjects

``` r
A_data <- metauni
rownames(A_data) <- A_data$studyid
A_data <- cbind(res_proportion$A_hat[rownames(A_data),], A_data)

p_subj <- ggplot(data=A_data, aes(x=A_data[,1], y=A_data[,2], color=delivery)) + 
  geom_point() + 
  labs(x='Component 1', y='Component 2', title='subject loading') 
print(p_subj)
```

<img src="man/figures/README-plot_sub_loading2-1.png" width="100%" />

### Plot the temporal loadings

``` r
p_time <- plot_time_loading(res_proportion, r=2) + 
  geom_line(size=1.5) + 
  labs(title='temporal loadings', x='days')
print(p_time)
```

<img src="man/figures/README-plot_time_loading2-1.png" width="100%" />

### Plot the feature loadings

``` r
pfeature <- ggplot(as.data.frame(res_proportion$B_hat), aes(x=PC1, y=PC2)) + 
  geom_point() + 
  labs(x='Component 1', y='Component 2', title='feature loading')
print(pfeature)
```

<img src="man/figures/README-plot_feature_loading2-1.png" width="100%" />

### Plot log ratio of top features

``` r
group <- unique(meta_table[,c("studyid", "delivery")])
plot_metafeature(res_proportion$metafeature_ratio, group) + xlab("Days of Life")
#> Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   
```

<img src="man/figures/README-plot_log_ratio2-1.png" width="100%" />

### Plot subject trajectories

``` r
group <- unique(meta_table[,c("studyid", "delivery")])
plot_metafeature(res_proportion$metafeature_aggregate, group) + xlab("Days of Life")
#> Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   
```

<img src="man/figures/README-plot_sub_traj2-1.png" width="100%" />

### Low-dimensional representation of samples

This is an alternative way to visualize the subject trajectories.
Instead of plotting against the time points, you can visualize the
samples in low-dimensional spaces.

``` r
tab_feat_obs <- res_proportion$metafeature_aggregate
colnames(tab_feat_obs)[2] <- 'studyid'
tab_feat_obs <- merge(tab_feat_obs, metauni)

reshape_feat_obs <- reshape(tab_feat_obs, 
                            idvar=c("studyid","timepoint") , 
                            v.names=c("value"), timevar="PC",
                            direction="wide")
colnames(reshape_feat_obs) <- 
  sub(".*value[.]", "",  colnames(reshape_feat_obs))
colnames(reshape_feat_obs)
#> [1] "studyid"   "timepoint" "delivery"  "PC1"       "PC2"

p_aggfeat_scatter <- ggplot(data=reshape_feat_obs, aes(x=PC1, y=PC2, 
                             color=timepoint, shape=delivery)) +
  geom_point() + scale_color_gradient(low = "#2b83ba", high = "#d7191c") + 
  labs(x='Component 1', y='Component 2', color='Day')
p_aggfeat_scatter
```

<img src="man/figures/README-plot_aggfeat_scatter2-1.png" width="100%" />

``` r

# subsetting timepoint between 3 months and 1 year
p_aggfeat_scatter2 <- ggplot(data=dplyr::filter(reshape_feat_obs, timepoint<365 & timepoint>30),
                             aes(x=PC1, y=PC2, 
                             color=delivery)) +
  geom_point() + 
  labs(x='Component 1', y='Component 2', color='Delivery Mode')
p_aggfeat_scatter2
```

<img src="man/figures/README-plot_aggfeat_scatter2-2.png" width="100%" />

## Run TEMPTED for General Form of Multivariate Longitudinal Data (Straightforward Way)

### Run TEMPTED

**IMPORTANT NOTE**: Different form the microbiome data, no features are
going to be filtered out by setting `threshold=1`, no transformation is
performed by setting `transform="none"`, and no log ratio is calculated
by setting `do_ratio=FALSE`.

``` r
res_processed <- tempted_all(processed_table,
                             meta_table$day_of_life,
                             meta_table$studyid,
                             threshold=1,
                             transform="none",
                             r=2,
                             smooth=1e-5,
                             do_ratio=FALSE)
#> Calculate the 1th Component
#> Convergence reached at dif=3.8157412154411e-05, iter=4
#> Calculate the 2th Component
#> Convergence reached at dif=2.44419512850139e-05, iter=7
```

### Low-dimensional representation of subjects

``` r
A_data <- metauni
rownames(A_data) <- A_data$studyid
A_data <- cbind(res_processed$A_hat[rownames(A_data),], A_data)

p_subj <- ggplot(data=A_data, aes(x=A_data[,1], y=A_data[,2], color=delivery)) + 
  geom_point() + 
  labs(x='Component 1', y='Component 2', title='subject loading') 
print(p_subj)
```

<img src="man/figures/README-plot_sub_loading3-1.png" width="100%" />

### Plot the temporal loadings

``` r
p_time <- plot_time_loading(res_processed, r=2) + 
  geom_line(size=1.5) + 
  labs(title='temporal loadings', x='days')
print(p_time)
```

<img src="man/figures/README-plot_time_loading3-1.png" width="100%" />

### Plot the feature loadings

``` r
pfeature <- ggplot(as.data.frame(res_processed$B_hat), aes(x=PC1, y=PC2)) + 
  geom_point() + 
  labs(x='Component 1', y='Component 2', title='feature loading')
print(pfeature)
```

<img src="man/figures/README-plot_feature_loading3-1.png" width="100%" />

### Plot subject trajectories

``` r
group <- unique(meta_table[,c("studyid", "delivery")])
plot_metafeature(res_processed$metafeature_aggregate, group) + xlab("Days of Life")
#> Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 |Multistart 1 of 1 /Multistart 1 of 1 |Multistart 1 of 1 |                   
```

<img src="man/figures/README-plot_sub_traj3-1.png" width="100%" />

### Low-dimensional representation of samples

``` r
tab_feat_obs <- res_processed$metafeature_aggregate
colnames(tab_feat_obs)[2] <- 'studyid'
tab_feat_obs <- merge(tab_feat_obs, metauni)

reshape_feat_obs <- reshape(tab_feat_obs, 
                            idvar=c("studyid","timepoint") , 
                            v.names=c("value"), timevar="PC",
                            direction="wide")
colnames(reshape_feat_obs) <- 
  sub(".*value[.]", "",  colnames(reshape_feat_obs))
colnames(reshape_feat_obs)
#> [1] "studyid"   "timepoint" "delivery"  "PC1"       "PC2"

p_aggfeat_scatter <- ggplot(data=reshape_feat_obs, aes(x=PC1, y=PC2, 
                             color=timepoint, shape=delivery)) +
  geom_point() + scale_color_gradient(low = "#2b83ba", high = "#d7191c") + 
  labs(x='Component 1', y='Component 2', color='Day')
p_aggfeat_scatter
```

<img src="man/figures/README-plot_aggfeat_scatter3-1.png" width="100%" />

``` r

# subsetting timepoint between 3 months and 1 year
p_aggfeat_scatter2 <- ggplot(data=dplyr::filter(reshape_feat_obs, timepoint<365 & timepoint>30),
                             aes(x=PC1, y=PC2, 
                             color=delivery)) +
  geom_point() + 
  labs(x='Component 1', y='Component 2', color='Delivery Mode')
p_aggfeat_scatter2
```

<img src="man/figures/README-plot_aggfeat_scatter3-2.png" width="100%" />

## Run TEMPTED in Customized Way

This is a breakdown of what `tempted_all()` does in each function it
wraps into.

### Format and Transform Data

The function `format_tempted()` transforms the sample-by-feature table
and the corresponding time points and subject IDs into the list format
that is accepted by `tempted()`. It filters features with a lot of
zeros, perform transformation of the features, and add pseudo count for
transformations that cannot handle zeros. When a subject has multiple
samples from the same time point, `format_tempted()` will only keep the
first sample in the data input.

**IMPORTANT NOTE**: For read count table as `feature_table`, set
`pseudo=0.5`. For compositional/proportion data as `feature_table`,
default setting is appropriate. For non-microbiome data, set
`transform="none"` and `threshold=1`.

``` r
# format the data frames into a list that can be used by TEMPTED
datlist <- format_tempted(count_table, meta_table$day_of_life, meta_table$studyid, threshold=0.95, pseudo=0.5, transform="clr")
length(datlist)
#> [1] 42
print(dim(datlist[[1]]))
#> [1] 796  29
```

### Run TEMPTED

The function `svd_centralize()` uses matrix SVD to fit a constant
trajectory (straight flat line) for all subject-feature pairs, and
remove it from the data. This step is optional. If it is not used, we
recommend the user to add 1 to the rank parameter `r` in `tempted()` and
in general the first component estimated by `tempted()` will reflect
this constant trajectory. In this example, we used `r=2`. Option
`smooth` allows the user to choose the smoothness of the estimated
temporal loading.

``` r
# centralize data using matrix SVD after log 
svd_tempted <- svd_centralize(datlist)
res_tempted <- tempted(svd_tempted$datlist, r = 2, resolution = 101, smooth=1e-5)
#> Calculate the 1th Component
#> Convergence reached at dif=3.531461408802e-05, iter=4
#> Calculate the 2th Component
#> Convergence reached at dif=4.37783267396658e-05, iter=7
# alternatively, you can replace the two lines above by the two lines below.
# the 2nd and 3rd component in the result below will be 
# similar to the 1st and 2nd component in the result above.
#svd_tempted <- NULL
#res_tempted <- tempted(datlist, r = 3, resolution = 101)
```

The output of `tempted()` can be used in the same way as the output of
`tempted_all()` to plot temporal loadings, subject loadings, and feature
loadings as in the previous sections of this document.

### Log ratio of top/bottom feature

The feature loadings can be used to rank features. The abundance ratio
of top ranking features over bottom ranking features corresponding to
each component can be a biologically meaningful marker. We provide an R
function `ratio_feature()` to calculate such the abundance ratio. The
abundance ratio is stored in the output variable `metafeature_ratio` in
the output. An TRUE/FALSE vector indicating whether the feature is in
the top ranking or bottom ranking is stored in variable `toppct` and
`bottompct`, respectively. Below are trajectories of the aggregated
features using observed data. Users can choose the percentage for the
cutoff of top/bottom ranking features through option `pct` (by default
`pct=0.05`). By default `absolute=TRUE` means the features are chosen if
they rank in the top `pct` percentile in the absolute value of the
feature loadings, and abundance ratio is taken between the features with
positive loadings over negative loadings. When `absolute=FALSE`, the
features are chose if they rank in the top `pct` percentile and has
positive loading or rank in the bottom `pct` percentile and has negative
loading. We also provide an option `contrast` allowing users to rank
features using linear combinations of the feature loadings from
different components.

``` r
datlist_raw <- format_tempted(count_table, meta_table$day_of_life, meta_table$studyid, 
                            threshold=0.95, transform='none')
contrast <- cbind(c(1/2,1), c(1/2,-1))
contrast
#>      [,1] [,2]
#> [1,]  0.5  0.5
#> [2,]  1.0 -1.0

ratio_feat <- ratio_feature(res_tempted, datlist_raw,
                        contrast=contrast, pct=0.1)
tab_feat_obs <- ratio_feat$metafeature_ratio
colnames(tab_feat_obs)[2] <- 'studyid'
tab_feat_obs <- merge(tab_feat_obs, metauni)
colnames(tab_feat_obs)
#> [1] "studyid"   "value"     "timepoint" "PC"        "delivery"
p_feat_obs <- ggplot(data=tab_feat_obs, 
                      aes(x=timepoint, y=value, group=studyid, color=delivery)) +
  geom_line() + facet_wrap(.~PC, scales="free", nrow=1) + xlab("Days of Life")

p_feat_obs_summary <- plot_metafeature(ratio_feat$metafeature_ratio, group, bws=30, nrow=1) + xlab("Days of Life")

grid.arrange(p_feat_obs, p_feat_obs_summary, nrow=2)
```

<img src="man/figures/README-plot_ratio_traj_breakdown-1.png" width="100%" />

### Subject trajectories

The feature loadings can be used as weights to aggregate features. The
aggregation can be done using the low-rank denoised data tensor, or the
original observed data tensor. We provide an R function
`aggregate_feature()` to perform the aggregation. The aggregated
features using low-rank denoised tensor is stored in variable
`metafeature_aggregate_est` and the aggregated features using observed
data is stored in variable `metafeature_aggregate`. Below are
trajectories of the aggregated features using observed data. Only the
features with absolute loading in the top pct percentile (by default set
to 100%, i.e. all features, through option `pct=1`) are used for the
aggregation. We also provide an option `contrast` allowing users to
aggregate features using linear combinations of the feature loadings
from different components.

``` r
## observed, by individual subject
contrast <- cbind(c(1/2,1), c(1/2,-1))
agg_feat <- aggregate_feature(res_tempted, svd_tempted, datlist, 
                              contrast=contrast, pct=1)
tab_feat_obs <- agg_feat$metafeature_aggregate
colnames(tab_feat_obs)[2] <- 'studyid'
tab_feat_obs <- merge(tab_feat_obs, metauni)
colnames(tab_feat_obs)
#> [1] "studyid"   "value"     "timepoint" "PC"        "delivery"
p_feat_obs <- ggplot(data=tab_feat_obs, 
                      aes(x=timepoint, y=value, group=studyid, color=delivery)) +
  geom_line() + facet_wrap(.~PC, scales="free", nrow=1) + xlab("Days of Life")

p_feat_obs_summary <- plot_metafeature(agg_feat$metafeature_aggregate, group, bws=30, nrow=1) + xlab("Days of Life")

grid.arrange(p_feat_obs, p_feat_obs_summary, nrow=2)
```

<img src="man/figures/README-plot_sub_traj_breakdown-1.png" width="100%" />

### Plot trajectory of top features

The feature loadings can also be used to investigate the driving
features behind each component. Here we focus on the 2nd component and
provide the trajectories of the top features using the estimated and
observed data tensor, respectively.

``` r
proportion_table <- count_table/rowSums(count_table)
feat_names <- c("OTU4447072", "OTU4467447")
# individual samples
tab_sel_feat <- cbind(proportion_table[,feat_names], meta_table)
tab_sel_feat <- reshape(tab_sel_feat, 
                        varying=list(feat_names), 
                        times=feat_names, 
                        timevar="feature", 
                        v.names="RA", 
                        ids=rownames(tab_sel_feat),
                        direction="long")
p_topfeat <- ggplot(data=tab_sel_feat) + 
  geom_point(aes(x=day_of_life, y=RA, color=delivery)) +
  facet_wrap(vars(feature)) + 
  labs(x="Day of Life", y="Relative Abundance") + 
  coord_trans(y="sqrt")
# summary plot
p_topfeat_summary <- plot_feature_summary(proportion_table[,feat_names], 
                     meta_table$day_of_life, 
                     meta_table$delivery, 
                     bws=30) + 
  labs(x="Day of Life", y="Relative Abundance")
grid.arrange(p_topfeat, p_topfeat_summary, nrow=2)
```

<img src="man/figures/README-plot_top_feature-1.png" width="100%" />

## Transferring TEMPTED result from training to testing data

### Split the example data into training and testing

Here we take thes subject with `studyid="2"` as the testing subject, and
the remaining subjects as training subjects.

``` r
id_test <- meta_table$studyid=="2"

count_train <- count_table[!id_test,]
meta_train <- meta_table[!id_test,]

count_test <- count_table[id_test,]
meta_test <- meta_table[id_test,]
```

### Run TEMPTED on training subjects

``` r
datlist_train <- format_tempted(count_train,
                                meta_train$day_of_life,
                                meta_train$studyid,
                                threshold=0.95,
                                pseudo=0.5,
                                transform="clr")

mean_svd_train <- svd_centralize(datlist_train, r=1)

res_tempted_train <- tempted(mean_svd_train$datlist,
r=2, smooth=1e-5)
#> Calculate the 1th Component
#> Convergence reached at dif=3.89952383939211e-05, iter=4
#> Calculate the 2th Component
#> Convergence reached at dif=4.45458808530784e-05, iter=7
```

### Obtain subject loading for testing subject

**IMPORTANT NOTE**: the testing samples should contain the features in
the training data.

``` r
# get the overlapping features between testing and training
count_test <- count_test[,rownames(datlist_train[[1]])[-1]]

# format testing data
datlist_test <- format_tempted(count_test,
                               meta_test$day_of_life,
                               meta_test$studyid,
                               threshold=1,
                               pseudo=0.5,
                               transform="clr")

# estimate the subject loading of the testing subject
sub_test <- est_test_subject(datlist_test, res_tempted_train, mean_svd_train)
sub_test
#>           PC1       PC2
#> 2 -0.07967744 0.1688863
```

Here we use logistic regression to illustrate how the subject loadings
of the testing data can be used.

``` r
# train logistic regression classifier on training subjects
metauni <- unique(meta_table[,c("studyid", "delivery")])
rownames(metauni) <- metauni$studyid
Atrain <- as.data.frame(res_tempted_train$A_hat)
Atrain$delivery <- metauni[rownames(Atrain),"delivery"]=="Cesarean"
glm_train <- glm(delivery ~ PC1+PC2,
                 data=Atrain, family=binomial(link="logit"))
summary(glm_train)
#> 
#> Call:
#> glm(formula = delivery ~ PC1 + PC2, family = binomial(link = "logit"), 
#>     data = Atrain)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)   
#> (Intercept)   -2.965      1.652  -1.795  0.07266 . 
#> PC1          -11.334      9.129  -1.242  0.21440   
#> PC2           16.865      5.529   3.050  0.00229 **
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> (Dispersion parameter for binomial family taken to be 1)
#> 
#>     Null deviance: 55.637  on 40  degrees of freedom
#> Residual deviance: 34.562  on 38  degrees of freedom
#> AIC: 40.562
#> 
#> Number of Fisher Scoring iterations: 6

# predict the label of testing subject "2", whose true label is "Cesarean"
predict(glm_train, newdata=as.data.frame(sub_test), type="response")
#>         2 
#> 0.6870014
```
