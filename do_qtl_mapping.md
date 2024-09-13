---
title: "QTL Mapping in Diversity Outbred Mice"
teaching: 120
exercises: 30
---

:::::::::::::::::::::::::::::::::::::: questions 

- How do I map traits in Diversity Outbred mice?
- How do I interpret the founder allele effects at a QTL peak?
- How do I perform association mapping in Diversity Outbred mice?
- How do I narrow down the set of candidate genes under a QTL?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Understand the data objects required for QTL mapping.
- Use the scan1 family of functions to perform different types of QTL scans.
- Interpret the founder allele effects.
- Find significant QTL peaks and identify candidate peaks under a QTL peak.
- Use auxiliary data to identify candidate genes under a QTL peak.

::::::::::::::::::::::::::::::::::::::::::::::::

## Introduction

This tutorial will take you through the process of mapping a QTL and searching 
for candidate genes for DNA damage from benzene exposure. The adverse health 
effects of benzene, including leukemia and aplastic anemia, were first studied 
in occupational settings in which workers were exposed to high benzene 
concentrations. Environmental sources of benzene exposure typically come from 
the petrochemical industry, however, a person’s total exposure can be increased 
from cigarettes, consumer products, gas stations, and gasoline powered engines 
or tools stored at home. 

Exposure thresholds for toxicants are often determined using animal models that
have limited genetic diversity, including standard inbred and outbred rats and 
mice. These animal models fail to capture the influence of genetic diversity on 
toxicity response, an important component of human responses to chemicals such 
as benzene. The Diversity Outbred (DO) mice reflect a level of genetic diversity
similar to that of humans.

The data comes from a toxicology study in which Diversity Outbred (DO) mice were
exposed to benzene via inhalation for 6 hours a day, 5 days a week for 4 weeks
[(French, J. E., et al. (2015) <i>Environ Health Perspect</i> 123(3): 237-245)](http://ehp.niehs.nih.gov/1408202/).
The study was conducted in two equally sized cohort of 300 male mice each, for a
total of 600 mice. They were then sacrificed and reticulocytes (red blood cell 
precursors) were isolated from bone marrow. The number of micro-nucleated 
reticulocytes, a measure of DNA damage, was then measured in each mouse. The 
goal is to map gene(s) that influence the level of DNA damage in the bone 
marrow. 

![Benzene Study Dosing](./figures/benzene_study_design.png){alt='Benzene study dosing showing 6 hours per day, 5 days per week of inhalation.'}

## DO Reference Data

As you work with DO data, you may need different reference files. These are
complied on a JAX website at: 
<https://www.jax.org/research-and-faculty/genetic-diversity-initiative/tools-data/diversity-outbred-reference-data>.
There are links to the markers and founder genotypes as well as the standard
colors that we use for the eight founder strains.

## Open the qtl_mapping Project.

If you have not opened the qtl_mapping project, do this now so that your
script and console will run in the correct directory. 

1. From the File menu, select "Open Project...".
1. Navigate to the project that you created in your Desktop/qtl_mapping folder
and open it.

This should reset your RStudio window and close any files.

## Create a New Markdown File

Create a new R Markdown file and save it in the "Desktop/scripts" folder as
"do_qtl_mapping.Rmd".

## Load Libraries

The primary library that we will use is [qtl2](https://github.com/kbroman/qtl2).
You should have `qtl2` installed from the Setup that you performed before the
workshop.


``` r
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggbeeswarm))
suppressPackageStartupMessages(library(qtl2))
suppressPackageStartupMessages(library(qtl2convert))
```

## Load and Explore the Data

The data for this tutorial has been saved as an R binary file that contains 
several data objects.  Load it in now by running the following command in your 
new script.

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: instructor

The code below is needed to load in the data objects to build the lesson on
Github. The `*.Rdata` file is over 200 Mb and Github has a file size limit of
100 Mb. The students should **NOT** run the next block.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: instructor

This is the only block that the students need to run.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


``` r
load("./data/qtl2_demo_grcm39.Rdata")
```

We loaded in four data objects from the *.Rdata file. Look at the Environment 
tab to see what was loaded. You should see an object called `pheno` which 
contains the phenotypes. There are two objects called `gmap` and `pmap`, which 
contain the marker positions in genetic and physical coordinates. Finally, there
is an object called `probs` which contains the founder allele probabilities for 
each mouse at each marker.

### Phenotypes

`pheno` is a data frame containing the phenotype data. Click on the blue circle
to the left of `pheno` in the Environment pane to view its contents.

::::::::::::::::::::::::::::::::::::: callout

The sample IDs **must** be in the rownames of `pheno` for qtl2 to work. 

:::::::::::::::::::::::::::::::::::::

`pheno` contains the sample ID, the study cohort, the concentration of 
benzene and several blood cell measurements. Note that the sample IDs are also 
stored in the rownames of `pheno`. 


In order to save time for this tutorial, we will only
map with 149 samples from the 100 ppm dosing group.

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 1: How many mice are there in the phenotype data?

Look in the Environment tab or use the Console to figure out how many mice
there are in the phenotype data.

:::::::::::::::::::::::: solution 

## Output
 
You can look at the number of observations in the Environment tab or you can
get the number of rows in `pheno` in the Console.
 

``` r
nrow(pheno)
```

``` output
[1] 598
```

There are 598 mice in the phenotype data.

:::::::::::::::::::::::

## Challenge 2: How many mice are there from each sex?

Use an R command to determine how many mice there are from each sex.

:::::::::::::::::::::::: solution 

## Output
 
You can look at the number of observations in the Environment tab or you can
get the number of rows in `pheno` in the Console.
 

``` r
count(pheno, sex)
```

``` output
   sex   n
1 male 598
```


All of the mice are male.

:::::::::::::::::::::::

## Challenge 3: How is the "prop.bm.mn.ret" phenotype distributed?

Make histogram of the "prop.bm.mn.ret" column and assess whether it is normally
distributed.

:::::::::::::::::::::::: solution 

## Output
 
You can look at the number of observations in the Environment tab or you can
get the number of rows in `pheno` in the Console.
 

``` r
ggplot(data = pheno, mapping = aes(prop.bm.mn.ret)) +
  geom_histogram() +
  labs(title = "Histogram of Bone Marrow MN-RETs") +
  theme(text = element_text(size = 20))
```

``` output
`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

``` warning
Warning: Removed 29 rows containing non-finite outside the scale range
(`stat_bin()`).
```

<img src="fig/do_qtl_mapping-rendered-unnamed-chunk-6-1.png" style="display: block; margin: auto;" />

The bone marrow MN-RET values do not look normally distributed.

:::::::::::::::::::::::

## Challenge 4: How is the log("prop.bm.mn.ret") phenotype distributed?

Take the log of the "prop.bm.mn.ret" column and make a histogram.Assess whether 
it is normally distributed.

:::::::::::::::::::::::: solution 

## Output
 
You can look at the number of observations in the Environment tab or you can
get the number of rows in `pheno` in the Console.
 

``` r
pheno |>
  mutate(log_mnret = log(prop.bm.mn.ret)) |>
  ggplot(mapping = aes(log_mnret)) +
    geom_histogram() +
  labs(title = "Histogram of log(Bone Marrow MN-RETs)") +
  theme(text = element_text(size = 20))
```

``` output
`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

``` warning
Warning: Removed 29 rows containing non-finite outside the scale range
(`stat_bin()`).
```

<img src="fig/do_qtl_mapping-rendered-unnamed-chunk-7-1.png" style="display: block; margin: auto;" />

The log of the bone marrow MN-RET values look more normally distributed.

:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::

Some researchers are concerned about the reproducibility of DO studies. The 
argument is that each DO mouse is unique and therefore can never be reproduced. 
But this misses the point of using the DO. While mice are the sampling unit, in 
QTL mapping we are sampling the founder alleles at each locus. An average of 
1/8th of the alleles should come from each founder at any given locus. Also, DO 
mice are a *population* of mice, not a single strain. While it is true that 
results in an individual DO mouse may not be reproducible, results at the 
population level should be reproducible. This is similar to the human population
in that results from one individual may not represent all humans, but results at
the population level should be reproducible.

The benzene inhalation study was conducted in two separate cohorts (termed 
*study* in the `pheno` object). Let's plot the proportion of micronucleated 
reticulocytes in bone marrow versus the benzene concentration for each study 
cohort.

:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: instructor

Remember, you just need to show the basic plot. Just type out the ggplot,
geom, log scale, and facet wrap lines in class.

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


``` r
pheno |>
  mutate(conc = factor(conc)) |>
  ggplot(mapping = aes(conc, prop.bm.mn.ret * 1000, color = conc)) +
    geom_violin(draw_quantiles = 0.5, linewidth = 1.25) +
    geom_beeswarm() +
    scale_y_log10() +
    scale_color_brewer(palette = "Reds") +
    facet_wrap(~study, ncol = 2) +
    labs(title = "Bone Marrow MN-RET by Study Cohort") +
    theme(text = element_text(size = 20), 
          legend.position = "none")
```

``` warning
Warning: Removed 29 rows containing non-finite outside the scale range
(`stat_ydensity()`).
```

``` warning
Warning: Removed 29 rows containing missing values or values outside the scale range
(`geom_point()`).
```

<img src="fig/do_qtl_mapping-rendered-unnamed-chunk-8-1.png" style="display: block; margin: auto;" />

As you can see, while individual mice have varying micronucleated 
reticulocytes, there is a dose-dependent increase in micronucleated reticulocytes
in both cohorts. This is an example of how results in the DO reproduce at the 
population level.

### Markers

The markers are the SNPs on the Mouse Universal Genotyping Array 
([MUGA](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4751547/)) 
that were used to reconstruct the haplotypes of the DO mice. This version of
the array had about 8,000 SNPs on it. The latest version contains over 
140,000 SNPs. We are using this smaller data set so that we can finish the
analysis in a timely manner during class.

The markers have been mapped to the latest mouse genome build (GRCm39) and
are provided by Karl Broman on Github at: 
<https://github.com/kbroman/MUGAarrays/tree/main/UWisc>. There are four
versions of the MUGA platforms. In this study, we used the MUGA and those
marker files are named as "muga_uwisc_vN.csv", where "N" is a version number.
We will use version 4, which you downloaded into your "data" directory 
during the Setup before the workshop.

Open the marker file now:


``` r
markers = read.csv("./data/muga_uwisc_v4.csv")
```


::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 5: Which columns contain the marker positions?

Look at the structure of the genetic or physical map in the Environment tab and
see which columns contain marker positions in some coordinate system.

:::::::::::::::::::::::: solution


``` r
str(markers)
```

``` output
'data.frame':	7854 obs. of  11 variables:
 $ marker        : chr  "JAX00240603" "UNC010001397" "UNC010515443" "UNC010001943" ...
 $ chr           : chr  "1" "1" "1" "1" ...
 $ bp_mm10       : int  3252796 3336839 3668628 3977130 4430623 4531029 5045931 5840130 6020779 6378154 ...
 $ bp_grcm39     : int  3323019 3407062 3738851 4047353 4500846 4601252 5116154 5910353 6091003 6448378 ...
 $ cM_cox        : num  0.15 0.172 0.194 0.197 0.202 ...
 $ strand        : chr  "minus" "minus" "minus" "plus" ...
 $ snp           : chr  "TC" "AC" "TC" "AC" ...
 $ unique        : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
 $ unmapped      : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
 $ probe         : chr  "AAAATATGATTCTTTTAATTAATGCATCAGTGATGAGAAGTGAGAGTGGC" "CATAGTGTCTGGTGAGAAGTCTGGAGTTATTCTAATAGGCCTGCCTTTAT" "CAGGAAATGATGCTGAGAAAGTGAGAAGTAGGAAAACGTGGAGAAAAATA" "TCTATTCCTATCACCTTGTACAAAGCTCAAGTCTTGTAAACCCCCCCCCC" ...
 $ strand_flipped: logi  TRUE FALSE FALSE FALSE FALSE FALSE ...
```

The columns labelled "bp_mm10", "bp_grcm39", and "cM_cox" contain marker
positions.

::::::::::::::::::::::::


## Challenge 6: Which chromosomes are in the markers?

Use an R command to determine which chromosomes are in `markers`.

:::::::::::::::::::::::: solution


``` r
count(markers, chr)
```

``` output
    chr   n
1     1 493
2    10 375
3    11 335
4    12 318
5    13 327
6    14 331
7    15 265
8    16 266
9    17 227
10   18 216
11   19 162
12    2 474
13    3 438
14    4 412
15    5 382
16    6 414
17    7 379
18    8 341
19    9 342
20    M   3
21    X 414
22    Y   5
23 <NA> 935
```

There are 23 chromosomes, including some markers with "NA". These are markers
that could not be uniquely mapped to one chromosome.

::::::::::::::::::::::::

## Challenge 6: What do you think the difference is between teh "bp_mm10" and
"bp_grcm39" columns?

Turn to your neighbor and discuss what the differences between the two column
are. Then share your ideas with the class.

:::::::::::::::::::::::: solution

There are different builds of the mouse genome. The Genome Resource Consortium 
(GRC) named mouse builds from build 1 to 39. The University of California at 
Santa Cruz developed their own naming system which conflicted with the GRC 
names. "GRCm38" was equivalent to "mm10". We are using the GRC naming 
conventions for genome builds.

::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::

[qtl2](https://github.com/kbroman/qtl2) needs the markers to be in a specific 
format. `qtl2` requires that the markers be in a list. The 
[qtl2convert](https://github.com/kbroman/qtl2convert) package provides functions
to make this easier.

First, we will filter the markers to retain ones on the autosomes and Chr X. We
don't use the mitochondrial or Y chromosomes in mapping, but they may be 
important covariates. The MUGA did not contain enough markers on Chr M and Y
to identify the contributing founder strain, but later versions of the MUGA, 
such as the GigaMUGA, do contain enough markers to identify founder 
contributions.


``` r
markers <- markers |>
             filter(chr %in% c(1:19, 'X'))
```

Next, we will create a column with GRCm39 Mb positions.


``` r
markers <- markers |>
             mutate(pos = bp_grcm39 / 1e6)
```

Finally, we will use the [qtl2convert::map_df_to_list()]()
function to create the list of markers that `qtl2` requires. We will call this
object the marker "map" to distinguish it from the marker data.frame.


``` r
map = qtl2convert::map_df_to_list(map = markers, pos_column = "pos")
```


### Founder Allele Probabilities (Genoprobs)

The `probs` object is a list with 20 elements. Each element of `probs` is a 
three-dimensional array containing the founder allele dosages for each sample at
each marker on one chromosome. These probabilities have been pre-calculated for 
you, so you can skip the step for calculating allele probabilities.


<!-- DMG: We need a figure showing the structure of the genoprobs. It should 
show a list of boxes. We can show the first three chromosomes and then ... -->

Let's look at the dimensions of `probs` for chromosome 1:


``` r
dim(probs[[1]])
```

``` output
[1] 590   8 479
```

`probs[[1]]` is a three-dimensional array containing the proportion of each 
founder haplotype at each marker for each DO sample.  The 590
samples are in the first dimension, the 8 founders in the second and the 
479 markers along chromosome 1 are in the third dimension.

Let's return to the `probs` object. Look at the contents of one sample on 
chromosome 1.


``` r
plot_genoprob(probs, map, ind = 1, chr = 1)
```

Uh oh! We have an error which says "Different numbers of positions in probs and 
map". This means that the number of markers in the genoprobs and the marker
map is different. In fact, this is quite important. The markers in `probs` and
`map` must match **exactly** on every chromosome.


::::::::::::::::::::::::::::::::::::: challenge

## Challenge 7: Align the markers in `probs` and `map`. 

Work with the person next to you to write a short script that will get the
intersection of the marker names in the `probs` and `map` for each 
chromosome, subset the markers in each object to only contain the common
markers, and make sure that the markers in `probs` are in the same order as in
`map`.

**Hint**: The markers names for each chromosome in `probs` are in dimention 3. 
i.e. dimnames(probs[[1]])[[3]] contains the marker names for chromosome 1.

:::::::::::::::::::::::: solution

We will write a loop that gets the common markers between probs and map for 
each chromosome and then subsets them.


``` r
# Verify that probs and map are the same length.
stopifnot(length(probs) == length(map))

# Verify that probs and map have their chromsomes in the same order.
stopifnot(names(probs) == names(map))

# Loop through each chromosome...
for(i in seq_along(probs)) {
  
  # Get the intersection of the markers in probs and map.
  common_markers <- intersect(dimnames(probs[[i]])[[3]], names(map[[i]]))
  
  # Subset the map to contain the common markers.
  map[[i]] <- map[[i]][common_markers]
  
  # Subset the probs to contain the common markers.
  probs[[i]] <- probs[[i]][,,names(map[[i]])]
  
  # Verify that the markers are identical in probs and map.
  stopifnot(dimnames(probs[[i]])[[3]] == names(map[[i]]))
  
} # for(i)
```

::::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::

Now that we have aligned the markers in `probs` and `map`, let's plot the 
allele probabilities for one sample on chromosome 1.


``` r
plot_genoprob(probs, map, ind = 1, chr = 1, main = "Founder Allele Probabilities")
```

<img src="fig/do_qtl_mapping-rendered-unnamed-chunk-18-1.png" style="display: block; margin: auto;" />

In the plot above, the founder contributions, which range between 0 and 1, are 
colored from white (= 0) to black (= 1.0). A value of ~0.5 is grey. The markers 
are on the X-axis and the eight founders (denoted by the letters A through H) on
the Y-axis. Starting at the left, we see that this sample has genotype BB 
because the row for B is black, indicating values of 1.0. Moving along the 
genome to the right, the genotype becomes CE where rows C and E are gray, 
followed by CD, FH, AG, GH, etc. The values at each marker sum to 1.0.  


### Calculating A Kinship Matrix

Next, we need to create a matrix that accounts for the kinship relationships 
between the mice. We do this by looking at the correlation between the founder 
haplotypes for each sample at each SNP. For each chromosome, we create a kinship
matrix using all markers *except* the ones on the current chromosome using the 
"loco" (leave-one-chromosome-out) method. Simulations suggest that mapping using
this approach increases the power to detect QTL.

::::::::::::::::::::::::::::::::::::: callout

The sample IDs must be in the rownames of `probs`. The sample IDs will be
copied to the row and column names in the kinship matrices.

:::::::::::::::::::::::::::::::::::::
           

``` r
K <- calc_kinship(probs = probs, type = "loco")
```

Kinship values between pairs of samples range between 0 (no relationship) and 
1.0 (completely identical). Let's look at the kinship matrix for the first
50 samples.


``` r
n_samples <- 50
heatmap(K[[1]][1:n_samples, 1:n_samples])
```

<img src="fig/do_qtl_mapping-rendered-kinship_probs-1.png" style="display: block; margin: auto;" />

The figure above shows kinship between all pairs of samples. Light yellow 
indicates low kinship and dark red indicates higher kinship. Orange values 
indicate varying levels of kinship between 0 and 1. The dark red diagonal of the
matrix indicates that each sample is identical to itself. The orange blocks 
along the diagonal may indicate close relatives (i.e. siblings or cousins).

<!-- DMG: STOPPED HERE -->



:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: instructor

Inline instructor notes can help inform instructors of timing challenges
associated with the lessons. They appear in the "Instructor View"

::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 1: Can you do it?

What is the output of this command?

```r
paste("This", "new", "lesson", "looks", "good")
```

:::::::::::::::::::::::: solution 

## Output
 
```output
[1] "This new lesson looks good"
```

:::::::::::::::::::::::::::::::::


## Challenge 2: how do you nest solutions within challenge blocks?

:::::::::::::::::::::::: solution 

You can add a line with at least three colons and a `solution` tag.

:::::::::::::::::::::::::::::::::
::::::::::::::::::::::::::::::::::::::::::::::::

## Figures

You can also include figures generated from R Markdown:


``` r
pie(
  c(Sky = 78, "Sunny side of pyramid" = 17, "Shady side of pyramid" = 5), 
  init.angle = 315, 
  col = c("deepskyblue", "yellow", "yellow3"), 
  border = FALSE
)
```

<div class="figure" style="text-align: center">
<img src="fig/do_qtl_mapping-rendered-pyramid-1.png" alt="pie chart illusion of a pyramid"  />
<p class="caption">Sun arise each and every morning</p>
</div>

Or you can use standard markdown for static figures with the following syntax:

`![optional caption that appears below the figure](figure url){alt='alt text for
accessibility purposes'}`

![You belong in The Carpentries!](https://raw.githubusercontent.com/carpentries/logo/master/Badge_Carpentries.svg){alt='Blue Carpentries hex person logo with no text.'}

::::::::::::::::::::::::::::::::::::: callout

Callout sections can highlight information.

They are sometimes used to emphasise particularly important points
but are also used in some lessons to present "asides": 
content that is not central to the narrative of the lesson,
e.g. by providing the answer to a commonly-asked question.

::::::::::::::::::::::::::::::::::::::::::::::::


## Math

One of our episodes contains $\LaTeX$ equations when describing how to create
dynamic reports with {knitr}, so we now use mathjax to describe this:

`$\alpha = \dfrac{1}{(1 - \beta)^2}$` becomes: $\alpha = \dfrac{1}{(1 - \beta)^2}$

Cool, right?


```r
sessionInfo()
```

::::::::::::::::::::::::::::::::::::: keypoints 

- Use `.md` files for episodes when you want static content
- Use `.Rmd` files for episodes when you need to generate output
- Run `sandpaper::check_lesson()` to identify any issues with your lesson
- Run `sandpaper::build_lesson()` to preview your lesson locally

::::::::::::::::::::::::::::::::::::::::::::::::

[r-markdown]: https://rmarkdown.rstudio.com/
