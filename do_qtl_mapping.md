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
library(tidyverse)
```

``` output
── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
✔ dplyr     1.1.4     ✔ readr     2.1.5
✔ forcats   1.0.0     ✔ stringr   1.5.1
✔ ggplot2   3.5.1     ✔ tibble    3.2.1
✔ lubridate 1.9.3     ✔ tidyr     1.3.1
✔ purrr     1.0.2     
── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()
ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
library(qtl2)
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

Make histogram of the "prop.bm.mn.ret" column and asses whether it is normally
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

:::::::::::::::::::::::

:::::::::::::::::::::::::::::::::::::




::::::::::::::::::::::::::::::::::::: challenge 

## Challenge 2: How many markers are there?

This will be a little trickier. You can look at the structure of the genetic or 
physical map in the Environment tab.

:::::::::::::::::::::::: solution

This is a task that is easiest to answer at the command line in the
Console. If you expand `gmap` in the Environment tab, you can see that it is a 
list. Each element in the list is a named, numeric vector. So we will get the
length of each list element and sum them.

```r
sum(sapply(gmap, length))
```

```output
6750
```

::::::::::::::::::::::::


:::::::::::::::::::::::::::::::::::::





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