---
title: Setup
---

## Software Setup

R is a programming language that is especially powerful for data exploration, 
visualization, and statistical analysis. To interact with R, we use RStudio. 

1. Install the latest version of R from [CRAN](https://cran.r-project.org/).

2. Install the latest version of RStudio [here](https://www.rstudio.com/products/rstudio/download/). 
Choose the free RStudio Desktop version for Windows, Mac, or Linux. 

3. Start RStudio. 

4. Install packages. 
    a. The [qtl2](https://github.com/rqtl/qtl2) package contains code for
    haplotype reconstruction, QTL mapping and plotting. 
    b. The [qtl2convert](https://github.com/rqtl/qtl2convert) package contains
    code for converting data objects from one format to another.
    c. Install qtl2 by copying and pasting the following code in the R console.

```r
install.packages(c("qtl2", "qtl2convert"))
```

Once the installation is complete, load the libraries to make sure that they 
installed correctly. 

```r
library(qtl2)
library(qtl2convert)
```

If the libraries don't load and you recieved errors during the installation,
please contact the workshop instructors before the workshop to help you.

## Project organization

1. Create a new project in your Desktop called `qtl_mapping`. 
- Click the `File` menu button, then `New Project`.
- Click `New Directory`. 
- Click `New Project`.
- Type `qtl_mapping` as the directory name. Browse to your Desktop to create the project there.
- Click the `Create Project` button.

2. Use the `Files` tab to create  a `data` folder to hold the data, a `scripts` folder to 
house your scripts, and a `results` folder to hold results. Alternatively, you can use the 
R console to run the following commands for step 2 only. You still need to create a 
project with step 1.

```r
dir.create("./data")
dir.create("./scripts")
dir.create("./results")
```

## Data Sets

For this course, we will have several data files which you will need to 
download to the "data" directory in the project folder on your Desktop.
Copy, paste, and run the following code in the RStudio console.

The first file contains the Diversity Outbred mapping data.

```r
download.file(url      = "https://thejacksonlaboratory.box.com/shared/static/wspizp2jgrtngvvw5ixredpu7627mh5w.rdata",
              destfile = "data/qtl2_demo_grcm39.Rdata",
              mode     = "wb")
```

Next, we need a database of the DO founder SNPs and gene positions.

```r
download.file(url      = "https://figshare.com/ndownloader/files/40157572",
              destfile = "data/fv.2021.snps.db3",
              mode     = "wb")
```
