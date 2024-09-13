################################################################################
# Run 10, 100, and 1000 permutations of a phenotype, map the phenotype, and
# get the maximum LOD from each permutation. Do this 1000 times for each of the
# permutation levels and show that the variance of the significance threshold
# decreases with increasing numbers of permutations.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2024-09-13
################################################################################

library(tidyverse)
library(qtl2)

out_dir = 'C:/Users/c-dgatti/Documents/classes/JAX/qtl-mapping/episodes/data'

# Read in the iron data.
iron  = read_cross2(file = system.file("extdata", "iron.zip", package="qtl2") )
map   = insert_pseudomarkers(map = iron$gmap, step = 1)
probs = calc_genoprob(cross = iron, map = map, error_prob = 0.002)


n_sim  = 1000
n_perm = c(10, 100, 1000)

for(i in seq_along(n_perm)) {

  print(n_perm[i])

  results = matrix(0, nrow = n_perm[i], ncol = n_sim)

  for(j in 1:n_sim) {

    print(paste('   ', j))

    perm_pheno = replicate(n = n_perm[i], expr = sample(log(iron$pheno[,'liver'])))
    dimnames(perm_pheno) = list(rownames(iron$pheno), 1:n_perm[i])
    lod = scan1(genoprobs = probs, pheno = perm_pheno)
    results[,j] = apply(lod, 2, max)

  } # for(j)

  saveRDS(results, 
          file = file.path(out_dir, paste0('sim_perm', n_perm[i], '.rds')))

} # for(i)


# Read in the data and make a boxplot of the significance threshold 
# estimates.
results_files = dir(out_dir, pattern = '^sim_perm', full.names = TRUE)
results = lapply(results_files, readRDS)
results = sapply(results, function(z) {
                   apply(z, 2, quantile, probs = 0.95)
                 })
colnames(results) = gsub('^sim_perm|\\.rds$', '', basename(results_files))

results = data.frame(results) %>%
            pivot_longer(cols = everything(), 
                         names_to = 'n_perm', values_to = 'threshold') %>%
            mutate(n_perm = str_replace(n_perm, 'X', ''))

png(file.path(out_dir, '..', 'figures', 'permutation_simulations.png'),
    width = 800, height = 800, res = 128)
p = results %>%
  ggplot(aes(n_perm, threshold)) +
    geom_boxplot() +
    labs(title = 'Variance of Sig. Thr. Estimate',
         x = 'Number of permutations',
         y = 'Significane Threshold') +
    theme(text = element_text(size = 24))
print(p)
dev.off()

