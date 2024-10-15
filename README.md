# phylosignalDB

<!-- badges: start -->

<!-- badges: end -->

The goal of phylosignalDB is to provide a unified method, called *M* statistic, for detecting phylogenetic signals in continuous traits, discrete traits, and multi-trait combinations. Blomberg and Garland (2002) provided a widely accepted statistical definition of the phylogenetic signal, which is the "tendency for related species to resemble each other more than they resemble species drawn at random from the tree". The *M* statistic strictly adheres to the definition of phylogenetic signal, formulating an index and developing a method of testing in strict accordance with the definition, instead of relying on correlation analysis or evolutionary models. The novel method equivalently expressed the textual definition of the phylogenetic signal as an inequality equation of the phylogenetic and trait distances and constructed the *M* statistic. Also, there are more distance-based methods under development.

## Installation

You can install the development version of phylosignalDB from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("anonymous-eco/phylosignalDB")
```

## Example

We prepared an ecological trait dataset for turtles and used the *M* statistic in package `phylosignalDB` to detect the phylogenetic signals. The dataset was derived from the recently published ReptTraits dataset (Oskyrko et al., 2024). The phylogeny of turtles was derived by pruning from the maximum clade credibility tree provided in Thomson et al. (2021). Only those species that are present in both the ReptTraits dataset and the turtle phylogenetic tree were selected. Ultimately, the dataset comprised 240 species, encompassing 5 morphology traits, 2 behaviour traits, 2 life history traits, 5 habitat variables, and 2 variables concerning species conservation status.

We will now demonstrate how to use package `phylosignalDB` to detect the phylogenetic signals in continuous traits, discrete traits, and combinations of multiple traits.

``` r
library(phylosignalDB)

## load data
data("turtles")

# Continuous trait
trait_df <- data.frame(M1 = turtles$traits$M1, row.names = turtles$traits$specie)
trait_dist <- gower_dist(x = trait_df)
phylosignal_M(trait_dist, turtles$phylo, reps = 99) # reps=999 better

# Nominal discrete trait
trait_df <- data.frame(B1 = turtles$traits$B1, row.names = turtles$traits$specie)
trait_dist <- gower_dist(x = trait_df, type = list(factor = 1))
phylosignal_M(trait_dist, turtles$phylo, reps = 99) # reps=999 better

# Ordinal discrete trait
trait_df <- data.frame(CS1 = turtles$traits$CS1, row.names = turtles$traits$specie)
trait_dist <- gower_dist(x = trait_df, type = list(ordered = 1))
phylosignal_M(trait_dist, turtles$phylo, reps = 99) # reps=999 better

# Multi-trait Combinations
trait_df <- data.frame(turtles$traits[, c("M1", "M2", "M3", "M4", "M5")],
                       row.names = turtles$traits$specie)
trait_dist <- gower_dist(x = trait_df, type = list(factor = c("M4", "M5")))
phylosignal_M(trait_dist, turtles$phylo, reps = 99) # reps=999 better
```

## References

Blomberg, S.P. & Garland, T., Jr (2002) Tempo and mode in evolution: phylogenetic inertia, adaptation and comparative methods. Journal of Evolutionary Biology, 15(6): 899-910.

Oskyrko, O., Mi, C., Meiri, S. & Du, W. (2024) ReptTraits: a comprehensive dataset of ecological traits in reptiles. Scientific Data, 11(1): 243.

Thomson, R.C., Spinks, P.Q. & Shaffer, H.B. (2021) A global phylogeny of turtles reveals a burst of climate-associated diversification on continental margins. Proceedings of the National Academy of Sciences, 118(7): e2012215118.
