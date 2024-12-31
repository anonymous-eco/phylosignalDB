#' Ecological Traits and Phylogeny of Turtles
#'
#' @description
#' An ecological trait dataset for turtles. The dataset was derived from
#' the recently published ReptTraits dataset (Oskyrko et al., 2024),
#' extracting species classified under the major group Testudines (comprising 361 species).
#' Only ecological traits with more than 50% of the species having trait records were retained.
#' The phylogeny of turtles was derived by pruning from the maximum clade credibility tree with
#' 288 tips provided in Thomson et al. (2021). Only those species that are present in both the
#' ReptTraits dataset and the turtle phylogenetic tree were selected. Ultimately, the dataset
#' comprised 240 species, encompassing 5 morphology traits, 2 behaviour traits, 2 life history
#' traits, 5 habitat variables, and 2 variables concerning species conservation status.
#' @usage
#' turtles
#' data("turtles")
#' @format
#' `turtles` is a list object with 3 components:
#' \describe{
#'   \item{traits}{The ecological traits of turtles as an object of class `data.frame`/`tibble`. }
#'   \item{phylo}{The phylogeny of turtles as an object of class `phylo`.}
#'   \item{traits_info}{The full names and id of ecological traits.}
#'   More details in Oskyrko et al. (2024) and Thomson et al. (2021).
#' }
#'
#' @references
#' Oskyrko, O., Mi, C., Meiri, S. & Du, W. (2024) ReptTraits: a comprehensive dataset of ecological traits in reptiles. Scientific Data, 11(1): 243.
#'
#' Thomson, R.C., Spinks, P.Q. & Shaffer, H.B. (2021) A global phylogeny of turtles reveals a burst of climate-associated diversification on continental margins. Proceedings of the National Academy of Sciences, 118(7): e2012215118.
#'
"turtles"
