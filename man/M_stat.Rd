% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/M_statistic.R
\name{M_stat}
\alias{M_stat}
\title{Calculate M statistic}
\usage{
M_stat(trait_dist = NULL, phy = NULL, auto_multi2di = TRUE)
}
\arguments{
\item{trait_dist}{A distance object of class \code{matrix} or \link{dist}.
Its row and column names should match the tip labels of the phylogenetic tree (\code{phy}).
The functions \code{\link[=gower_dist]{gower_dist()}} and \code{\link[cluster:daisy]{cluster::daisy()}} can be used to calculate distances using trait data.}

\item{phy}{A phylogenetic tree of class \link{phylo}.}

\item{auto_multi2di}{A logical switch, \code{TRUE} or \code{FALSE}. Default is \code{TRUE},
then function \code{\link[=multi2di]{multi2di()}} in \code{ape} package will be called to make the phylogeney (tree)
be dichotomous if the tree (\code{phy}) contains some polytomies.}
}
\value{
A value that lies between 0 and 1, inclusive.
}
\description{
\code{M_stat} calculates the value of M statistic as a measurement of the strength of
the phylogenetic signal for the trait(s). The trait(s) could be continuous, discrete, or multi-variable.
Blomberg and Garland (2002) provided a widely accepted statistical definition of
the phylogenetic signal, which is the "tendency for related species to resemble
each other more than they resemble species drawn at random from the tree".
The M statistic strictly adheres to the definition of phylogenetic signal,
formulating an index and developing a method of testing in strict accordance
with the definition, instead of relying on correlation analysis or evolutionary models.
The novel method equivalently expressed the textual definition of the phylogenetic signal
as an inequality equation of the phylogenetic and trait distances and constructed the M statistic.
}
\examples{
data("turtles")
# Continuous trait
trait_df <- data.frame(M1 = turtles$traits$M1, row.names = turtles$traits$specie)
trait_dist <- gower_dist(x = trait_df)
M_stat(trait_dist, turtles$phylo)

# Nominal discrete trait
trait_df <- data.frame(B1 = turtles$traits$B1, row.names = turtles$traits$specie)
trait_dist <- gower_dist(x = trait_df, type = list(factor = 1))
M_stat(trait_dist, turtles$phylo)

# Ordinal discrete trait
trait_df <- data.frame(CS1 = turtles$traits$CS1, row.names = turtles$traits$specie)
trait_dist <- gower_dist(x = trait_df, type = list(ordered = 1))
M_stat(trait_dist, turtles$phylo)

# Multi-trait Combinations
trait_df <- data.frame(turtles$traits[, c("M1", "M2", "M3", "M4", "M5")],
                       row.names = turtles$traits$specie)
trait_dist <- gower_dist(x = trait_df, type = list(factor = c("M4", "M5")))
M_stat(trait_dist, turtles$phylo)

}
\references{
Blomberg, S.P. & Garland, T., Jr (2002) Tempo and mode in evolution: phylogenetic inertia, adaptation and comparative methods. Journal of Evolutionary Biology, 15(6): 899-910.
}
\seealso{
\code{\link[=M_rand_perm]{M_rand_perm()}} \code{\link[=phylosignal_M]{phylosignal_M()}}
}
