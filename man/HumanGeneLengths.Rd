\name{HumanGeneLengths}
\alias{HumanGeneLengths}
\docType{data}

\title{
Human Gene Lengths
}

\description{
Lengths in number of nucleotides (\code{gene.length}) for 20,290 human genes
}

\usage{data(HumanGeneLengths)}
\format{
  A data frame with 20,290 observations on the following variable.
  \describe{
    \item{\code{gene.length}}{a numeric vector}
  }
}

%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}

\source{
Hubbard, T., D. Andrews, M. Caccamo, G. Cameron, Y. Chen, M. Clamp, L. Clarke, G. Coates, T. Cox, F. Cunningham, V. Curwen, T. Cutts, T. Down, R. Durbin, X. M. Fernandez-Suarez, J. Gilbert, M. Hammond, J. Herrero, H. Hotz, K. Howe, V. Iyer, K. Jekosch, A. Kahari, A. Kasprzyk, D. Keefe, S. Keenan, F. Kokocinsci, D. London, I. Longden, G. McVicker, C. Melsopp, P. Meidl, S. Potter, G. Proctor, M. Rae, D. Rios, M. Schuster, S. Searle, J. Severin, G. Slater, D. Smedley, J. Smith, W. Spooner, A. Stabenau, J. Stalker, R. Storey, S. Trevanion, A. Ureta-Vidal, J. Vogel, S. White, C. Woodwark, and E. Birne. 2005. Ensembl 2005. \emph{Nucleic Acids Research} 33: D447-D453.
}

\references{
\url{http://www.ncbi.nlm.nih.gov/pmc/articles/PMC540092/}

\url{http://www.ensembl.org/}
}

\examples{
data(HumanGeneLengths)
str(HumanGeneLengths)
histogram(~gene.length, HumanGeneLengths, subset=gene.length<15000)

if (require(manipulate)){
  ## Resampled mean, sd, SE
  genes <- function(n){
    nreps = 1000
    sample.mean <- numeric(nreps) # vector for sample means
    sample.se <- numeric(nreps)   # vector for sample standard errors
    sample.sd <- numeric(nreps)   # vector for sample standard deviations
    for (i in 1:nreps){
      random.sample <- sample(HumanGeneLengths$gene.length, size = n)
      sample.mean[i] <- mean(random.sample)
      sample.sd[i] <- sd(random.sample)
      sample.se[i] <- se(random.sample)
    }
    hist.bins <- hist(sample.mean, breaks = 30, plot = FALSE)
    hist(sample.mean,
      breaks = 30, xlim = c(1000, 4000),
      xlab = "Sample mean length (nucleotides)",
      col = "red", main = "")
    abline(v = mean(sample.mean), col = "blue", lwd = 2)
    text(x = 3200, y = 0.6 * max(hist.bins$counts), 
      pos = 4, cex = 1.25,
      paste("n = ", n, 
        "\nmean = ", round(mean(sample.mean), digits = 1), 
        "\nsd = ", round(mean(sample.sd), digits = 1), 
        "\nse = ", round(mean(sample.se), digits = 1), sep = ""))
}
manipulate(genes(nreps), nreps = slider(1, 500))
}
}

\keyword{datasets}
