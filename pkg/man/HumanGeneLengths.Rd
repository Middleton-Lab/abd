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
  A data frame with 20290 observations on the following variable.
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

# Subset to only include genes with less than 15000 nucleotides
genes.under.15k <- subset(HumanGeneLengths, gene.length < 15000)

# Remove default space between the origin and the axes
par(xaxs = "i", yaxs = "i")

hist(genes.under.15k$gene.length,
  breaks = 30,
  ylim = c(0, 3000),
  xlab = "Gene length (number of nucleotides)",
  col = "red")

\dontrun{
require(ggplot2)
p <- ggplot(genes.under.15k, aes(gene.length))
p + geom_histogram(fill = "red") +
  scale_x_continuous("Gene length (number of nucleotides)") +
  scale_y_continuous("Frequency")
}


# Population mean and standard deviation
mean(HumanGeneLengths$gene.length)
sd(HumanGeneLengths$gene.length)


# Random sample of 100 gene lengths
set.seed(1234) # For repeatability
random.gene.lenghts <- sample(HumanGeneLengths$gene.length, size = 100)
# Note that random.gene.lengths is a vector, rather than a data.frame

par(xaxs = "i", yaxs = "i")
hist(random.gene.lenghts,
  breaks = 30,
  xlab = "Gene length (number of nucleotides)",
  col = "red")

# Sample mean and standard deviation
mean(random.gene.lenghts)
sd(random.gene.lenghts)


# Sampling distribution of mean gene length
set.seed(4321)
nreps <- 1000    # Number of iterations
sample.mean <- numeric(nreps)   # initialize a vector
                                # to hold the sample means

for (i in 1:nreps){
  random.sample <- sample(HumanGeneLengths$gene.length, size = 100)
  sample.mean[i] <- mean(random.sample)
  }

par(xaxs = "i", yaxs = "i")
hist(sample.mean,
  breaks = 30,
  xlab = "Sample mean length (nucleotides)",
  col = "red")

# Comparison of the distribution of sample means and standard errors
# for different sample sizes
set.seed(6)
par(xaxs = "i", yaxs = "i")
par(mfrow = c(3, 1))
nreps = 10000
for (n in c(20, 100, 500)){
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

  par(xaxs = "i", yaxs = "i")
  hist(sample.mean,
    breaks = 30,
    xlim = c(1000, 5000),
    xlab = "Sample mean length (nucleotides)",
    col = "red",
    main = "")
  abline(v = mean(sample.mean), col = "blue", lwd = 2)
  text(x = 3000, y = 0.7 * max(hist.bins$counts), 
    pos = 4,
    paste("n = ", n, 
      "\nmean = ", round(mean(sample.mean), digits = 2), 
      "\nsd = ", round(mean(sample.sd), digits = 2), 
      "\nse = ", round(mean(sample.se), digits = 2), sep = ""))
}

}

\keyword{datasets}
