#' fpop
#' @description fpop is an extension of Fpop, an exact data segmentation algorithm based on functionnal pruning. This extension implements a penalty that depends on the size of the segments.
#' @param y a vector of data ordered according to an attribute
#' @param alpha a constant used in the calculation of the penalty (recommended: 10+2.5*log(length(y)))
#' @param wt a vector of weight linked to the data
#' @examples
#' path_to_Fpop <- path.package("Fpop")
#' df <- read.table(paste0(path_to_Fpop,"/wellLogData.txt"), header = TRUE)
#' n <- length(df$y)
#' beta <- 2.5
#' alpha <- 10 + beta * log(n)
#' res <- Fpop::fpop(df$y, beta, alpha)
#' res$changepoints
#' @return a list with the changepoints and the number of intervals/candidates at each step

fpop <- function(y, alpha,  muMinLocal=0, muMaxLocal=0, wt=-1, nbThreads=1)
{
        if (wt==-1)
        {
            wt = rep(1,length(y))
        }
        if (!muMinLocal & !muMaxLocal)
        {
            muMinLocal <- min(y)
            muMaxLocal <- max(y)
        }
        return (fpop_cpp(y, alpha, muMinLocal, muMaxLocal, wt, nbThreads))
}