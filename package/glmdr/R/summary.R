#' Summarizing Generalized Linear Model Fits
#' 
#' These functions are all \code{\link{methods}} for class \code{glmdr} or \code{summary.glmdr} objects.
#' @param object an object of class \code{"glmdr"}, usually, a result of a call to \code{\link{glmdr}}.
#' @param x an object of class \code{"summary.glmdr"}, usually, a result of a call to \code{summary.glmdr}.
#' @param correlation logical; if \code{TRUE}, the correlation matrix of the estimated parameters is returned and printed.
#' @param digits the number of significant digits to use when printing.
#' @param symbolic.cor logical. If \code{TRUE}, print the correlations in 
#' a symbolic form (see \code{\link{symnum}}) rather than as numbers.
#' @param signif.stars logical. If \code{TRUE}, \sQuote{significance stars} are printed for each coefficient.
#' @param \dots further arguments passed to or from other methods.
#' @return \code{summary.glmdr} returns an object of class \code{"summary.glmdr"}, a list with components
#'  \item{overview}{the note explaining whether the MLE is in the OM or LCM
#'  and whether the MLE model is completely degenerate.}
#'  \item{type}{one of \code{"original"}, \code{"lcm"}, or \code{"degenerate"}.
#'  A machine readable form of the overview.}
#'  \item{linearity}{The linearity vector (see Details).  Not present if \code{type} is \code{"original"}.}
#'  \item{summary}{the result of calling \code{\link[stats]{summary.glm}}
#'  on the fit done by \code{\link[stats]{glm}} on either the OM or the LCM
#'  as the case may be.  Not present if \code{type} is \code{"degenerate"}.}
#' @usage \method{summary}{glmdr}(object, correlation = FALSE, symbolic.cor = FALSE, \dots)
#' \method{print}{summary.glmdr}(x, digits = max(3, getOption("digits") - 3),
#' symbolic.cor = x$symbolic.cor,signif.stars = getOption("show.signif.stars"), \dots)
#' @details    These functions call \code{\link[stats]{summary.glm}} and \code{\link[stats]{print.summary.glm}} to do their work.
#' See the details for those functions.
#' 
#' There are three cases.
#' \itemize{
#' \item The MLE exists in the original model (OM), in which case
#' a note is printed saying this and otherwise the output is the
#' same as if one had called \code{\link[stats]{glm}} instead of
#' \code{\link{glmdr}}.  This is the case where \code{\link[stats]{glm}}
#' does the right thing. 
#' \item The MLE in the Barndorff-Nielsen completion is completely degenerate,
#' concentrated at the observed value of the response vector, in which case
#' a note is printed saying this and there is no other output since a
#' completely degenerate model has no identifiable parameters.
#' \item The MLE in the Barndorff-Nielsen completion is not completely
#' degenerate, the limiting conditional model (LCM) conditions on
#' the cases such that \code{object$linearity} is \code{FALSE}, in which case
#' a note is printed saying this and the output is the
#' same as if one had called \code{\link[stats]{glm}} with argument
#' \code{subset} that produces this conditioning.}
#' @seealso \code{\link{glmdr}, \link{glm}}, \code{\link{summary}}.
#' @export
#' @export summary.glmdr
#' @examples
#' ## For examples see example(glmdr)
summary.glmdr <- function(object, correlation = FALSE, symbolic.cor = FALSE,
    ...)
{
    call <- match.call()

    if (is.null(object$linearity)) {
        mycall <- call
        mycall$object <- object$om
        mycall[[1L]] <- quote(base::summary)
        sout <- eval(mycall, parent.frame())
        result <- list(overview =
            "MLE exists in the conventional sense in the original model",
            type = "original",
            summary = sout)
    } else if (any(object$linearity)) {
        mycall <- call
        mycall$object <- object$lcm
        mycall[[1L]] <- quote(base::summary)
        sout <- eval(mycall, parent.frame())
        result <- list(overview =
            cat(paste("MLE exists in Barndorff-Nielsen completion",
            "it is conditional on components of the response",
            "corresponding to object$linearity == FALSE being",
            "conditioned on their observed values", sep="\n")),
            type = "lcm",
            summary = sout,
            linearity = object$linearity)
    } else {
        result <- list(overview =
            cat(paste("MLE exists in Barndorff-Nielsen completion",
            "it is completely degenerate",
            "the MLE says the response actually observed is the only",
            "possible value that could ever be observed",sep="\n")),
            type = "degenerate",
            linearity = object$linearity)
    }
    class(result) <- "summary.glmdr"
    result
}

print.summary.glmdr <- function(x, digits = max(3, getOption("digits") - 3),
    symbolic.cor = x$symbolic.cor,
    signif.stars = getOption("show.signif.stars"), ...)
{
    call <- match.call()

    cat("\n")
    for (o in x$overview)
        cat(o, "\n")
    cat("\n")

    if (x$type != "degenerate") {

        if (x$type == "original") {
            cat("GLM summary for original model\n\n")
        } else {
            cat("GLM summary for limiting conditional model\n\n")
        }

        call$x <- x$summary
        call[[1L]] <- quote(base::print)
        eval(call, parent.frame())
    }

    invisible(x)
}

