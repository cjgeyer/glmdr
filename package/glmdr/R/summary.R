
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
            c("MLE exists in Barndorff-Nielsen completion",
            "it is conditional on components of the response",
            "corresponding to object$linearity == FALSE being",
            "conditioned on their observed values"),
            type = "lcm",
            summary = sout,
            linearity = object$linearity)
    } else {
        result <- list(overview =
            c("MLE exists in Barndorff-Nielsen completion",
            "it is completely degenerate",
            "the MLE says the response actually observed is the only",
            "possible value that could ever be observed"),
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

