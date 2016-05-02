
glmdr <- function(formula, conditioner = NULL,
    family = c("poisson", "binomial"),
    data, subset, na.action, start = NULL,
    offset, control = list(...), model = TRUE,
    x = FALSE, y = TRUE, contrasts = NULL, ...)
{
    call <- match.call()
    family <- match.arg(family)

    stopifnot(inherits(formula, "formula"))
    stopifnot(is.null(conditioner) || inherits(conditioner, "formula"))

    # make data frame containing all variables
    if (missing(data)) {
        mf <- lm(big, method = "model.frame")
    } else {
        mf <- lm(big, data = data, method = "model.frame")
    }


    return(NULL)
}

