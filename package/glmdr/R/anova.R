
anova.glmdr <- function(object, ..., tolerance = .Machine$double.eps^0.75,
    test = c("LRT", "Rao")) {

    test <- match.arg(test)

    objectlist <- list(object, ...)

    if (length(objectlist) < 2)
        stop("must compare two or more models")
    # We know anova.lm and anova.glm do something with just one
    # but what they do is often statistical nonsense, just something
    # for naive users to make fools of themselves with

    if (! all(sapply(objectlist, function(x) inherits(x, "glmdr"))))
        stop("some arguments not of class \"glmdr\"")

    nmodels <- length(objectlist)

    # check that all models have the same response
    responses <- lapply(objectlist, function(o) o$om$y)
    sameresp <- sapply(responses, function(y) identical(y, responses[[1]]))
    if (! all(sameresp))
        stop("not all models have same response vector")

    # check that all models have the same family
    families <- sapply(objectlist, function(o) o$om$family$family)
    samefam <- sapply(families, function(y) identical(y, families[1]))
    if (! all(samefam))
        stop("not all models have same family")
    fam <- unique(families)

    # check that models are nested, that is,
    # column spaces of model matrices are nested vector spaces

    modelmatrices <- lapply(objectlist, function(o) o$om$x)
    qrlist <- lapply(modelmatrices, qr)
    for (i in 2:nmodels) {
        modmat1 <- modelmatrices[[i - 1]]
        qr2 <- qrlist[[i]]
        resid1 <- qr.resid(qr2, modmat1)
        norm.modmat1 <- apply(modmat1, 2, function(x) sqrt(sum(x^2)))
        norm.resid1 <- apply(resid1, 2, function(x) sqrt(sum(x^2)))
        if (any(norm.resid1 / norm.modmat1 > tolerance))
            stop("model matrices not nested\nmodel ",
                i - 1, " and model ", i)
    }

    # also check that differences of offset vectors are contained
    # in column space of larger model

    modeloffsets <- list()
    for (i in 1:nmodels) {
        offs <- objectlist[[i]]$om$offset
        if (is.null(offs))
            offs <- 0 * objectlist[[i]]$om$y
        modeloffsets[[i]] <- offs
    }

    for (i in 2:nmodels) {
        off1 <- modeloffsets[[i - 1]]
        off2 <- modeloffsets[[i]]
        offdiff <- off1 - off2
        if (all(offdiff == 0)) next
        qr2 <- qrlist[[i]]
        residoffdiff <- qr.resid(qr2, offdiff)
        norm.offdiff <- sqrt(sum(offdiff^2))
        norm.residoffdiff <- sqrt(sum(offdiffresid^2))
        if (norm.residoffdiff / norm.offdiff > tolerance)
            stop("model offsets not nested\nmodel ",
                i - 1, " and model ", i)
    }

    # if we get to here, models are nested
    # so can do the tests

    for (i in 2:nmodels) {
        linearity <- objectlist[[i - 1]]$linearity
        if (is.null(linearity)) {
            # MLE for null hypothesis is in original model
            # do conventional test
        } else if (all(linearity)) {
            # MLE for null hypothesis is completely degenerate
            # do trivial test, test statistic = 0, df = 0, p-value = 1
        } else {
            # MLE for null hypothesis is partially but not completely
            # degenerate
            #
            # must refit alternative hypothesis with same conditioning
            # as null

            modmat2 <- objectlist[[i]]$om$x
            resp2 <- objectlist[[i]]$om$y
            offs2 <- objectlist[[i]]$om$offset
            if (is.null(offs2)) {
                gout2 <- glm(resp2 ~ 0 + modmat2, family = fam,
                    subset = linearity)
            } else {
                gout2 <- glm(resp2 ~ 0 + modmat2, family = fam,
                    subset = linearity, offset = offs2)
            }
        }
    }
}

