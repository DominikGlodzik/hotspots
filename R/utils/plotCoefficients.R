plotCoefficients <- function(nb.fit, nbfit.ref=NULL, main.text, model.variable.names, myxlim=c(0,1.8)) {

    coeff.df <- data.frame(coeff=nb.fit$coefficients, se=coef(summary(nb.fit))[, "Std. Error"])
    coeff.df <- coeff.df[2:nrow(coeff.df),]


    
    if (is.null(nbfit.ref)) {
        coeff.order <- order(coeff.df$coeff)
    } else {
        nbfit.ref$coefficients <- nbfit.ref$coefficients[2:length( nbfit.ref$coefficients)]
        coeff.order <- order(nbfit.ref$coefficients)
    }


    coeff.df<- coeff.df[coeff.order, ]
    coeff.df$labels <- model.variable.names[as.character(rownames(coeff.df))]
    bin.mean <- exp(nb.fit$coefficients[1])
    bin.theta <- summary(nb.fit)$theta
    bp <- barplot(exp(coeff.df$coeff), names.arg=coeff.df$labels,  border=NA, horiz=TRUE, main=paste(main.text, '\n mean ', signif(bin.mean,2), ', theta ', round(bin.theta ,2)), las=2, xlim=myxlim)
    abline(v=1)
    segments( exp(coeff.df$coeff-coeff.df$se), bp, exp(coeff.df$coeff+coeff.df$se), bp, lwd=2)

    coeff.df
}
