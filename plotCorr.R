#' Function that computes the correlations between a vector and a each column
#' of a matrix, and plots the significant correlations.
#'
#' @param vect A numeric vector.
#' @param vars A numeric data frame, with number of rows equal to the length of
#' the vector vect.
#' @param type Correlation method to pass to stats::cor.test function.
#' Possible values : "pearson", "kendall", or "spearman".
#' @param seuil Significance cut-off value. Default is 0.05.
#' @param filename To save plot on disk, provide filename without the extension.
#' If NULL, the plot is only returned but not saved.
#' @param plottitle Plot title. If NULL, a default title will be used.
#' @param color.positive Colour of the points in the plot that correspond to
#' positive significant correlation coefficient. Default is "darkgrey".
#' @param color.negative Colour of the points in the plot that correspond to
#' negative significant correlation coefficient. Default is "black".
#' @return The correlation coefficients.
plotCorr = function(vect,
                    vars,
                    type=c("pearson", "kendall", "spearman"),
                    seuil=.05,
                    filename=NULL,
                    plottitle=NULL,
                    color.positive="darkgrey",
                    color.negative="black") {
    ### Specify correlation type
    type = match.arg(type)
    ### Utilitary function
    simpleCap = function(x) {
        s = strsplit(x, " ")[[1]]
        paste(toupper(substring(s, 1,1)), substring(s, 2),
              sep="", collapse=" ")
    }
    ### Correlation test
    cortest = apply(vars, 2, cor.test, y=vect, method=type)
    cortest_pvalues = sapply(cortest, function(a) glance(a)$p.value)
    cortest_values = sapply(cortest, function(a) glance(a)$estimate)
    ### Keep only significant variables
    vars_select = colnames(vars)[which(cortest_pvalues < seuil)]
    values_select = sort(abs(cortest_values[vars_select]))
    sign.values_select = sign(cortest_values[vars_select][names(values_select)])
    ### Prepare plot
    dframe = data.frame(factor(1:length(values_select)),
                        values_select,
                        factor(sign.values_select, levels=c(-1, 1)))
    names(dframe) = c("Vars", "Abs.cor", "Signe.cor")
    tickcols = sapply(rownames(dframe), function(a)
        ifelse(dframe[a, "Signe.cor"]==-1, color.negative, color.positive))
    ### Plot correlation
    if (is.null(plottitle)) {
        plottitle=paste0("Significantly correlated variables (",
                         length(values_select), ")")
    }
    gp = ggplot(dframe, aes(x=Abs.cor, y=Vars, color=Signe.cor)) +
        geom_point(size=4) +
        scale_color_manual(values=c(color.negative, color.positive),
                           guide=guide_legend(title="Correlation sign"),
                           drop=FALSE) +
        scale_y_discrete(breaks=1:nrow(dframe), labels=rownames(dframe)) +
        scale_x_continuous(name=paste("Absolute value of the correlation of",
                                      simpleCap(type)),
                           breaks=seq(floor(min(values_select)/0.05)*0.05,
                                      ceiling(max(values_select)/0.05)*0.05,
                                      by=0.05),
                           limits=c(floor(min(values_select)/0.05)*0.05,
                                    ceiling(max(values_select)/0.05)*0.05)) +
        theme(legend.position='bottom',
              axis.text.y=element_text(size=10, hjust=0.5, color=tickcols),
              axis.text.x=element_text(size=11),
              plot.title=element_text(size=18, hjust=0.5),
              axis.title.y=element_text(size=14),
              axis.title.x=element_text(size=14),
              legend.text=element_text(size=15),
              legend.title=element_text(size=18)) +
        ggtitle(plottitle)
    print(gp)
    ### Save to disk
    if (!is.null(filename)) {
        png(paste0(filename, ".png"))
        print(gp)
        dev.off()
    }
    
    return(cortest_values)
}
