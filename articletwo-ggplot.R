plot.alpha <- function(diversity,
                       alphabreaks = seq(-1, 0, 0.2),
                       alphalabels = c("(-1.0, -0.8]", "(-0.8, -0.6]", "(-0.6, -0.4]", "(-0.4, -0.2]", "(-0.2, -0.0]")) {
    diversity <- diversity %>%
        mutate(estimate = ifelse(response == "HYPERTENSION", log(alpha.effect), alpha.effect),
               estimate_fac = cut(estimate,
                                  breaks = alphabreaks,
                                  labels = alphalabels,
                                  include.lowest = TRUE),
               Qstar = ifelse(alpha.p < 0.05, '*', ' '))

    plot.alpha.diversity <- ggplot(diversity, aes(x = reorder(Name, -deseqplotorder(Name)),
                                                  y = 1,
                                                  fill = estimate_fac)) +
        geom_bar(stat="identity", color="black") +
        coord_flip() +
        scale_fill_brewer(palette="Blues",
                          name = 'Regression coefficient \nin linear model \nfor Shannon index',
                          drop=FALSE,
                          direction = -1) +
        theme_classic(20) +
        geom_point(aes(Name, y=0.5, shape=Qstar), show.legend=FALSE, color='black', size=14) +
        scale_shape_manual(name="",
                           values=c('*'='*', ' '=' '),
                           labels=c("*"="significant at\nFDR 0.05", ' '=' '),
                           breaks=c("*", ' ')) +
        xlab("") +
        ylab("") +
        theme(legend.position="right",
              legend.title=element_text(size = 14),
              legend.direction="vertical",
              legend.text=element_text(size = floor(0.75*20)),
              legend.box.margin=margin(t = 0, unit = "cm"),
              legend.margin=margin(0,0,0,0,"pt"),
              legend.justification="left",
              axis.line.y = element_line(linetype="blank"),
              axis.line.x = element_line(linetype="blank"),
              legend.title.align = 0,
              legend.text.align = 0,
              axis.text.y = element_text(color="black"),
              axis.text.x = element_text(color="black"))+
        scale_y_continuous(breaks=c(0.5), labels=c(expression(alpha)))

    plot.alpha.diversity.gtable <- ggplot_gtable(ggplot_build(plot.alpha.diversity))
    plot.alpha.diversity.legend <- legend_extractor(plot.alpha.diversity.gtable)
    plot.alpha.diversity.yaxis <- plot.alpha.diversity.gtable$grobs[[3]]

    plot.alpha.diversity <- plot.alpha.diversity +
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              axis.text.y=element_blank())

    plot.gtable.alpha.diversity.stripped <- ggplot_gtable(ggplot_build(plot.alpha.diversity))
    return(list(plot = plot.gtable.alpha.diversity.stripped$grobs[[6]],
                yaxis = plot.alpha.diversity.yaxis,
                xaxis = plot.alpha.diversity.gtable$grobs[[7]],
                legend = plot.alpha.diversity.legend))
}

plot.beta <- function(diversity,
                      fig.textsize = 15,
                      ymax = 0.0010) {
    diversity <- diversity %>% mutate(Qstar = ifelse(beta.p < 0.05, "*", "-"))
    plot.beta.diversity <- ggplot(diversity,
                                  aes(x = reorder(Name, deseqplotorder(Name)), y = beta.R2)) +
        geom_bar(stat="identity", color = "black", fill = "gray96") +
        geom_point(aes(y = beta.R2 - 0.00007, shape=Qstar), show.legend=FALSE, color='black', size=14) +
        scale_shape_manual(values=c('*'='*', '-'='')) +
        coord_flip() +
        theme_classic(fig.textsize) +
        ylab("Coefficient of determination for\nbeta diversity") +
        scale_y_continuous(limits = c(0, ymax),
                           breaks = seq(0, ymax, 0.0005),
                           labels = function(x) sprintf("%.2f%%", x*100)) +
        theme(legend.position="right",
              legend.direction="vertical",
              legend.box.margin=margin(0,0,0,0,"pt"),
              legend.margin=margin(0,0,0,0,"pt"),
              legend.justification="left",
              legend.title.align = 0,
              legend.text.align = 0)

    plot.gtable.beta.diversity <- ggplot_gtable(ggplot_build(plot.beta.diversity))
    plot.xlab.beta.diversity <- xlabb_extractor(plot.gtable.beta.diversity)

    plot.beta.diversity <- plot.beta.diversity  +
        theme(axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x  =  element_blank(),
              axis.title.y = element_blank(),
              legend.position = "none")
    
    plot.gtable.beta.diversity.stripped <- ggplot_gtable(ggplot_build(plot.beta.diversity))

    return(list(plot = plot.gtable.beta.diversity.stripped$grobs[[6]],
                xaxis = plot.gtable.beta.diversity$grobs[[7]],
                xlab = plot.xlab.beta.diversity))
}

myforestplot <- function(data,
                         xlab = "",
                         ylab = "",
                         interceptone = FALSE,
                         subtitle = "",
                         my_y_scale = scale_y_continuous()) {
    ggplot(data = data,
           aes(x = term, y = estimate, ymin = conf.low, ymax = conf.high)) +
        geom_point() +
        { if (interceptone) geom_hline(aes(fill=fake_salt), yintercept=1, linetype=2) } +
        geom_errorbar(width=0.1)+ 
        ggplot2::theme_minimal()  +
        theme(plot.title=element_text(size=16),
              axis.title=element_text(size=10),
              strip.background = element_blank()) +
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.ticks.x = element_line(),
              axis.ticks.y = element_line(),
              axis.line = element_line(colour = "black")) +
        scale_x_discrete(labels=seq(3)) +
        my_y_scale + 
        xlab(xlab) +
        ylab(ylab) +
        labs(subtitle = subtitle)
}


gvlma.table.plot <- function(model) {
    { sink("/dev/null");
        gvlma.model <- gvlma::gvlma(model)
        gvlma.modek.summary <- summary(gvlma.model)
        ; sink(); }
    ret <- list("table" = gvlma.modek.summary,
                 "plot" = gvlma.model)
    return(ret)
}

plot.diversities <- function(diversities) {
    alpha.diversity.min <- plot.alpha(diversities$min)
    alpha.diversity.max <- plot.alpha(diversities$max)

    beta.diversity.min <- plot.beta(diversities$min)
    beta.diversity.max <- plot.beta(diversities$max)
    
    gs4 <- list(alpha.diversity.min$yaxis,
                alpha.diversity.min$plot,
                beta.diversity.min$plot,
                alpha.diversity.min$legend,
                alpha.diversity.min$xaxis,
                beta.diversity.min$xaxis,
                beta.diversity.min$xlab,
                text_grob("Age- and sex-adjusted model", size = 18),
                alpha.diversity.max$yaxis,
                alpha.diversity.max$plot,
                beta.diversity.max$plot,  
                alpha.diversity.max$xaxis,
                beta.diversity.max$xaxis,
                beta.diversity.max$xlab,
                text_grob("Multivariable-adjusted model", size = 18))
 
    g4 <- arrangeGrob(grobs = gs4,
                      layout_matrix = rbind(
                          c(8,    8,    8, NA,   NA, NA, NA),
                          c(NA, 1,    2,    NA, 3, NA, NA),
                          c(NA, 1,    2,    NA, 3, NA, 4),
                          c(NA, NA, 5,    NA,6, NA, NA),
                          c(NA, NA, NA, NA,7, NA, NA),
                          c(15,    15,    15, NA,   NA, NA, NA),
                          c(NA, 9,    10,  NA,  11, NA, NA),
                          c(NA, 9,    10,  NA,  11, NA, NA),
                          c(NA, NA, 12,   NA, 13, NA, NA),
                          c(NA, NA, NA, NA,14, NA, NA)),
                      widths=c(0.7, 4, 0.7, 0.1, 5, 0.7, 3),
                      heights=rep(c(1.5, 1, 6, 1, 2), 2))
    return(g4)
}

deseqheatmap <- function(df, legend.name = "Log Base Two\nFold Change") {
    ggplot(df,
           aes(x = reorder(Name, deseqplotorder(Name)),
               y = ordered(Feature, levels=rev(sort(unique(Feature)))),
               fill = log2FoldChange)) +
        geom_tile(aes(fill = 1)) +
        geom_tile(colour="white") +
        scale_fill_distiller(palette = "RdBu",
                             name = legend.name,
                             limits=c(-1, 1),
                             breaks = c(-1, -0.5, 0.5, 0, 1),
                             na.value="grey30") +
        coord_fixed(ratio=1) +
        xlab("") +
        ylab("") +
        theme_classic() +  
        theme(axis.text.x = element_text(angle=90, size=10, hjust=1.0, vjust=0.5),
              axis.text.y = element_text(size=10),
              legend.title = element_text(size=10),
              legend.position = "right")
}

standardize_geomtile <- function(p) {
    g <- ggplot_build(p)
    return(g)
}

deseqplotorder <- function(x) {
    y <- seq(length(x))
    y[x == "Systolic BP"] = 0
    y[x == "Diastolic BP"] = 1
    y[x == "Pulse pressure"] = 2
    y[x == "MAP"] = 3
    y[x == "Hypertension"] = 4
    y
}

saltplot <- function(dset, ymax = 2) {
    ggplot(data = dset, aes(x = NA., y = Lactobacillus.Bacteria)) +
        geom_point(aes(color = SEX), size = 0.05, position = "jitter") +
        geom_smooth(aes(group = SEX, color = SEX), method='lm', formula= y~x, colour = "black", size = 0.4) +
        scale_x_continuous(limits=c(20, 205), expand = c(0, 0)) +
        scale_y_log10(limits = c(100, 10000), expand = c(0, 0)) +
        scale_colour_manual(name = "SEX",
                            labels = c("Male", "Female"),
                            breaks=c("0", "1"),
                            values = c("blue", "red")) +
        xlab("24-hour urinary sodium [mmol]") +
        ylab("Normalized Lactobacillus abundance count") +
        guides(colour = guide_legend(override.aes = list(size=4, alpha = 0.4))) +
        theme_classic() +
        theme(plot.title = element_blank(),
              legend.title = element_blank(),
              legend.position = c(0.9,0.95))
}

myoutlier <- function(list) {
    list < quantile(list, 0.25) - IQR(list) * 1.5 |
        list > quantile(list, 0.75) + IQR(list) * 1.5 |
        (list > quantile(list, 0.40) & list < quantile(list, 0.60)) # buggy jitterdodge
}

saltboxplot <- function(pseq, dds) {
    df.mu <- assays(dds)[["mu"]] %>%
        as.data.frame %>%
        tibble::rownames_to_column("rowname") %>%
        gather(sampleid, value, -rowname) %>%
        spread(rowname, value) %>%
        rename_all(replace.brackets)

    df.countmu <- full_join(meta(pseq) %>%
                            tibble::rownames_to_column("sampleid"),
                            df.mu %>% select(sampleid, Lactobacillus.Bacteria),
                            by = "sampleid") %>%
        mutate(group = factor(dplyr::ntile(Lactobacillus.Bacteria, 4))) %>%
        group_by(group) %>%
        mutate(outlier = myoutlier(NA.)) %>%
        ungroup

    ggplot(data = df.countmu, aes(x = group, y = NA.)) +
        geom_boxplot(lwd=0.1, outlier.shape = NA, notch = FALSE, fill = "gray96") +
        geom_jitter(alpha = 0.05, size = 0.7, width = 0.3, shape = 20) +
        scale_x_discrete(labels=c("1" = "Q1",
                                  "2" = "Q2",
                                  "3" = "Q3",
                                  "4" = "Q4")) +
        scale_y_continuous(expand = c(0, 0), lim = c(0, 220)) +
        xlab("Lactobacillus abundance") +
        ylab("24-hour urinary sodium [mmol]") +
        theme_classic() +
        theme(strip.background = element_blank(),
              strip.placement = "outside")
}

pcoaplot <- function(pcoa) {
    pcoa.vectors <- pcoa$vectors %>% as.data.frame
    eig.fracs <- pcoa$values$Relative_eig

    pcoa.df <- merge(pcoa.vectors, meta(pseq.species), by=0, all = TRUE) %>%
        mutate_at(vars(starts_with("Axis.")), .funs = funs(myscale(.)))

    axis_labeller <- function(variable, value){
        axis.labels <- list("Axis.2" = sprintf("PCoA Axis 2 (%.1f%%)", 100*eig.fracs[2]),
                            "Axis.3" = sprintf("PCoA Axis 3 (%.1f%%)", 100*eig.fracs[3]))
    }
    
    pcoa.plot <- ggplot(data=pcoa.df %>% gather(key, Axis, Axis.2, Axis.3),
                        aes(x=Axis.1, y=Axis, color = HYPERTENSION)) +
        facet_wrap(~key, scales = "free_y", strip.position = "left", labeller = axis_labeller) +
        geom_jitter(shape = ".") +
        xlab(sprintf("PCoA Axis 1 (%.1f%%)", 100*eig.fracs[1])) +
        ylab(NULL) +
        scale_colour_manual(name = "State and treatment",
                            labels = c("Normotensive", "Hypertensive"),
                            breaks=c(0, 1),
                            values = c("blue", "red")) +   
        scale_x_continuous(breaks=seq(-5, 5, 1)) +
        scale_y_continuous(breaks=seq(-5, 5, 1)) +
        theme_classic() +
    guides(colour = guide_legend(override.aes = list(size=4, shape = 19))) +
        coord_cartesian(clip = 'off') +
        theme(text = element_text(size=10),
              legend.position = c(0.35, 0.95),
              legend.title = element_blank(),
              legend.text = element_text(size = 10),
              legend.key.size = unit(5, "mm"),
              legend.background = element_blank(),
              strip.background = element_blank(),
              strip.placement = "outside",
              aspect.ratio = 1,
              strip.text = element_text(size = 12),
              axis.title.x = element_text(size = 12))
}

pca.kde2d <- function(pcoa, htn = 0, n = 500) {
    pcoa.vectors <- pcoa$vectors %>% as.data.frame
    eig.fracs <- pcoa$values$Relative_eig

    pcoa.df <- merge(pcoa.vectors, meta(pseq.species), by=0, all = TRUE) %>%
        mutate_at(vars(starts_with("Axis.")), .funs = funs(myscale(.)))

    kde2d(pcoa.df$Axis.1, pcoa.df$Axis.2, n = n)
}
