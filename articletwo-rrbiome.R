# Helper functions for articletwo

`%difference%` <- function(a, b) {
    a[!a %in% b]
}

`%intersect%` <- function(a, b) {
    intersect(a, b)
}


`%union%` <- function(a, b) {
    c(a, b)
}

replace.brackets <- function (genus) gsub('(.*) \\(+(.+)\\)', '\\1.\\2', genus)

renametaxa <- function (names, star = "*") {
    suffix <- ifelse(grepl("BacteriaPlasmid", names), star, "")
    base <- case_when(grepl("_.*\\.Bacteria", names) ~ gsub('(.).*_(.*).Bacteria.*', '\\1. \\2', names),
                      grepl("\\.Bacteria", names) ~ gsub('(.*)\\.Bacteria.*', '\\1', names),
                      grepl("_", names) ~ gsub('(.).*_(.*) \\(+(.+)\\)', '\\1. \\2', names),
                      TRUE ~  gsub('(.*) \\(+(.+)\\)', '\\1', names))
    paste0(base, suffix)
}

mycbind <- function(list) {
    max.length <- c(list) %>% lapply(length) %>% unlist %>% max
    lapply(list, function(x) {
        length(x) <- max.length
        x}) %>%
        do.call(cbind, .)
}

list.partition <- function(list, parts = 3) {
    partition <- floor(parts*seq(0, length(list)-1)/length(list))
    split(list, partition)
}

rename.genus  <- function (genus, markword = "Plasmid", mark = "*") {
    name <- gsub('(.*) \\(.*\\)', '\\1', genus)
    star <- ifelse(grepl(markword, genus), mark, "")
    paste0(gsub("_", " ", name), star)
}

myscale <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}



meta.merge.alphadiversity <- function(pseq, index = "shannon") {
    alphadiversity  <- microbiome::alpha(pseq, index = index)
    base::merge(meta(pseq), alphadiversity, by=0, all=TRUE) %>%
        dplyr::rename(Sample_ID = Row.names) %>%
        dplyr::mutate(diversity_shannon = scale(diversity_shannon),
                      MAP = oMAP,
                      PULSEPRESSURE = oPULSEPRESSURE,
                      SYSTM = oSYSTM,
                      DIASM = oDIASM)
}

calculate.beta.matrix <- function(pseq) {
    compositional.profile <- microbiome::transform(pseq, 'compositional')
    otu <- microbiome::abundances(compositional.profile)
    bray.dist.m <- vegan::vegdist(t(otu), method="bray")
    dist.matrix <- as.matrix(bray.dist.m)
    attr(dist.matrix, "method") <- "bray"
    dist.matrix
}
calculateglm <- function(dset,
                         responses = list(model_1 = "MAP", model_2 = "SYSTM", model_3 = "DIASM",
                                          model_4 = "PULSEPRESSURE", model_5 = "HYPERTENSION"),
                         min_n_for_continuous = 10,
                         covariates = c(),
                         filterstr = ".") {
    glmlist <- lapply(responses, function(response) {
        is.logistic <- length(unique(pull(dset, response))) < min_n_for_continuous
        fo.family <- ifelse(is.logistic, stats::binomial, stats::gaussian)
        fo <- sprintf("%s ~ %s", response, paste0(covariates, collapse = "+"))
        stats::glm(formula = as.formula(fo), family = fo.family, data = dset) %>%
            broom::tidy(conf.int = TRUE, conf.level = 0.95, exponentiate = is.logistic) %>%
                dplyr::filter(grepl(filterstr, term)) %>%
                dplyr::mutate(response = response, fo = fo)
    })
    data.table::rbindlist(glmlist, id = "model_name") %>%
        as.data.frame
}


calculateadonis <- function(dset,
                            matrix,
                            responses = list(model_1 = "MAP", model_2 = "SYSTM", model_3 = "DIASM",
                                          model_4 = "PULSEPRESSURE", model_5 = "HYPERTENSION"),
                            covariates = ".",
                            npermutations = 99,
                            maxcores = 100) {
    mclapply(responses, function(response) {
        fo <- sprintf("matrix ~ %s + %s", paste(covariates, collapse = " + "), response)
        ad <- adonis(formula = as.formula(fo), data = dset, permutations = npermutations)
        ad$aov.tab %>%
            as.data.frame %>%
            tibble::rownames_to_column(var = "term") %>%
            dplyr::mutate(response = response)
      }, mc.cores = min(maxcores, length(responses)))
}

merge.pheno.abu <- function(pseq,
                             core = FALSE,
                             core.detection = 0.001,
                             core.prevalence = 0.01) {
     pseq.comp <- microbiome::transform(pseq, "compositional")
     if (core) {
         pseq.comp  <- microbiome::core(pseq.comp, detection = core.detection, prevalence = core.prevalence)
    }
     pheno <- microbiome::meta(phyloseq::sample_data(pseq.comp))
     sample.ids <- rownames(pheno)
     abud.rel <- abundances(pseq.comp)
     rownames(abud.rel) <- sapply(rownames(abud.rel), replace.brackets)
     cbind(pheno, t(abud.rel))
}


legend_extractor <- function (tmp) {
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}

ylabl_extractor <- function(tmp) {
  yl <- which(grepl('axis.title.y.left', sapply(tmp$grobs, function(x) x$name)))
  tmp$grobs[[yl]]
}

xlabb_extractor <- function(tmp) {
  xl <- which(grepl('axis.title.x.bottom', sapply(tmp$grobs, function(x) x$name)))
  tmp$grobs[[xl]]
}

vectortolist <- function(c) {
  l <- as.list(c)
  names(l) <- l
  l
}

getdescriptions <- function() {
    names.dset <- readRDS("mydata/phfinrisk_metadatadesc.rds") %>%
        dplyr::filter(Ignored.Covariate.in.Cross.Sectional.Analysis.Aaro20181115==0)
    bind_rows(names.dset,
              data.frame(Covariate = c("MAP", "PULSEPRESSURE", "HYPERTENSION", "SEX"),
                         Category = rep("Physical", 4),
                         Name = c("MAP", "Pulse pressure", "Hypertension", "Female"),
                         Desc = c("Mean arterial presure, 2./3.*DIASM + 1./3.*SYSTM)",
                                  "Pulse pressure, SYSTM - DIASM",
                                  "Hypertension, SYSTM >=140 or DIAS >= 90 or BP mediaction",
                                  "Sex is female, True for female, False for male")))
}

mygrep <- function(..., word, ignorecase = TRUE, complement = FALSE) {
    c(...)[xor(grepl(word, c(...), ignore.case = ignorecase), (complement == TRUE))]
}


myin <- function(x,y, complement = FALSE) {
    x[xor(x %in% y, complement)]
}

cc <- function(...) {
    c2l(c(...))
}


c2l <- function(...) {
    l <- as.list(c(...))
    names(l) <- c(...)
    l
}

pub.p <- function(p, Nfdr = FALSE, tiny = TRUE) {
    p <- as.numeric(p)
    if (Nfdr) p <- p.adjust(p, method="BH", n = Nfdr)
    case_when(tiny & p < 0.001 ~ "<0.001",
              p < 0.001 ~ sprintf("%.3e", p),
              p >= 0.999 ~ "1",
              TRUE ~ sprintf("%.3f", p))
}

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

mykable <- function(x, ...) {
  capture.output(x <- print(x))
  knitr::kable(x, ...)
}

myinstall.packages <- function(...) {
    list.of.packages <- c(...)
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if (length(new.packages) == 0) { return(TRUE) }
    for (package in new.packages) {
        message(sprintf("Installing: %s", package))
        myinstall.packages(gtools::getDependencies(package))
        install.packages(package)
    }
}


mydropna <- function(...) {
    c(...)[!is.na(c(...))]
}


diversities <- function(pseq, vars, betadiversity, names.dset) {
    lapply(c2l(names(vars)), function(var) {
        alphadiversity <- calculateglm(meta.merge.alphadiversity(pseq),
                                       covariates = c("diversity_shannon", vars[[var]]),
                                       filterstr = "shannon") %>%
            select(response, alpha.p = p.value, alpha.effect = estimate,
                   alpha.low = conf.low, alpha.high = conf.high)
        betadiversity <- betadiversity[[var]] %>%
            map_df(., ~as.data.frame(.x %>%
                                     dplyr::filter(term %in% var.BP) %>%
                                     select(response, beta.R2 = R2, beta.p=`Pr(>F)`)))
        full_join(alphadiversity, betadiversity, by = "response") %>%
            merge(names.dset %>% select(Covariate, Name), by.x = "response", by.y = "Covariate")
    })
}

diversities.tidy <- function(diversity) {
    map_df(diversity, ~as.data.frame(.x), .id = "covariates") %>%
        mutate(alpha = sprintf("%.2f (%.2f–%.2f)", alpha.effect, alpha.low, alpha.high),
               alpha.p = pub.p(alpha.p),
               beta.R2 = sprintf("%.3f%%", 100*beta.R2),
               beta.p = pub.p(beta.p)) %>%
        select(Name, covariates, alpha, alpha.p, beta.R2, beta.p) 
}

calculate.betadiversity <- function(pseq, matrix, vars, npermutations = 999, maxcores = 10) {
    mclapply(vars, function(var) {
        calculateadonis(dset = meta(pseq),
                        matrix = matrix,
                        covariates = var,
                        npermutations = npermutations)
    }, mc.cores = min(maxcores, length(vars)))
}

deseqresults <- function(modellist,
                         names.dset,
                         p = 0.05,
                         hypertensionvars = c("HYPERTENSION")) {
    vars <- c2l(names(modellist))
    lapply(c2l(vars), function(x, models = modellist) {
        name <- ifelse(x %in% hypertensionvars, sprintf("HYPERTENSION_1_vs_0", x), x)
        DESeq2::results(models[[x]], name = name) %>%
            as.data.frame %>%
            tibble::rownames_to_column("Feature") }) %>%    
        map_df(., ~as.data.frame(.x), .id="name") %>%
        dplyr::mutate(qval = p.adjust(pvalue, method="BH"),
                      Feature = renametaxa(Feature)) %>%
        dplyr::filter(qval < p) %>%
        merge(names.dset %>% select(Covariate, Name), by.x ="name", by.y = "Covariate")
}

prune_lactobacillus <- function(pseq, transform = "compositional", drop = TRUE) {
    microbiome::transform(pseq, transform = transform) %>%
        prune_taxa(c("g_Lactobacillus"), .) %>%
        psmelt %>%
        spread(OTU, Abundance) %>%
        dplyr::mutate(lacto.prevalent = ifelse(g_Lactobacillus > 0.1/100, 1, 0),
                      duna.present = ifelse(is.na(NA.), 0, 1)) %>%
        { if (drop == TRUE) dplyr::filter(., duna.present == 1) else . }
}

pseq_prevalence <- function(pseq,
                            taxa = "Lactobacillus (Bacteria)",
                            limit = 0.1/100,
                            fun = mean) {
    pseq %>%
        microbiome::transform(transform = "compositional") %>%
        abundances %>%
        as.data.frame %>%
        .[taxa, ] %>%
        { ifelse(. > limit, 1, 0) } %>%
        fun
}

pseq_prevalence_format <- function(pseq, taxa) {
    sprintf("%.1f%%", 100*pseq_prevalence(pseq, taxa))
}

mytableone <- function(dset, variables) {
    { if (!missing(variables)) dplyr::select(dset, variables) else dset } %>%
        select(-Sample) %>% 
        gather(key,value) %>%
        group_by(key) %>%
        summarise_all(funs(N = n(),
                           Mean = mean(as.numeric(value), na.rm = TRUE),
                           ones = sum(value == 1, na.rm = TRUE),
                           sd = sd(as.numeric(value), na.rm = TRUE),
                           NAs = sum(is.na(value))))  %>%
        mutate_if(is.numeric, round, 3)
}

deseq.list <- function(var.CL, var.CL.min) {
    l <- lapply(c2l(var.CL), function(x) {
        unique(c(var.CL.min, x))
    })
    l[!duplicated(l)]
}

coretaxa <- function(pseq, detection = 0.1/100, prevalence = 1/100) {
    microbiome::transform(pseq, "compositional") %>%
        core(detection = detection, prevalence = prevalence) %>%
        taxa
}

deseq.formula <- function(..., fo = "~ %s") {
    as.formula(sprintf(fo, paste0(c(...), collapse = "+")))
}

pseqsubset <- function(pseq, coretaxa, saltsubset = TRUE, FUN = identity) {
    exists.coretaxa <- !missing(coretaxa)
    participants <- pseq %>% meta %>% tibble::rownames_to_column("sampleid") %>%
        { if (saltsubset) dplyr::filter(., !is.na(dUNA)) else . } %>% FUN %>% pull(sampleid)
    prune_samples(participants, pseq) %>%
        { if (exists.coretaxa) prune_taxa(coretaxa, .) else . }
}

myDESeq <- function(pseq, vars, coretaxa, coreterm = ".", saltsubset = FALSE, FUN = identity) {
    mydds.pseq <- pseqsubset(pseq, coretaxa = coretaxa, saltsubset = saltsubset, FUN = FUN)
    mydds.formula <- deseq.formula(vars)
    mydds.data <- phyloseq_to_deseq2(mydds.pseq, mydds.formula)
    mydds.size <- estimateSizeFactors(mydds.data)
    mydds.dispersion <- mydds.size[mygrep(word = coreterm, taxa(mydds.pseq)),] %>%
        estimateDispersions(., fitType="parametric")
    nbinomWaldTest(mydds.dispersion)
}

myDESeq.tidy <- function(dds, plimit = 1.0) {
    results(dds, name = "dUNA", tidy = TRUE) %>%
        dplyr::mutate(qval = p.adjust(pvalue, method="BH"),
                      row = renametaxa(row),
                      lfc_se = sprintf("%.3f±%.3f", log2FoldChange, lfcSE),
                      qval = pub.p(qval, tiny = FALSE)) %>%
        filter(qval < p.limit) %>%
        select(row, lfc_se, qval)
}

myDESeq.subsets <- function(var) {
    list("female" = list(fun = function(x) subset(x, SEX == 1),
                         exclude = "SEX"),
         "male" = list(fun = function(x) subset(x, SEX == 0),
                       exclude = "SEX"),
         "drug" = list(fun = function(x) subset(x, ANYDRUG == 1),
                       exclude = c("BL_USE_RX_C03",
                                   "BL_USE_RX_C07",
                                   "BL_USE_RX_C08",
                                   "BL_USE_RX_C09")),
         "no-drug" = list(fun = function(x) subset(x, ANYDRUG == 0),
                          exclude = c("BL_USE_RX_C03",
                                      "BL_USE_RX_C07",
                                      "BL_USE_RX_C08",
                                      "BL_USE_RX_C09")))
}

cleanemptyfactors <- function(df, factors, FUN = function(x) filter(x, .x == 1)) {
    exclude <- lapply(c2l(factors), function(x) df[[x]] %>% factor %>% nlevels) %>%
        map_df(~as.data.frame(.x), .id = "var") %>%
        FUN %>%
        pull(var)
    factors %difference% exclude
}

deseq.abundances <- function(dds) {
    assays(dds)[["mu"]] %>%
        as.data.frame %>%
        tibble::rownames_to_column("rowname") %>%
        gather(sampleid, value, -rowname) %>%
        spread(rowname, value) %>%
        rename_all(replace.brackets) 
}

loop.lm <- function(dset,
                    response,
                    loops,
                    covariates = c()) {
    models <- lapply(c2l(loops), function(loop) {
        fo <- sprintf("%s ~ %s", response, paste(c(loop, covariates), collapse = " + "))
        ret <- stats::lm(formula = as.formula(fo), data = dset, na.action = na.omit)
        ret$call <- as.formula(fo)
        ret
    })
}

loop.binomial <- function(dset,
                    response,
                    loops,
                    covariates = c()) {
    stopifnot(!missing(dset), !missing(response), !missing(loops))
    mclapply(c2l(loops), function(loop) {
        fo <- sprintf("%s ~ %s", response, paste(c(loop, covariates), collapse = " + "))
        ret <- stats::glm(formula = as.formula(fo),
                          family=binomial(link='logit'),
                          data = dset,
                          na.action = na.omit)
        ret$call <- as.formula(fo)
        ret
    }, mc.cores = min(length(loops), 8))
}

loop.results <- function(..., filterstr = "^K[0-9]*", exponentiate = FALSE) {
    purrr::map_df(c(...), ~tidy(.x, exponentiate = exponentiate)) %>%
        dplyr::filter(grepl(filterstr, term)) %>%
        dplyr::mutate(conf.low = estimate - qnorm(0.975) * std.error,
                      conf.high = estimate + qnorm(0.975) * std.error,
                      qval = p.adjust(p.value, method="BH"))
}

mychisq.test <- function(x, y) {
    table(x, y) %>% chisq.test() 
}

mykeggget <- function(list, index = "DEFINITION") {
    mclapply(list, function(x) keggGet(x)[[1]][[index]], mc.cores = 8) %>%
        unlist %>%
        gsub(" \\[.*\\]", "", .)
}
