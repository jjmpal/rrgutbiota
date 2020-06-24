import_filter_data <- function(file, included = c("MEN",
                                          "BL_AGE",
                                          "SYSTM",
                                          "DIASM",
                                          "BMI",
                                          "CURR_SMOKE",
                                          "PREVAL_DIAB",
                                          "Q57X",
                                          "BL_USE_RX_C09",
                                          "BL_USE_RX_C03",
                                          "BL_USE_RX_C07",
                                          "BL_USE_RX_C08",
                                          "NA.",
                                          "BP_TREAT",
                                          "PREVAL_CHD",
                                          "PREVAL_CR_ANYCANC")) {
    pseq.full <- readRDS(file)
    pseq.meta <- pseq.full %>%
        meta %>%
        as.data.frame %>%
        phenotype.filter(included = included) %>%
        phenotype.definitions 
    phyloseq::sample_data(pseq.full) <- phyloseq::sample_data(pseq.meta)
    return(pseq.full)
}

phenotype.extra <- function(extra, fextra = "pheno/2015_60_Salomaa_Jain_dataFR02_FU17_2020-02-19.txt.gz") {
    read_tsv(fextra, col_types = cols(
                         .default = col_double(),
                         CDT = col_character(),
                         Sample_ID = col_character(),
                         Barcode = col_character(),
                         WESID = col_character(),
                         WGSID = col_character(),
                         BATCH = col_character(),
                         FID = col_character(),
                         CRP = col_character(),
                         K_TPKS = col_character(),
                         K_VKS = col_character(),
                         K_M1 = col_character(),
                         K_M2 = col_character(),
                         APOE_BATCH = col_character())) %>%
        select(one_of(extra)) %>%
        filter(!is.na(Barcode))
}

phenotype.filter <- function(df,
                            included,
                             extra = c("Barcode", "U_VIRTSA_YHT"),
                             allowna = c("NA.",
                                         "U_VIRTSA_YHT",
                                         "BP_TREAT",
                                         "PREVAL_CHD",
                                         "PREVAL_CR_ANYCANC")) {
    df %>%
        tibble::rownames_to_column(var = "rowname") %>%
        dplyr::select(rowname, included) %>%
        left_join(phenotype.extra(extra), by = c("rowname" = "Barcode")) %>%
        tidyr::drop_na(myin(included, allowna, complement = TRUE))
}

phenotype.definitions <- function(df) {
    df %>%
        dplyr::mutate(ANYDRUG = factor(ifelse(BL_USE_RX_C03  == 1 | BL_USE_RX_C07 == 1 |
                                              BL_USE_RX_C08 == 1  | BL_USE_RX_C09 == 1, 1, 0)),
                      ANYEXERCICE = factor(ifelse(Q57X == 1, 0, 1)),
                      PULSEPRESSURE = SYSTM - DIASM,
                      HYPERTENSION = factor(ifelse(SYSTM >= 140 | DIASM >= 90 | ANYDRUG == 1, 1, 0)),
                      SEX = factor(ifelse(MEN == "Female", 1, 0)),
                      MAP = 2./3.*DIASM + 1./3.*SYSTM) %>%
        dplyr::mutate(oSYSTM = SYSTM,
                      oDIASM = DIASM,
                      oPULSEPRESSURE = PULSEPRESSURE,
                      oMAP = MAP,
                      SYSTM = myscale(SYSTM),
                      DIASM = myscale(DIASM),
                      PULSEPRESSURE = myscale(PULSEPRESSURE),
                      MAP = myscale(MAP),
                      SEX = as.factor(SEX),
                      Q57X = factor(Q57X, ordered = FALSE),
                      dUNA = NA.*U_VIRTSA_YHT/1000,
                      BL_USE_RX_C03 = as.factor(BL_USE_RX_C03),
                      BL_USE_RX_C07 = as.factor(BL_USE_RX_C07),
                      BL_USE_RX_C08 = as.factor(BL_USE_RX_C08),
                      BL_USE_RX_C09 = as.factor(BL_USE_RX_C09)) %>%
        dplyr::select(-MEN) %>%
        tibble::remove_rownames() %>%
        tibble::column_to_rownames(var = "rowname")
}
