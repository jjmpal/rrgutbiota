getfactorvariables <- function(df, vars) {
    classes <- lapply(c2l(vars), function(x) class(df[[x]]))
    classes[classes == "factor"] %>% names
}

characteristics.names <- function(onlyvars = FALSE) {
    if (onlyvars == TRUE)
        return(c("BL_AGE", "SEX", "BMI", "oSYSTM", "oDIASM",
                 "oPULSEPRESSURE", "oMAP", "HYPERTENSION",  
                 "CURR_SMOKE", "PREVAL_DIAB",  "Q57X", "BL_USE_RX_C03",
                 "BL_USE_RX_C07", "BL_USE_RX_C08", "BL_USE_RX_C09"))
    list("BL_AGE" = "Age, y (SD)",
         "SEX" = "Female, N (%)",
         "BMI" = "BMI, kg/m² (SD)",
         "oSYSTM" = "Systolic BP, mmHg (SD)",
         "oDIASM" = "Diastolic BP, mmHg (SD)",
         "oPULSEPRESSURE" = "Pulse pressure, mmHg (SD)",
         "oMAP" = "Mean arterial pressure, mmHg (SD)",
         "HYPERTENSION" = "Hypertension, N (%)",
         "CURR_SMOKE" = "Current smoker, N (%)",
         "PREVAL_DIAB" = "Diabetes mellitus, N (%)",
         "Q57X" = "Exercise, N (%)",
         "ANYDRUG" = "Antihypertensive medication, N (%)",
         "BL_USE_RX_C03" = "  Diuretics, N (%)",
         "BL_USE_RX_C07" = "  Beta blockers, N (%)",
         "BL_USE_RX_C08" = "  Calcium channel blockers, N (%)",
         "BL_USE_RX_C09" = "  RAS blockers, N (%)",
         "1"  =  "  Light",
         "2"  =  "  Moderate",
         "3"  =  "  Heavy",
         "4"  = "  Competitive")
}

characteristicsTable <- function(dset, strata) {
    nostrata <- missing(strata)
    characteristics <- tableone::CreateTableOne(
                                     strata = strata,
                                     data = dset,
                                     vars = characteristics.names(TRUE),
                                     factorVars = getfactorvariables(dset,
                                                                     characteristics.names(TRUE)),
                                     test = !missing(strata))
    print(characteristics,
                      exact = "stage",
                      quote = FALSE,
                      noSpaces = TRUE,
                      printToggle = FALSE,
                      digits = 1,
                      pDigits = 3,
                      contDigits=1)  %>%
        as.data.frame %>%
        tibble::rownames_to_column(var = "rowname") %>%
        format(justify = "left", trim = TRUE) %>%
        mutate(rowname = characteristics.names()[gsub("^ *([A-Za-z_0-9]+).*", "\\1", rowname)]) %>%
        { if (!nostrata) select(., -test) else . }
}

characteristics <- function(dset, tableone.names, tableone.factors, extras = list()) {
    title <- "Characteristics"
    overall <- paste0("Cases, n=", dim(dset)[1])
    tableobject <- tableone::CreateTableOne(data = dset, vars = names(tableone.names), factorVars = tableone.factors)
    tablecsv <- print(tableobject,
                      exact = "stage",
                      quote = FALSE,
                      noSpaces = TRUE,
                      printToggle = FALSE,
                      digits = 1,
                      pDigits = 3,
                      contDigits=1)

    tableone.fullnames <- c(tableone.names, extras)
    tablecsv %>%
        as.data.frame %>%
        tibble::rownames_to_column(var = "rowname") %>%
        dplyr::filter(row_number() > 1) %>%
        format(justify = "left", trim = TRUE) %>%
        rowwise() %>%
        mutate(id = gsub("^ *([A-Za-z_0-9]+).*", "\\1", rowname)) %>%
        mutate(present = id %in%  names(tableone.fullnames)) %>%
        mutate(rowname = ifelse(present == TRUE, tableone.fullnames[[id]], rowname)) %>%
        select(rowname, Overall)
}

typologyformatter <- function(data, font = 12, typology, left = c(1), hleft = c(1)) {
  flex <- flextable(data = data) %>%
    flextable::theme_booktabs() %>%
    flextable::border(border = fp_border(width=0), part="body") %>%
    flextable::border(border = fp_border(width=0), part="header") %>%
    flextable::border(part="header", border.bottom = fp_border(width=1))

  if (!missing(typology)) {
      flex <- flex %>%
          set_header_df(mapping = typology, key = "col_keys") %>%
          merge_h(part = "header") %>%
          flextable::border(part="header", border.bottom = fp_border(width=1))
      if (missing(hleft)) {
          hleft <- c(2)
      }
  }

  flex %>%
      flextable::border(i = nrow(data), part="body", border.bottom = fp_border(width=1)) %>%
      flextable::bold(bold = FALSE, part = "header") %>%
      flextable::bold(bold = FALSE, part = "body") %>%
      flextable::fontsize(size = font, part = "header") %>%
      flextable::fontsize(size = font, part = "body") %>%
      flextable::align(align = "center", part = "all") %>%
      flextable::align(align = "left", part = "header", j = left, i = hleft) %>%
      flextable::align(align = "left", part = "body", j = left)
}

myspread <- function(ret, list = c2l("lfc_se", "p.value"), term  = "Feature", key = "Name") {
    lapply(list, function(x) ret %>%
                             select(term, key, x) %>%
                             spread(key, x) %>%
                             rename_at(vars(-term), funs(paste0(., "_", x)))) %>%
        Reduce(function(...) full_join(..., by = term), .) %>%
        dplyr::select(term, noquote(order(colnames(.))))
}
