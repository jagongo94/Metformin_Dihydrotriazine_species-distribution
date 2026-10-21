# Metformin_Dihydrotriazine_species-distribution
Metformin dihydrotriazine conjugates species distribution from masst results and HNRC reanalysis
# =============================================================================
#  Metformin-carbonyl conjugates in humans - figure codebook, panels a to e
#  Publication-ready R
# -----------------------------------------------------------------------------
#  Two independent datasets, five panels, one script.
#
#  COHORT  panels a-d        merged_MS_metadata_matched_metformin_exposure.xlsx
#      Long format, one row per molecule x sample injection, carrying the MS
#      peak area and the full clinical metadata of that sample.
#        a  detection rate per conjugate, metformin-exposed vs unexposed
#        b  detection by reported diabetes, with odds ratio and P value
#        c  per-sample detection heatmap, log(peak area + 1)
#        d  chord diagram, conjugate to metabolic risk factor
#
#  REPOSITORY  panel e       conjugate_library_matches.txt
#      Tab-separated library-match export from public MS/MS repositories.
#        e  occurrence of each conjugate by host class
#
#  OUTPUT   every panel as .svg, .pdf and .png, the plotted numbers as .csv,
#           and panel_statistics.txt.
#
#  All headings, thresholds and category lists live in the CONFIG block.
#  Nothing below it writes text or decides which compounds appear.
#
#  DEFINITIONS, fixed once so every panel agrees
#    Sample       one biological sample, keyed on ID. Repeat injections collapse
#                 to the maximum area, so the cohort is n = 314, not 315.
#    Detected     Area_sub_corrected > 0.
#    Odds ratio   ad/bc; Haldane-Anscombe (+0.5 to every cell) only when a cell
#                 is zero, otherwise the plain estimate.
#    P value      Fisher's exact test, two-sided, raw and uncorrected.
#    Panel e fill log(number of library matches + 1), natural log.
#
#  Retention-time confirmation (the asterisks in panel e) is not in either data
#  file - it comes from authentic standards run in-house and is listed in
#  RT_VALIDATED below.
#
#  Panel d needs the 'circlize' package. Without it the script reports the fact,
#  still writes the chord link table, and carries on.
# =============================================================================

## ---- 0. Packages ------------------------------------------------------------
# install.packages(c("readxl","readr","dplyr","tidyr","stringr","ggplot2",
#                    "svglite","scales","circlize"))
suppressPackageStartupMessages({
  library(readxl); library(readr); library(dplyr); library(tidyr)
  library(stringr); library(ggplot2); library(svglite); library(scales)
})

## =============================================================================
##  1. CONFIGURATION
## =============================================================================
PANELS <- c("a", "b", "c", "d", "e")   # subset this to rebuild only some panels

COHORT_FILE <- "merged_MS_metadata_matched_metformin_exposure.xlsx"
COHORT_SHEET <- 1
REPO_FILE   <- "conjugate_library_matches.txt"

OUT_DIR     <- "figures"
DEVICES     <- c("svg", "pdf", "png")
PNG_DPI     <- 600
BASE_SIZE   <- 9
BASE_FAMILY <- ""          # "" = device default sans; set "Helvetica"/"Arial" if installed

## Thresholds for risk factors recorded as measurements rather than flags
HBA1C_CUTOFF <- 6.5        # % glycosylated haemoglobin
BMI_CUTOFF   <- 30         # kg/m2, obesity

## Panel sizes, inches
SIZE <- list(a = c(7.0, 5.6), b = c(6.4, 5.6), c = c(10.5, 4.6),
             d = c(7.2, 7.2), e = c(11.0, 3.4))

## Panel e: what one unit of "count" means.
##   "spectra" one row of the export = one library match   (used for the published panel)
##   "files"   distinct source files; "usi" distinct spectra by identifier
COUNT_UNIT <- "spectra"
PRE_FILTER      <- FALSE   # the published panel used the export as issued
MIN_COSINE      <- 0.80
MIN_MATCH_PEAKS <- 5
DROP_UNCLASSIFIED <- TRUE  # rows with no assignable host class

## ---- 1b. Figure text --------------------------------------------------------
FIG <- list(
  a = list(tag = "a",
           title     = "Conjugates are detected almost exclusively in metformin users",
           subtitle  = "Fisher's exact test, exposed vs unexposed. *** P < 0.001  ** P < 0.01  * P < 0.05",
           x         = "Detection rate (%)", y = NULL,
           exposed   = "Metformin-exposed",
           unexposed = "Unexposed"),
  b = list(tag = "b",
           title    = "Detection by reported diabetes",
           subtitle = "Odds ratio and two-sided Fisher's exact P value, uncorrected.",
           cols     = c("Yes", "No", "Odds\nratio", "P value"),
           fill     = "Odds ratio"),
  c = list(tag = "c",
           title    = "Per-sample detection",
           subtitle = "Samples with at least one conjugate detected.",
           x = NULL, y = NULL,
           fill     = "log(area + 1)"),
  d = list(tag = "d",
           title    = "Conjugate detection and metabolic risk factors",
           subtitle = "Ribbon width is the number of samples positive for both.",
           legend   = "Risk factors"),
  e = list(tag = "e",
           title    = "Metformin-carbonyl conjugates in public untargeted metabolomics data",
           subtitle = "Library matches across public MS/MS repositories, by host class.",
           x = NULL, y = NULL,
           fill     = "Log(count + 1)",
           star     = "Further validated by\nretention time")
)

## ---- 1c. Categories ---------------------------------------------------------
## Risk factors: display name -> colour. Order sets the arc order in panel d.
RISK <- c("HbA1c > 6.5"        = "#C0267E",
          "High blood sugar"   = "#2B4FA2",
          "High triglycerides" = "#E43FA0",
          "Abdominal obesity"  = "#8B3FC4",
          "Diabetes"           = "#0E7C7B",
          "Low HDL"            = "#4BC6D6",
          "Insulin"            = "#1E8A4C",
          "BMI"                = "#7E1A1A")

## Panel e rows, top to bottom
TAXON_ORDER <- c("Homo sapiens", "Rodent", "Plant", "Microbe")

## Panel e columns: straight-chain aldehydes by carbon number, then unsaturated
## and long-chain, small oxo-acids, terpenoids, aromatics, steroids.
COMPOUND_ORDER <- c(
  "Acetaldehyde", "Propanal", "Butyraldehyde", "Crotonaldehyde", "Pentanal",
  "Hexanal", "Heptanal", "Octanal", "Nonanal", "Decanal", "Undecanal",
  "Dodecanal", "Tridecanal", "Tetradecanal", "Pentadecanal", "Hexadecanal",
  "Heptadecanal", "Octadecanal", "9-Octadecenal", "Icosanal",
  "Glyceraldehyde", "6-oxo-Hexanoic acid", "Citral", "Safranal",
  "Trans-Cinnamaldehyde", "3-Phenylacetaldehyde", "3-Phenylpropionaldehyde",
  "Indole-3-Carbaldehyde", "Vanillin", "Syringaldehyde",
  "3-oxo-Lithocholicacid", "3-oxo-Deoxycholicacid", "3-oxo-Cholestanone",
  "3-oxo-Cholest-5-ene", "3-oxo-sitostanone")

## Conjugates confirmed against an authentic standard by retention time.
## Not derivable from either data file - maintained by hand.
RT_VALIDATED <- c(
  "Acetaldehyde", "Pentanal", "Hexanal", "Heptanal", "Octanal", "Nonanal",
  "Decanal", "Undecanal", "Dodecanal", "Tridecanal", "Tetradecanal",
  "Pentadecanal", "Hexadecanal", "Heptadecanal", "Octadecanal",
  "9-Octadecenal", "Icosanal", "Citral", "Safranal", "Indole-3-Carbaldehyde",
  "Vanillin", "3-oxo-Lithocholicacid", "3-oxo-Deoxycholicacid",
  "3-oxo-Cholestanone", "3-oxo-Cholest-5-ene", "3-oxo-sitostanone")

## Spellings in the repository export that do not normalise onto COMPOUND_ORDER
COMPOUND_ALIASES <- c(
  "phenylacetaldehyde"      = "3-Phenylacetaldehyde",
  "3-phenylpropionaldehyde" = "3-Phenylpropionaldehyde",
  "indole-3-carbaldehyde"   = "Indole-3-Carbaldehyde",
  "trans-cinnamaldehyde"    = "Trans-Cinnamaldehyde",
  "3-oxo-cholestane"        = "3-oxo-Cholestanone",
  "3-oxo-cholest-5-ene"     = "3-oxo-Cholest-5-ene",
  "3-oxo-lithocholicacid"   = "3-oxo-Lithocholicacid",
  "3-oxo-deoxycholicacid"   = "3-oxo-Deoxycholicacid",
  "3-oxo-sitostanone"       = "3-oxo-sitostanone",
  "6-oxo-hexanoic-acid"     = "6-oxo-Hexanoic acid")

## ---- 1d. Colours ------------------------------------------------------------
COL_EXPOSED   <- "#B5202C"
COL_UNEXPOSED <- "#4E8FD1"
OR_LOW        <- "#F2F6FC"; OR_HIGH  <- "#2F5FA8"
AREA_LOW      <- "#FBE9E7"; AREA_HIGH <- "#C0392B"
FILL_LIMITS   <- c(0, 10)   # panel e, deliberately wider than the data
FILL_MIDPOINT <- 5
FILL_LOW      <- "#4A77B4"; FILL_MID <- "#F7F7F7"; FILL_HIGH <- "#D6604D"
FILL_NA       <- "#FAFAFA"  # host class searched, nothing found

## =============================================================================
##  2. SHARED HELPERS
## =============================================================================
dir.create(OUT_DIR, showWarnings = FALSE)
want <- function(p) p %in% PANELS

theme_pub <- function(base_size = BASE_SIZE, base_family = BASE_FAMILY) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(plot.title    = element_text(size = base_size + 2, face = "bold", hjust = 0,
                                       margin = margin(b = 2)),
          plot.subtitle = element_text(size = base_size - 1, colour = "#555555", hjust = 0,
                                       margin = margin(b = 6)),
          plot.tag      = element_text(size = base_size + 5, face = "bold"),
          plot.tag.position = c(0.005, 0.99),
          axis.line     = element_line(linewidth = 0.4),
          axis.ticks    = element_line(linewidth = 0.4),
          axis.text     = element_text(colour = "black"),
          legend.position = "bottom",
          legend.title  = element_text(size = base_size - 1),
          legend.text   = element_text(size = base_size - 1),
          plot.background = element_rect(fill = "white", colour = NA),
          plot.margin   = margin(6, 8, 4, 6))
}

## Panel e is a grid rather than a plot with axes, so it drops the axis furniture.
theme_grid <- function(base_size = BASE_SIZE) {
  theme_pub(base_size) +
    theme(panel.grid = element_blank(), axis.line = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = base_size - 1.5),
          axis.text.y = element_text(size = base_size, colour = "black"),
          legend.justification = "left", legend.box = "horizontal",
          legend.key.height = unit(0.32, "lines"), legend.key.width = unit(1.5, "lines"))
}

save_panel <- function(p, name, wh) {
  for (dev in DEVICES) {
    f <- file.path(OUT_DIR, paste0(name, ".", dev))
    if (dev == "svg")      ggsave(f, p, width = wh[1], height = wh[2], device = svglite::svglite)
    else if (dev == "pdf") ggsave(f, p, width = wh[1], height = wh[2], device = cairo_pdf)
    else                   ggsave(f, p, width = wh[1], height = wh[2], dpi = PNG_DPI, bg = "white")
  }
  message("  wrote ", name)
}

LOG <- character(0)
log_line <- function(...) { s <- paste0(...); LOG <<- c(LOG, s); message(s) }

stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***",
                     ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

## Odds ratio and Fisher P from a 2x2: a,b = group 1 detected / not; c,d = group 2.
or_fisher <- function(a, b, c, d) {
  zero <- any(c(a, b, c, d) == 0)
  or <- if (zero) ((a + .5) * (d + .5)) / ((b + .5) * (c + .5)) else (a * d) / (b * c)
  list(or = or, p = fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE))$p.value,
       haldane = zero)
}

num <- function(x) suppressWarnings(as.numeric(as.character(x)))

normalise_key <- function(x) {
  x %>% str_replace("^[Mm]etformin[-_ ]*", "") %>%
    str_replace_all("[_ ]+", "-") %>% str_replace_all("-+", "-") %>%
    str_to_lower() %>% str_trim()
}

## Host class from NCBI taxonomy fields. Extend the patterns as the search widens.
classify_host <- function(taxonomy, kingdom, family, species) {
  tx <- paste(taxonomy, kingdom, family, species)
  case_when(
    str_detect(tx, regex("homo sapiens", ignore_case = TRUE))                 ~ "Homo sapiens",
    str_detect(tx, regex("\\bMus\\b|\\bRattus\\b|Muridae|Cricetidae|Rodentia",
                         ignore_case = TRUE))                                 ~ "Rodent",
    str_detect(tx, regex("Viridiplantae|Streptophyta|Magnoliopsida",
                         ignore_case = TRUE))                                 ~ "Plant",
    str_detect(tx, regex("Bacteria|Archaea|Fungi|Firmicutes|Proteobacteria|Bacteroidetes|Actinobacteria",
                         ignore_case = TRUE))                                 ~ "Microbe",
    TRUE                                                                      ~ NA_character_)
}

## =============================================================================
##  3. COHORT DATA  (panels a-d)
## =============================================================================
if (any(want(c("a", "b", "c", "d")))) {

  raw <- read_excel(COHORT_FILE, sheet = COHORT_SHEET, guess_max = 20000,
                    .name_repair = "minimal")
  names(raw) <- make.unique(names(raw))
  insulin_cols <- grep("^on_medication_.*insulin", names(raw), value = TRUE, ignore.case = TRUE)

  dat <- raw %>%
    transmute(
      id = as.character(ID), molecule = as.character(Molecule),
      area = num(Area_sub_corrected),
      exposure = as.character(metformin_exposure),
      diabetes = as.character(diabetes),
      hba1c = num(glycosylated_hemoglobin_a1c), bmi = num(bmi),
      abd_obesity = as.character(abdominal_obesity_for_mets),
      high_sugar  = as.character(high_blood_sugar_glucose_for_mets),
      high_trig   = as.character(high_trig_hypertriglyceridemia_for_mets),
      low_hdl     = as.character(low_hdl_dislipidemia_for_mets),
      insulin = if (length(insulin_cols))
                  apply(raw[insulin_cols], 1, function(v)
                    any(!is.na(v) & !(as.character(v) %in% c("0", "NA", "")))) else FALSE) %>%
    filter(!is.na(molecule), !is.na(id)) %>%
    mutate(area = ifelse(is.na(area), 0, area)) %>%
    group_by(id, molecule) %>% slice_max(area, n = 1, with_ties = FALSE) %>% ungroup() %>%
    mutate(detected = area > 0)

  log_line(sprintf("cohort: %d samples x %d conjugates after collapsing repeat injections",
                   n_distinct(dat$id), n_distinct(dat$molecule)))

  ## Conjugates ordered once, by detection rate among the exposed; panels a, b
  ## and c all reuse this order.
  MOL_ORDER <- dat %>% filter(exposure == "Yes") %>%
    group_by(molecule) %>% summarise(rate = mean(detected), .groups = "drop") %>%
    arrange(desc(rate), molecule) %>% pull(molecule)
  dat <- mutate(dat, molecule = factor(molecule, levels = MOL_ORDER))

  n_exp <- n_distinct(dat$id[dat$exposure == "Yes"])
  n_unx <- n_distinct(dat$id[dat$exposure == "No"])

  panel_a_tbl <- dat %>% filter(exposure %in% c("Yes", "No")) %>%
    group_by(molecule) %>%
    summarise(det_exp = sum(detected & exposure == "Yes"),
              det_unx = sum(detected & exposure == "No"), .groups = "drop") %>%
    rowwise() %>%
    mutate(rate_exp = 100 * det_exp / n_exp, rate_unx = 100 * det_unx / n_unx,
           P = or_fisher(det_exp, n_exp - det_exp, det_unx, n_unx - det_unx)$p,
           sig = stars(P)) %>% ungroup()
}

## ---- Panel a ----------------------------------------------------------------
if (want("a")) {
  write_csv(panel_a_tbl, file.path(OUT_DIR, "panel_a_detection_rates.csv"))
  lab_exp <- sprintf("%s (n = %d)", FIG$a$exposed, n_exp)
  lab_unx <- sprintf("%s (n = %d)", FIG$a$unexposed, n_unx)

  p_a <- panel_a_tbl %>%
    select(molecule, !!lab_exp := rate_exp, !!lab_unx := rate_unx, sig) %>%
    pivot_longer(c(!!lab_exp, !!lab_unx), names_to = "grp", values_to = "rate") %>%
    mutate(grp = factor(grp, levels = c(lab_exp, lab_unx)),
           molecule = factor(molecule, levels = rev(MOL_ORDER)),
           lab = ifelse(grp == lab_exp, sprintf("%.0f%% %s", rate, sig),
                        sprintf("%.0f%%", rate))) %>%
    ggplot(aes(rate, molecule, fill = grp)) +
    geom_col(position = position_dodge(width = 0.78), width = 0.72) +
    geom_text(aes(label = lab, colour = grp), position = position_dodge(width = 0.78),
              hjust = -0.1, size = 2.3, fontface = "bold", show.legend = FALSE) +
    scale_fill_manual(values = setNames(c(COL_EXPOSED, COL_UNEXPOSED), c(lab_exp, lab_unx)),
                      name = NULL) +
    scale_colour_manual(values = setNames(c(COL_EXPOSED, COL_UNEXPOSED), c(lab_exp, lab_unx))) +
    scale_x_continuous(labels = label_percent(scale = 1), expand = expansion(c(0, 0.18))) +
    labs(title = FIG$a$title, subtitle = FIG$a$subtitle, x = FIG$a$x, y = FIG$a$y,
         tag = FIG$a$tag) +
    theme_pub() + theme(legend.position = c(0.72, 0.18),
                        legend.background = element_rect(fill = "white", colour = NA))
  save_panel(p_a, "panel_a_detection_rate", SIZE$a)
}

## ---- Panel b ----------------------------------------------------------------
if (want("b")) {
  dia <- filter(dat, diabetes %in% c("YES", "NO"))
  n_dia_y <- n_distinct(dia$id[dia$diabetes == "YES"])
  n_dia_n <- n_distinct(dia$id[dia$diabetes == "NO"])

  panel_b_tbl <- dia %>% group_by(molecule) %>%
    summarise(yes = sum(detected & diabetes == "YES"),
              no  = sum(detected & diabetes == "NO"), .groups = "drop") %>%
    rowwise() %>%
    mutate(res = list(or_fisher(yes, n_dia_y - yes, no, n_dia_n - no)),
           OR = res$or, P = res$p, haldane = res$haldane) %>%
    ungroup() %>% select(-res) %>% arrange(match(molecule, MOL_ORDER))
  write_csv(panel_b_tbl, file.path(OUT_DIR, "panel_b_odds_ratios.csv"))

  cn <- c(sprintf("Yes\n(n = %d)", n_dia_y), sprintf("No\n(n = %d)", n_dia_n),
          FIG$b$cols[3], FIG$b$cols[4])
  p_b <- panel_b_tbl %>%
    transmute(molecule, !!cn[1] := as.character(yes), !!cn[2] := as.character(no),
              !!cn[3] := sprintf("%.1f", OR),
              !!cn[4] := format(P, digits = 2, scientific = TRUE), orv = OR) %>%
    pivot_longer(all_of(cn), names_to = "col", values_to = "txt") %>%
    mutate(col = factor(col, levels = cn),
           molecule = factor(molecule, levels = rev(MOL_ORDER)),
           fillv = ifelse(col == cn[3], log10(orv), NA_real_)) %>%
    ggplot(aes(col, molecule)) +
    geom_tile(aes(fill = fillv), colour = "white", linewidth = 0.6) +
    geom_text(aes(label = txt), size = 2.4) +
    scale_fill_gradient(low = OR_LOW, high = OR_HIGH, na.value = "white",
                        name = FIG$b$fill, labels = function(x) round(10^x),
                        guide = guide_colourbar(title.position = "top", barheight = 0.3,
                                                barwidth = 5)) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(title = FIG$b$title, subtitle = FIG$b$subtitle, x = NULL, y = NULL, tag = FIG$b$tag) +
    theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x.top = element_text(face = "bold", size = BASE_SIZE - 1))
  save_panel(p_b, "panel_b_odds_ratio_table", SIZE$b)
  log_line(sprintf("  b  diabetes yes n = %d, no n = %d%s", n_dia_y, n_dia_n,
                   if (any(panel_b_tbl$haldane))
                     sprintf("; Haldane correction for %d conjugates with a zero cell",
                             sum(panel_b_tbl$haldane)) else ""))
}

## ---- Panel c ----------------------------------------------------------------
if (want("c")) {
  hit_ids <- dat %>% group_by(id) %>% filter(any(detected)) %>% pull(id) %>% unique()
  panel_c_tbl <- dat %>% filter(id %in% hit_ids) %>%
    mutate(sample = str_remove(str_remove(id, "(?i)^sample[_\\.]?x?"), "^_"),
           value  = ifelse(detected, log(area + 1), NA_real_)) %>%
    select(sample, molecule, area, value)
  write_csv(panel_c_tbl, file.path(OUT_DIR, "panel_c_sample_heatmap.csv"))

  p_c <- panel_c_tbl %>%
    mutate(molecule = factor(molecule, levels = rev(MOL_ORDER))) %>%
    ggplot(aes(sample, molecule, fill = value)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    scale_fill_gradient(low = AREA_LOW, high = AREA_HIGH, na.value = "white",
                        name = FIG$c$fill,
                        guide = guide_colourbar(title.position = "top", barheight = 0.3,
                                                barwidth = 5)) +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
    labs(title = FIG$c$title, subtitle = FIG$c$subtitle, x = FIG$c$x, y = FIG$c$y,
         tag = FIG$c$tag) +
    theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = BASE_SIZE - 3))
  save_panel(p_c, "panel_c_sample_heatmap", SIZE$c)
  log_line(sprintf("  c  %d samples with at least one detection, %d detections total",
                   length(hit_ids), sum(dat$detected)))
}

## ---- Panel d ----------------------------------------------------------------
if (want("d")) {
  risk_long <- dat %>% filter(detected) %>%
    transmute(id, molecule = as.character(molecule),
              `Diabetes`           = diabetes == "YES",
              `HbA1c > 6.5`        = !is.na(hba1c) & hba1c > HBA1C_CUTOFF,
              `High blood sugar`   = high_sugar == "1",
              `Abdominal obesity`  = abd_obesity == "1",
              `Low HDL`            = low_hdl == "1",
              `High triglycerides` = high_trig == "1",
              `Insulin`            = insulin,
              `BMI`                = !is.na(bmi) & bmi >= BMI_CUTOFF) %>%
    pivot_longer(-c(id, molecule), names_to = "risk", values_to = "pos") %>%
    mutate(pos = !is.na(pos) & pos) %>% filter(pos) %>%
    count(molecule, risk, name = "n")
  write_csv(risk_long, file.path(OUT_DIR, "panel_d_chord_links.csv"))
  log_line(sprintf("  d  %d conjugate-risk links, total weight %d",
                   nrow(risk_long), sum(risk_long$n)))

  if (requireNamespace("circlize", quietly = TRUE)) {
    mol_order <- panel_a_tbl %>% filter(molecule %in% risk_long$molecule) %>%
      arrange(desc(rate_exp)) %>% pull(molecule) %>% as.character()
    grid_col <- c(setNames(rep("#BDBDBD", length(mol_order)), mol_order), RISK)
    mat <- risk_long %>%
      mutate(molecule = factor(molecule, levels = mol_order)) %>%
      pivot_wider(names_from = risk, values_from = n, values_fill = 0) %>%
      arrange(molecule) %>%
      select(molecule, all_of(names(RISK)))        # arcs follow the RISK order, not alphabetical
    m <- as.matrix(mat[, -1]); rownames(m) <- as.character(mat$molecule)

    for (dev in DEVICES) {
      f <- file.path(OUT_DIR, paste0("panel_d_chord.", dev))
      if (dev == "svg")      svglite::svglite(f, width = SIZE$d[1], height = SIZE$d[2])
      else if (dev == "pdf") cairo_pdf(f, width = SIZE$d[1], height = SIZE$d[2])
      else                   png(f, width = SIZE$d[1], height = SIZE$d[2], units = "in",
                                 res = PNG_DPI, bg = "white")
      circlize::circos.clear()
      circlize::circos.par(start.degree = 90, gap.degree = 1.2,
                           track.margin = c(0.01, 0.01), points.overflow.warning = FALSE)
      circlize::chordDiagram(m, grid.col = grid_col, transparency = 0.35, directional = 0,
                             annotationTrack = "grid",
                             preAllocateTracks = list(track.height = 0.22))
      circlize::circos.track(track.index = 1, panel.fun = function(x, y) {
        s <- circlize::get.cell.meta.data("sector.index")
        circlize::circos.text(circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1], s,
                              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                              cex = 0.5, col = if (s %in% names(RISK)) RISK[[s]] else "black")
      }, bg.border = NA)
      title(main = FIG$d$title, cex.main = 0.9, adj = 0, line = 0)
      circlize::circos.clear(); dev.off()
      message("  wrote panel_d_chord.", dev)
    }
  } else {
    log_line("  d  SKIPPED: package 'circlize' is not installed. Run ",
             "install.packages(\"circlize\") and re-run; the link table ",
             "panel_d_chord_links.csv is already written.")
  }
}

## =============================================================================
##  4. REPOSITORY DATA  (panel e)
## =============================================================================
if (want("e")) {
  repo <- read_tsv(REPO_FILE, show_col_types = FALSE, progress = FALSE)
  stopifnot(all(c("NAME_clean", "NCBITaxonomy", "filename", "USI") %in% names(repo)))

  if (PRE_FILTER) {
    n0 <- nrow(repo)
    repo <- filter(repo, cosine_raw >= MIN_COSINE, `Matching Peaks` >= MIN_MATCH_PEAKS)
    log_line(sprintf("  e  match-quality filter removed %d of %d rows", n0 - nrow(repo), n0))
  }

  COMPOUND_LOOKUP <- setNames(COMPOUND_ORDER, normalise_key(COMPOUND_ORDER))
  rdat <- repo %>%
    mutate(key = normalise_key(NAME_clean),
           compound = dplyr::coalesce(COMPOUND_ALIASES[key], COMPOUND_LOOKUP[key]),
           host = classify_host(NCBITaxonomy,
                                if ("NCBIKingdom" %in% names(repo)) NCBIKingdom else NA,
                                if ("NCBIFamily"  %in% names(repo)) NCBIFamily  else NA,
                                if ("NCBISpecies" %in% names(repo)) NCBISpecies else NA))

  ## Fail loudly rather than silently dropping a compound from the panel.
  unmapped <- sort(unique(rdat$NAME_clean[is.na(rdat$compound)]))
  if (length(unmapped))
    stop("NAME_clean values missing from COMPOUND_ORDER / COMPOUND_ALIASES:\n  ",
         paste(unmapped, collapse = "\n  "))

  n_unclassified <- sum(is.na(rdat$host))
  if (n_unclassified)
    log_line(sprintf("  e  %d of %d rows have no assignable host class%s",
                     n_unclassified, nrow(rdat),
                     if (DROP_UNCLASSIFIED) " and are dropped" else " and are kept as NA"))
  if (DROP_UNCLASSIFIED) rdat <- filter(rdat, !is.na(host))

  counts <- rdat %>% group_by(host, compound) %>%
    summarise(count = switch(COUNT_UNIT,
                             spectra = dplyr::n(),
                             files   = dplyr::n_distinct(filename),
                             usi     = dplyr::n_distinct(USI),
                             stop("COUNT_UNIT must be spectra, files or usi")),
              .groups = "drop")

  ## Complete the grid so a searched-but-empty host class still draws as a row.
  grid <- expand_grid(host = TAXON_ORDER, compound = COMPOUND_ORDER) %>%
    left_join(counts, by = c("host", "compound")) %>%
    mutate(host = factor(host, levels = rev(TAXON_ORDER)),
           compound = factor(compound, levels = COMPOUND_ORDER),
           log_count = ifelse(is.na(count), NA_real_, log(count + 1)),
           rt = ifelse(!is.na(count) & compound %in% RT_VALIDATED, "yes", NA_character_))
  write_csv(arrange(grid, host, compound) %>%
              select(host, compound, count, log_count, retention_time_validated = rt),
            file.path(OUT_DIR, "panel_e_occurrence_counts.csv"))

  p_e <- ggplot(grid, aes(compound, host)) +
    geom_tile(aes(fill = log_count), colour = "white", linewidth = 0.35) +
    geom_point(aes(shape = rt), size = 1.5, stroke = 0.45, colour = "#333333", na.rm = TRUE) +
    scale_fill_gradient2(name = FIG$e$fill, low = FILL_LOW, mid = FILL_MID, high = FILL_HIGH,
                         midpoint = FILL_MIDPOINT, limits = FILL_LIMITS,
                         breaks = pretty_breaks(3), na.value = FILL_NA,
                         guide = guide_colourbar(title.position = "top", order = 1,
                                                 ticks.colour = "white", frame.colour = NA)) +
    scale_shape_manual(name = NULL, values = c(yes = 8), labels = FIG$e$star,
                       na.translate = FALSE,
                       guide = guide_legend(order = 2, override.aes = list(size = 2))) +
    scale_x_discrete(labels = function(x) paste0("Metformin-", x), expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    coord_fixed(ratio = 1) +
    labs(title = FIG$e$title, subtitle = FIG$e$subtitle, tag = FIG$e$tag) +
    theme_grid()
  save_panel(p_e, "panel_e_occurrence_heatmap", SIZE$e)
  log_line(sprintf("  e  %d matches over %d conjugates; count unit = %s; host classes with data: %s",
                   nrow(rdat), n_distinct(rdat$compound), COUNT_UNIT,
                   paste(sort(unique(as.character(counts$host))), collapse = ", ")))
  log_line(sprintf("  e  retention-time confirmed: %d of %d conjugates",
                   sum(COMPOUND_ORDER %in% RT_VALIDATED), length(COMPOUND_ORDER)))
}

## =============================================================================
##  5. ORGANISED DATA TABLES  (one formatted sheet per panel, plus a README)
## =============================================================================
##  Every number plotted in panels a-e, laid out for a reader rather than for a
##  parser: block banners over related columns, units in the headers, frozen
##  panes and filters. The flat .csv files written above stay as the
##  machine-readable copy; this workbook is the one to attach as source data.

if (requireNamespace("openxlsx", quietly = TRUE)) {
  library(openxlsx)
  XL <- file.path(OUT_DIR, "Figure_source_data.xlsx")
  wb <- createWorkbook()
  addWorksheet(wb, "README", tabColour = "#555555")   # created first so it sorts to the front

  BLOCK_COL <- c("#3A6EA5", "#B5202C", "#4E8FD1", "#1F7A5C", "#6B4C9A", "#B06A2C")
  st_title <- createStyle(fontName = "Arial", fontSize = 12, textDecoration = "bold")
  st_sub   <- createStyle(fontName = "Arial", fontSize = 9, textDecoration = "italic", fontColour = "#666666")
  st_hdr   <- createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold",
                          fgFill = "#EDEDED", border = "TopBottomLeftRight",
                          borderColour = "#BFBFBF", wrapText = TRUE,
                          halign = "center", valign = "bottom")
  st_body  <- createStyle(fontName = "Arial", fontSize = 10, border = "TopBottomLeftRight",
                          borderColour = "#D9D9D9")
  st_key   <- createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold",
                          border = "TopBottomLeftRight", borderColour = "#D9D9D9")

  ## blocks: list(list(name = "...", n = <how many columns it spans>, col = <index>))
  add_sheet <- function(sheet, title, subtitle, blocks, df, widths,
                        key_cols = 1, tab_col = "#3A6EA5") {
    addWorksheet(wb, sheet, tabColour = tab_col)
    writeData(wb, sheet, title, startRow = 1); addStyle(wb, sheet, st_title, 1, 1)
    writeData(wb, sheet, subtitle, startRow = 2); addStyle(wb, sheet, st_sub, 2, 1)
    mergeCells(wb, sheet, cols = 1:max(8, ncol(df)), rows = 2)
    j <- 1
    for (b in blocks) {
      writeData(wb, sheet, b$name, startRow = 4, startCol = j)
      addStyle(wb, sheet, createStyle(fontName = "Arial", fontSize = 10,
                                      textDecoration = "bold", fontColour = "#FFFFFF",
                                      fgFill = BLOCK_COL[b$col], halign = "center",
                                      border = "TopBottomLeftRight",
                                      borderColour = "#BFBFBF"),
               rows = 4, cols = j:(j + b$n - 1), gridExpand = TRUE)
      if (b$n > 1) mergeCells(wb, sheet, cols = j:(j + b$n - 1), rows = 4)
      j <- j + b$n
    }
    writeData(wb, sheet, df, startRow = 5, headerStyle = st_hdr, borders = "all",
              borderColour = "#D9D9D9", keepNA = FALSE)
    addStyle(wb, sheet, st_body, rows = 6:(5 + nrow(df)), cols = 1:ncol(df),
             gridExpand = TRUE, stack = TRUE)
    if (key_cols > 0)
      addStyle(wb, sheet, st_key, rows = 6:(5 + nrow(df)), cols = 1:key_cols,
               gridExpand = TRUE, stack = TRUE)
    setColWidths(wb, sheet, cols = seq_along(widths), widths = widths)
    setRowHeights(wb, sheet, rows = 5, heights = 42)
    freezePane(wb, sheet, firstActiveRow = 6, firstActiveCol = key_cols + 1)
    addFilter(wb, sheet, rows = 5, cols = 1:ncol(df))
  }

  readme <- data.frame(Sheet = character(), Contents = character())
  add_readme <- function(sh, txt) readme <<- rbind(readme, data.frame(Sheet = sh, Contents = txt))

  ## ---- Fig a ----------------------------------------------------------------
  if (want("a")) {
    d <- panel_a_tbl %>% arrange(match(molecule, MOL_ORDER)) %>%
      transmute(Conjugate = as.character(molecule),
                `Detected` = det_exp, `Not detected` = n_exp - det_exp,
                `Detection rate (%)` = round(rate_exp, 1),
                `Detected ` = det_unx, `Not detected ` = n_unx - det_unx,
                `Detection rate (%) ` = round(rate_unx, 1),
                `Fisher P` = P, Significance = sig)
    add_sheet("Fig a  detection rate",
              "Figure panel a - conjugate detection by metformin exposure",
              sprintf(paste("One row per conjugate. Metformin-exposed n = %d, unexposed n = %d,",
                            "on %d samples. Detected means Area_sub_corrected > 0. P is a",
                            "two-sided Fisher exact test of the 2x2, uncorrected."),
                      n_exp, n_unx, n_exp + n_unx),
              list(list(name = "Conjugate", n = 1, col = 1),
                   list(name = sprintf("Metformin-exposed (n = %d)", n_exp), n = 3, col = 2),
                   list(name = sprintf("Unexposed (n = %d)", n_unx), n = 3, col = 3),
                   list(name = "Exposed vs unexposed", n = 2, col = 4)),
              d, c(34, 11, 13, 14, 11, 13, 14, 13, 13), tab_col = "#B5202C")
    add_readme("Fig a  detection rate",
               "Detection counts and rates per conjugate in metformin-exposed vs unexposed samples, with the Fisher exact P value plotted as asterisks.")
  }

  ## ---- Fig b ----------------------------------------------------------------
  if (want("b")) {
    d <- panel_b_tbl %>%
      transmute(Conjugate = as.character(molecule),
                `Detected` = yes, `Not detected` = n_dia_y - yes,
                `Detected ` = no, `Not detected ` = n_dia_n - no,
                `Odds ratio` = round(OR, 1),
                `Odds ratio method` = ifelse(haldane, "Haldane-Anscombe (+0.5)", "plain ad/bc"),
                `Fisher P` = P, Significance = stars(P))
    add_sheet("Fig b  diabetes odds ratio",
              "Figure panel b - conjugate detection by reported diabetes",
              sprintf(paste("One row per conjugate. Diabetes reported yes n = %d, no n = %d.",
                            "One sample with diabetes = NA is excluded. The odds ratio is the",
                            "plain ad/bc, except where a cell is zero, when the",
                            "Haldane-Anscombe correction is used. P is two-sided Fisher, uncorrected."),
                      n_dia_y, n_dia_n),
              list(list(name = "Conjugate", n = 1, col = 1),
                   list(name = sprintf("Diabetes reported: Yes (n = %d)", n_dia_y), n = 2, col = 2),
                   list(name = sprintf("Diabetes reported: No (n = %d)", n_dia_n), n = 2, col = 3),
                   list(name = "Association", n = 4, col = 4)),
              d, c(34, 11, 13, 11, 13, 12, 24, 13, 13), tab_col = "#2F5FA8")
    add_readme("Fig b  diabetes odds ratio",
               "2x2 counts, odds ratio and Fisher exact P for each conjugate against reported diabetes. The odds ratio column is the one shaded blue in the figure.")
  }

  ## ---- Fig c ----------------------------------------------------------------
  if (want("c")) {
    w <- panel_c_tbl %>%
      mutate(molecule = factor(as.character(molecule), levels = MOL_ORDER)) %>%
      select(molecule, sample, value) %>%
      pivot_wider(names_from = sample, values_from = value) %>%
      arrange(molecule)
    d <- bind_cols(
      tibble(Conjugate = as.character(w$molecule),
             `Samples detected in` = rowSums(!is.na(w[, -1]))),
      w[, -1])
    add_sheet("Fig c  per-sample areas",
              "Figure panel c - log(peak area + 1) per sample",
              sprintf(paste("Rows are conjugates in figure order; columns are the %d samples with",
                            "at least one detection. Values are the natural log of the",
                            "background-subtracted peak area plus one. An empty cell means not",
                            "detected in that sample, not a missing measurement."),
                      ncol(w) - 1),
              list(list(name = "Conjugate", n = 2, col = 1),
                   list(name = "Sample ID", n = ncol(d) - 2, col = 3)),
              d, c(34, 17, rep(11, ncol(d) - 2)), key_cols = 2, tab_col = "#C0392B")
    add_readme("Fig c  per-sample areas",
               "Conjugate x sample matrix of log(peak area + 1). Empty = not detected.")
  }

  ## ---- Fig d ----------------------------------------------------------------
  if (want("d")) {
    w <- risk_long %>%
      mutate(molecule = factor(molecule, levels = MOL_ORDER)) %>%
      pivot_wider(names_from = risk, values_from = n, values_fill = 0) %>%
      arrange(molecule) %>%
      select(molecule, all_of(names(RISK)))        # keep the figure's risk-factor order
    d <- bind_cols(
      tibble(Conjugate = as.character(w$molecule),
             `Total links` = rowSums(w[, -1])), w[, -1])
    add_sheet("Fig d  risk factor links",
              "Figure panel d - conjugate to metabolic risk factor",
              sprintf(paste("Each value is the number of samples in which that conjugate was",
                            "detected AND that risk factor was positive; it is the ribbon width",
                            "in the chord diagram. HbA1c positive is > %.1f%%, BMI positive is",
                            ">= %d kg/m2, Insulin is any insulin medication recorded, and the",
                            "remaining factors are the cohort's own metabolic-syndrome flags."),
                      HBA1C_CUTOFF, BMI_CUTOFF),
              list(list(name = "Conjugate", n = 2, col = 1),
                   list(name = "Risk factor", n = ncol(d) - 2, col = 5)),
              d, c(34, 12, rep(15, ncol(d) - 2)), key_cols = 2, tab_col = "#0E7C7B")
    add_readme("Fig d  risk factor links",
               "Conjugate x risk-factor matrix of co-positive sample counts; the ribbon widths of the chord diagram.")
  }

  ## ---- Fig e ----------------------------------------------------------------
  if (want("e")) {
    ## One row per conjugate keeps the counts numeric and gives the retention-time
    ## flag a column of its own, rather than a text row that would coerce the matrix.
    d <- grid %>% select(host, compound, count) %>%
      mutate(host = factor(as.character(host), levels = TAXON_ORDER)) %>%
      pivot_wider(names_from = host, values_from = count) %>%
      arrange(match(as.character(compound), COMPOUND_ORDER)) %>%
      transmute(Conjugate = paste0("Metformin-", compound),
                across(all_of(TAXON_ORDER)),
                `Total matches` = rowSums(across(all_of(TAXON_ORDER)), na.rm = TRUE),
                `Retention-time validated` =
                  ifelse(as.character(compound) %in% RT_VALIDATED, "yes", ""))
    add_sheet("Fig e  repository counts",
              "Figure panel e - library matches in public repositories by host class",
              sprintf(paste("One row per conjugate. Each value is the number of library matches",
                            "(%s) in that host class; an empty cell means the class was searched",
                            "and nothing was found. The figure plots log(count + 1). The last",
                            "column marks conjugates confirmed against an authentic standard by",
                            "retention time, which is recorded in the script, not in either data file."),
                      COUNT_UNIT),
              list(list(name = "Conjugate", n = 1, col = 1),
                   list(name = "Library matches by host class", n = length(TAXON_ORDER) + 1, col = 6),
                   list(name = "Confirmation", n = 1, col = 4)),
              d, c(36, rep(14, length(TAXON_ORDER)), 14, 24), key_cols = 1, tab_col = "#4A77B4")
    add_readme("Fig e  repository counts",
               "Host class x conjugate matrix of library-match counts, plus the retention-time confirmation row.")
  }

  ## ---- README ---------------------------------------------------------------
  writeData(wb, "README", "Source data for figure panels a to e", startRow = 1)
  addStyle(wb, "README", st_title, 1, 1)
  info <- data.frame(Item = c("Cohort data", "Repository data", "Sample", "Detected",
                              "Odds ratio", "P value", "Panel e fill", "Generated"),
                     Detail = c(
                       COHORT_FILE, REPO_FILE,
                       "One biological sample, keyed on ID. Repeat injections collapse to the maximum peak area.",
                       "Area_sub_corrected > 0.",
                       "ad/bc, with the Haldane-Anscombe correction (+0.5 to every cell) only where a cell is zero.",
                       "Fisher exact test, two-sided, raw and uncorrected for multiple comparisons.",
                       "Natural log of the match count plus one.",
                       format(Sys.time(), "%Y-%m-%d %H:%M")))
  writeData(wb, "README", info, startRow = 3, headerStyle = st_hdr)
  writeData(wb, "README", readme, startRow = 3 + nrow(info) + 2, headerStyle = st_hdr)
  setColWidths(wb, "README", cols = 1:2, widths = c(26, 118))

  saveWorkbook(wb, XL, overwrite = TRUE)

  ## openxlsx declares a drawing and a vmlDrawing relationship on every sheet but
  ## writes neither file when the sheet has no shapes. Excel and LibreOffice
  ## ignore the dangling reference; stricter readers (openpyxl, some R and
  ## Python parsers a reviewer may use) refuse the file. Drop the references
  ## whose target is genuinely absent.
  repair_xlsx <- function(path) {
    tmp <- file.path(tempdir(), paste0("xlsxrepair_", as.integer(Sys.time())))
    unlink(tmp, recursive = TRUE); dir.create(tmp, recursive = TRUE)
    utils::unzip(path, exdir = tmp)
    present <- list.files(file.path(tmp, "xl", "drawings"))
    dropped <- 0L
    for (f in list.files(file.path(tmp, "xl", "worksheets", "_rels"), full.names = TRUE)) {
      x <- paste(readLines(f, warn = FALSE), collapse = "")
      for (tag in regmatches(x, gregexpr("<Relationship[^>]*/>", x))[[1]]) {
        tgt <- sub('.*Target="\\.\\./([^"]+)".*', "\\1", tag)
        if (grepl("^drawings/", tgt) && !(basename(tgt) %in% present)) {
          x <- sub(tag, "", x, fixed = TRUE); dropped <- dropped + 1L
        }
      }
      writeLines(x, f)
    }
    for (f in list.files(file.path(tmp, "xl", "worksheets"), pattern = "\\.xml$",
                         full.names = TRUE)) {
      x <- paste(readLines(f, warn = FALSE), collapse = "")
      x <- gsub("<drawing [^>]*/>|<legacyDrawing [^>]*/>", "", x)
      writeLines(x, f)
    }
    owd <- setwd(tmp); on.exit(setwd(owd), add = TRUE)
    out <- file.path(owd, basename(path))
    if (!startsWith(path, "/")) out <- file.path(owd, path) else out <- path
    unlink(out)
    utils::zip(out, list.files(".", recursive = TRUE, all.files = TRUE, no.. = TRUE),
               flags = "-q -X")
    dropped
  }
  dropped <- tryCatch(repair_xlsx(normalizePath(XL)), error = function(e) NA_integer_)

  log_line(sprintf("  organised source-data workbook: %s (%d sheets%s)", basename(XL),
                   length(names(wb)),
                   if (is.na(dropped)) "; drawing-reference repair skipped"
                   else sprintf("; %d dangling drawing references removed", dropped)))
} else {
  log_line("  openxlsx not installed; the .csv tables were written but the formatted ",
           "workbook was skipped. Run install.packages(\"openxlsx\").")
}

## =============================================================================
##  6. STATISTICS LOG
## =============================================================================
if (want("a")) {
  log_line(sprintf("  a  exposed n = %d, unexposed n = %d", n_exp, n_unx))
  for (i in seq_len(nrow(panel_a_tbl))) with(panel_a_tbl[i, ],
    log_line(sprintf("     %-34s exposed %2d/%d = %4.1f%%   unexposed %2d/%d = %4.1f%%   P = %.2e %s",
                     molecule, det_exp, n_exp, rate_exp, det_unx, n_unx, rate_unx, P, sig)))
}
writeLines(c(paste("Generated", format(Sys.time(), "%Y-%m-%d %H:%M")),
             paste("Cohort source:    ", COHORT_FILE),
             paste("Repository source:", REPO_FILE),
             paste("Panels built:     ", paste(PANELS, collapse = ", ")), "", LOG),
           file.path(OUT_DIR, "panel_statistics.txt"))
message("done -> ", normalizePath(OUT_DIR))
