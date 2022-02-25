library(decoupleR)
library(purrr)

# Base directories definition --------------------------------------------------

# Inputs
input_dir <- system.file("testdata", "inputs", package = "decoupleR")

# Outputs
output_dir <- file.path("inst", "testdata", "outputs")

# Specific directories creation ------------------------------------------------
# Here you need to extend the vector with the name of your new statistic.
available_statistics = c(
    'udt',
    'mdt',
    'aucell',
    'wmean',
    'wsum',
    'ulm',
    'mlm',
    'viper',
    'gsva',
    'ora',
    'fgsea'
)

out <- available_statistics %>%
    file.path(output_dir, .) %>%
    setNames(object = ., nm = basename(.)) %>%
    as.list()

sapply(out, dir.create, showWarnings = TRUE)

# Collect individual default outputs files to decouple() test.

out_default <- stringr::str_glue(
    "output-{available_statistics}_dorothea_default.rds"
) %>%
    map2(.x = out, .y = ., file.path)

decouple_dir <- file.path(output_dir, "decouple")
dir.create(decouple_dir, showWarnings = TRUE)

# Load data to generated outputs -----------------------------------------------

emat <- file.path(input_dir, "input-expr_matrix.rds") %>%
    readRDS()

dorothea_genesets <- file.path(input_dir, "input-dorothea_genesets.rds") %>%
    readRDS()

#----- run_udt() ---------------------------------------------------------------
run_udt(emat, dorothea_genesets, .source='tf') %>%
    saveRDS(out_default$udt)

#----- run_mdt() ---------------------------------------------------------------
run_mdt(emat, dorothea_genesets, .source='tf', trees=1000) %>%
    saveRDS(out_default$mdt)

#----- run_aucell() ------------------------------------------------------------
run_aucell(emat, dorothea_genesets, .source='tf', nproc=1) %>%
    saveRDS(out_default$aucell)

#----- run_wmean() -------------------------------------------------------------
run_wmean(emat, dorothea_genesets, .source='tf') %>%
    saveRDS(out_default$wmean)

#----- run_wsum() --------------------------------------------------------------
run_wsum(emat, dorothea_genesets, .source='tf') %>%
    saveRDS(out_default$wsum)

#----- run_ulm() ---------------------------------------------------------------
run_ulm(emat, dorothea_genesets, .source='tf') %>%
    saveRDS(out_default$ulm)

#----- run_mlm() ---------------------------------------------------------------
run_mlm(emat, dorothea_genesets, .source='tf') %>%
    saveRDS(out_default$mlm)

#----- run_viper() -------------------------------------------------------------
run_viper(emat, dorothea_genesets, .source='tf') %>%
    saveRDS(out_default$viper)

#----- run_gsva() --------------------------------------------------------------
run_gsva(emat, dorothea_genesets, .source='tf') %>%
    saveRDS(out_default$gsva)

#----- run_ora() ---------------------------------------------------------------
run_ora(emat, dorothea_genesets, .source='tf', n_up=300, n_bottom=300) %>%
    saveRDS(out_default$ora)

#----- run_fgsea() -------------------------------------------------------------
run_fgsea(emat, dorothea_genesets, .source='tf') %>%
    saveRDS(out_default$fgsea)

# decouple() --------------------------------------------------------------
# This section should be kept at the end of the file
# and should be executed every time a new statistic
# is added or any entry of the default models is modified.

map_dfr(out_default, readRDS) %>%
    dplyr::arrange(.data$statistic, .data$source, .data$condition) %>%
    saveRDS(file.path(decouple_dir, "output-decouple_dorothea_default.rds"))
