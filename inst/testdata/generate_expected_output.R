library(decoupleR)
library(purrr)

# Base directories definition --------------------------------------------------
input_dir <- file.path("inst", "testdata", "inputs")
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
    "output-{available_statistics}.rds"
) %>%
    map2(.x = out, .y = ., file.path)

decouple_dir <- file.path(output_dir, "decouple")
dir.create(decouple_dir, showWarnings = TRUE)

# Load data to generated outputs -----------------------------------------------
mat <- readRDS(file.path(input_dir, 'mat.rds'))
net <- readRDS(file.path(input_dir, 'net.rds'))

#----- run_udt() ---------------------------------------------------------------
run_udt(mat, net, minsize = 0) %>%
    saveRDS(out_default$udt)

#----- run_mdt() ---------------------------------------------------------------
run_mdt(mat, net, minsize=0, trees=1000) %>%
    saveRDS(out_default$mdt)

#----- run_aucell() ------------------------------------------------------------
run_aucell(mat, net, minsize=0, nproc=1, aucMaxRank=3) %>%
    saveRDS(out_default$aucell)

#----- run_wmean() -------------------------------------------------------------
run_wmean(mat, net, minsize=0) %>%
    saveRDS(out_default$wmean)

#----- run_wsum() --------------------------------------------------------------
run_wsum(mat, net, minsize=0) %>%
    saveRDS(out_default$wsum)

#----- run_ulm() ---------------------------------------------------------------
run_ulm(mat, net, minsize=0) %>%
    saveRDS(out_default$ulm)

#----- run_mlm() ---------------------------------------------------------------
run_mlm(mat, net, minsize=0) %>%
    saveRDS(out_default$mlm)

#----- run_viper() -------------------------------------------------------------
run_viper(mat, net, minsize=0) %>%
    saveRDS(out_default$viper)

#----- run_gsva() --------------------------------------------------------------
run_gsva(mat, net, minsize=0) %>%
    saveRDS(out_default$gsva)

#----- run_ora() ---------------------------------------------------------------
run_ora(mat, net, minsize=0, n_up=3, n_bottom=3) %>%
    saveRDS(out_default$ora)

#----- run_fgsea() -------------------------------------------------------------
run_fgsea(mat, net, minsize=0) %>%
    saveRDS(out_default$fgsea)

# decouple() --------------------------------------------------------------
# This section should be kept at the end of the file
# and should be executed every time a new statistic
# is added or any entry of the default models is modified.

map_dfr(out_default, readRDS) %>%
    dplyr::arrange(.data$statistic, .data$source, .data$condition) %>%
    saveRDS(file.path(decouple_dir, "output-decouple.rds"))
