# see expected toy call

    Code
      partial_decouple(show_toy_call = TRUE, include_time = FALSE) %>% {
        TRUE
      }
    Message <simpleMessage>
      run_scira(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
      run_pscira(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
      run_mean(mat = mat, network = dorothea_genesets, .source = tf, .target = target,  .likelihood = NULL)
      run_viper(mat = mat, network = dorothea_genesets, .source = tf, .target = target,  verbose = FALSE)
      run_gsva(mat = mat, network = dorothea_genesets, .source = tf, .target = target,  verbose = FALSE)
      run_ora(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
    Output
      [1] TRUE

