# see expected toy call

    Code
      partial_decouple(show_toy_call = TRUE, include_time = FALSE) %>% {
        TRUE
      }
    Message <message>
      run_udt(mat = mat, network = dorothea_genesets, .source = tf, .target = target,
      *   minsize = 5)
      run_mdt(mat = mat, network = dorothea_genesets, .source = tf, .target = target,
      *   trees = 1000, minsize = 5)
      run_aucell(mat = mat, network = dorothea_genesets, .source = tf, .target = target,
      *   nproc = 1, minsize = 5)
      run_wmean(mat = mat, network = dorothea_genesets, .source = tf, .target = target,
      *   minsize = 5)
      run_wsum(mat = mat, network = dorothea_genesets, .source = tf, .target = target,
      *   minsize = 5)
      run_ulm(mat = mat, network = dorothea_genesets, .source = tf, .target = target,
      *   minsize = 5)
      run_viper(mat = mat, network = dorothea_genesets, .source = tf, .target = target,
      *   minsize = 5)
      run_gsva(mat = mat, network = dorothea_genesets, .source = tf, .target = target,
      *   minsize = 5)
      run_ora(mat = mat, network = dorothea_genesets, .source = tf, .target = target,
      *   n_up = 300, n_bottom = 300, minsize = 5)
      run_fgsea(mat = mat, network = dorothea_genesets, .source = tf, .target = target,
      *   minsize = 5)
      run_mlm(mat = mat, network = dorothea_genesets, .source = tf, .target = target,
      *   minsize = 5)
    Output
      [1] TRUE

