# see expected toy call

    Code
      partial_decouple(show_toy_call = TRUE, include_time = FALSE) %>% {
        TRUE
      }
    Message <message>
      run_udt(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
      run_mdt(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
      run_aucell(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
      run_wmean(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
      run_wsum(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
      run_ulm(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
      run_viper(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
      run_gsva(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
      run_ora(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
      run_fgsea(mat = mat, network = dorothea_genesets, .source = tf, .target = target,
      *   force_ties = FALSE)
    Warning <simpleWarning>
      
            FGSEA: Ties were detected, NAs will be returned.
            To force ties use force_ties = T, but results might not be reproducible.
      
            FGSEA: Ties were detected, NAs will be returned.
            To force ties use force_ties = T, but results might not be reproducible.
      
            FGSEA: Ties were detected, NAs will be returned.
            To force ties use force_ties = T, but results might not be reproducible.
      
            FGSEA: Ties were detected, NAs will be returned.
            To force ties use force_ties = T, but results might not be reproducible.
    Message <message>
      run_mlm(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
    Output
      [1] TRUE

