# see expected toy call

    Code
      partial_decouple(show_toy_call = TRUE, include_time = FALSE)
    Message <simpleMessage>
      run_scira(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
      run_pscira(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
      run_mean(mat = mat, network = dorothea_genesets, .source = tf, .target = target,  .likelihood = NULL)
      run_viper(mat = mat, network = dorothea_genesets, .source = tf, .target = target,  verbose = FALSE)
      run_gsva(mat = mat, network = dorothea_genesets, .source = tf, .target = target,  verbose = FALSE)
      run_ora(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
    Output
      # A tibble: 140 x 11
         run_id statistic tf    condition   score p_value estimate conf.low conf.high
         <chr>  <chr>     <chr> <chr>       <dbl>   <dbl>    <dbl>    <dbl>     <dbl>
       1 1      scira     FOXO4 GSM2753335  0.472      NA       NA       NA        NA
       2 1      scira     FOXO4 GSM2753336 -1.29       NA       NA       NA        NA
       3 1      scira     FOXO4 GSM2753337  1.02       NA       NA       NA        NA
       4 1      scira     FOXO4 GSM2753338 -0.251      NA       NA       NA        NA
       5 1      scira     NFIC  GSM2753335  2.24       NA       NA       NA        NA
       6 1      scira     NFIC  GSM2753336  1.76       NA       NA       NA        NA
       7 1      scira     NFIC  GSM2753337 -1.93       NA       NA       NA        NA
       8 1      scira     NFIC  GSM2753338 -2.21       NA       NA       NA        NA
       9 1      scira     SMAD3 GSM2753335 -0.495      NA       NA       NA        NA
      10 1      scira     SMAD3 GSM2753336 -0.234      NA       NA       NA        NA
      # ... with 130 more rows, and 2 more variables: method <chr>, alternative <chr>

