# see expected toy call

    Code
      partial_decouple(show_toy_call = TRUE, include_time = FALSE)
    Message <simpleMessage>
      run_scira(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
    Message <simpleMessage>
      run_pscira(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
    Message <simpleMessage>
      run_mean(mat = mat, network = dorothea_genesets, .source = tf, .target = target,  .likelihood = NULL)
    Message <simpleMessage>
      run_viper(mat = mat, network = dorothea_genesets, .source = tf, .target = target,  verbose = FALSE)
    Message <simpleMessage>
      run_gsva(mat = mat, network = dorothea_genesets, .source = tf, .target = target,  verbose = FALSE)
    Message <simpleMessage>
      run_ora(mat = mat, network = dorothea_genesets, .source = tf, .target = target)
    Output
      # A tibble: 132 x 11
         run_id statistic tf    condition score p_value estimate conf.low conf.high
         <chr>  <chr>     <chr> <chr>     <dbl>   <dbl>    <dbl>    <dbl>     <dbl>
       1 1      scira     FOXO4 GSM27533~ 2.06       NA       NA       NA        NA
       2 1      scira     FOXO4 GSM27533~ 2.05       NA       NA       NA        NA
       3 1      scira     FOXO4 GSM27533~ 2.32       NA       NA       NA        NA
       4 1      scira     FOXO4 GSM27533~ 2.31       NA       NA       NA        NA
       5 1      scira     NFIC  GSM27533~ 0.903      NA       NA       NA        NA
       6 1      scira     NFIC  GSM27533~ 0.978      NA       NA       NA        NA
       7 1      scira     NFIC  GSM27533~ 0.679      NA       NA       NA        NA
       8 1      scira     NFIC  GSM27533~ 0.748      NA       NA       NA        NA
       9 1      scira     SMAD3 GSM27533~ 0.548      NA       NA       NA        NA
      10 1      scira     SMAD3 GSM27533~ 0.786      NA       NA       NA        NA
      # ... with 122 more rows, and 2 more variables: method <chr>, alternative <chr>

