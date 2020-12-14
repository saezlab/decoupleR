# see expected toy call

    Code
      partial_decouple(show_toy_call = TRUE) %>% select(-statistic_time)
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
    Output
      # A tibble: 112 x 6
         run_id statistic tf    condition  score p_value
         <chr>  <chr>     <chr> <chr>      <dbl>   <dbl>
       1 1      scira     FOXO4 GSM2753335 2.06       NA
       2 1      scira     FOXO4 GSM2753336 2.05       NA
       3 1      scira     FOXO4 GSM2753337 2.32       NA
       4 1      scira     FOXO4 GSM2753338 2.31       NA
       5 1      scira     NFIC  GSM2753335 0.903      NA
       6 1      scira     NFIC  GSM2753336 0.978      NA
       7 1      scira     NFIC  GSM2753337 0.679      NA
       8 1      scira     NFIC  GSM2753338 0.748      NA
       9 1      scira     SMAD3 GSM2753335 0.548      NA
      10 1      scira     SMAD3 GSM2753336 0.786      NA
      # ... with 102 more rows

