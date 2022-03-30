#' DoRothEA gene regulatory network.
#'
#' @description
#' Wrapper to access DoRothEA gene regulatory network. DoRothEA is a 
#' comprehensive resource containing a curated collection of transcription 
#' factors (TFs) and their target genes. Each interaction is weighted by its 
#' mode of regulation (either positive or negative) and by its confidence level
#' 
#' @param organism Which organism to use. Only human and mouse are available.
#' @param levels List of confidence levels to return. Goes from A to D, A 
#' being the most confident and D being the less.
#' @param weight_dict Dictionary of values to divide the mode of regulation 
#' (-1 or 1), one for each confidence level. Bigger values will generate 
#' weights close to zero.
#' 
#' @export
#' @examples 
#' dorothea <- get_dorothea(organism='human', levels=c('A', 'B', 'C'))
get_dorothea <- function(organism='human', levels=c('A', 'B', 'C'),
                         weight_dict = list('A'= 1, 'B'= 2, 'C'= 3, 'D'= 4)){
  
  # Process organism
  organism <- tolower(organism)
  if (!organism %in% c('human','mouse')){
    stop("organism can only be human or mouse.")
  }
  if(organism=='human'){
    organism <- 9606
  } else {
    organism <- 10090
  }
  
  # Get Dorothea
  do <- OmnipathR::import_dorothea_interactions(organism = organism,
                                                dorothea_levels = c('A','B','C','D'),
                                                genesymbols=T) %>%
    
    # Filter columns
    dplyr::select('source_genesymbol', 'target_genesymbol', 'is_stimulation', 'is_inhibition',
                  'consensus_direction', 'consensus_stimulation', 'consensus_inhibition',
                  'dorothea_level') %>%
    
    # Remove duplicates
    dplyr::distinct(.data$source_genesymbol, .data$dorothea_level, .data$target_genesymbol, .keep_all = TRUE) %>%
    
    # Get bets confidence if more than one
    dplyr::mutate(dorothea_level=unlist(map(.data$dorothea_level, function(lvl){
      stringr::str_split(lvl, ';')[[1]][[1]]
    }))) %>%
    
    # Define mor
    mutate(
      mor=ifelse(
        .data$is_stimulation & .data$is_inhibition,
        ifelse(.data$consensus_stimulation, 1, -1),
        ifelse(.data$is_stimulation, 1,ifelse(.data$is_inhibition, -1, 1))
      )
    ) %>%
    
    # Weight mor by confidence
    mutate(mor=.data$mor / unlist(map(.data$dorothea_level, function(lvl){weight_dict[[lvl]]}))) %>%
    
    # Filter columns
    dplyr::select('source_genesymbol', 'dorothea_level', 'target_genesymbol', 'mor')
  
  # Rename
  colnames(do) <- c('source', 'confidence', 'target', 'mor')
  
  # Filter by levels
  do <- do %>% dplyr::filter(.data$confidence %in% levels)
  
  return(do)
}


#' Shows available resources in Omnipath. For more information visit the 
#' official website for [Omnipath](https://omnipathdb.org/).
#' 
#' @export
#' @examples 
#' decoupleR::show_resources()
show_resources <- function(){
  return(OmnipathR::get_annotation_resources())
}


#' Wrapper to access resources inside Omnipath.
#' This wrapper allows to easily query different prior knowledge resources.
#' To check available resources run `decoupleR::show_resources()`. For more 
#' information visit the official website for [Omnipath](https://omnipathdb.org/).
#' 
#' @param name Name of the resource to query.
#' 
#' @export
#' @examples 
#' df <- decoupleR::get_resource('SIGNOR')
get_resource <- function(name){
  if (!name %in% show_resources()){
    stop(stringr::str_glue('{name} is not a valid resource. Please, run 
                         decoupleR::show_resources() to see the list of 
                         available resources.'))
  }
  df <- OmnipathR::import_omnipath_annotations(resources = name) %>%
    OmnipathR::pivot_annotations(.)
  return(df)
}


#' Pathway RespOnsive GENes for activity inference (PROGENy).
#' 
#' Wrapper to access PROGENy model gene weights. Each pathway is defined with a 
#' collection of target genes, each interaction has an associated p-value and 
#' weight. The top significant interactions per pathway are returned.
#' 
#' @param organism Which organism to use. Only human and mouse are available.
#' @param top Number of genes per pathway to return.
#' 
#' @importFrom utils head
#' 
#' @export
#' @examples 
#' progeny <- get_progeny(organism='human', top=500)
get_progeny <- function(organism='human', top=500){
  # Process organism
  organism <- tolower(organism)
  if (!organism %in% c('human','mouse')){
    stop("organism can only be human or mouse.")
  }
  p <- get_resource('PROGENy') %>%
    dplyr::distinct(.data$pathway, .data$genesymbol, .keep_all = TRUE) %>%
    dplyr::mutate(weight=as.double(.data$weight), p_value=as.double(.data$p_value)) %>%
    dplyr::select(.data$genesymbol, .data$p_value, .data$pathway, .data$weight) %>%
    dplyr::group_by(.data$pathway) %>%
    dplyr::group_split() %>%
    purrr::map(function(df){
      df %>%
        dplyr::arrange(.data$p_value) %>%
        head(top)
    }) %>%
    dplyr::bind_rows() %>%
    dplyr::select(.data$pathway, .data$genesymbol, .data$weight, .data$p_value)
  colnames(p) <- c('source', 'target', 'weight', 'p_value')
  
  if (organism=='mouse'){
    p$target <- stringr::str_to_sentence(p$target)
  }
  
  return(p)
}

