# Title     : HTML Functions
# Objective : Functions for creating pretty HTML links
# Created by: mumichae
# Created on: 6/21/21

build_link_list <- function(file_paths, captions=NULL) {
  if (is.null(captions)) {
    captions <- file_paths
  }
  file_link <- paste0('\n* [', captions , '](', file_paths,
                      '){target="_blank"}\n', collapse = ' ')
  file_link
}

display_text <- function(caption='', links) {
  captions <- ifelse(
    caption != '',
    paste0('**', caption, '**', names(links)),
    caption
  )
  paste0(captions, '\n', links, collapse = '\n')
}
