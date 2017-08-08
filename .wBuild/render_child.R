# R functions
#
# author: daniel bader
#-------------------------------

#' render_child
#' 
#' Analog to knit_child(). Takes an R notebook and spins and knits it to an HTML document.
#' The HTML code is written to STDOUT.
#' 
render_child= function(file, quiet= TRUE, knit_fct= 'knit_child', return=FALSE, knit_flag=FALSE ){
	library(knitr)
	
	textin= readLines(file)
	textout= get(knit_fct)(
			text= spin(text= textin, knit=knit_flag, report=F, format='Rmd'), 
			quiet= quiet
	)
	if(return){
		return(textout)
	}else{
		writeLines(textout)
	}
	rm(list=c('textin', 'textout'))
}


#' rmarkdown_render_on_demand
#' 
#' Function to render a Rscript file on demand, w.r.t. a boolean if clause
#' 
rmarkdown_render_on_demand <- function(script, spin_script, quiet=TRUE){
	
	# generate prot report
	if(spin_script){
		render_child(file=script, quiet= quiet)
	}
}


#' enumerate_htmlwidgets
#' 
#' The input file is scanned for htmlwidget container. If present there IDs will be enumerated, 
#' starting by 1 to avoid ID conflicts, i.e. two container having the same ID.
#' 
enumerate_htmlwidgets= function(infile, outfile= infile, prefix='hw'){
	
	file_by_lines <- readLines(infile)
	counter=1
	
	for(i in 1:length(file_by_lines)){
		line = file_by_lines[i]
		line= sub('id=\"htmlwidget-[0-9]{1,4}', paste0('id=\"', prefix, counter), line)
		
		if(grepl('data-for=\"htmlwidget-[0-9]', line)){
			line= sub('data-for=\"htmlwidget-[0-9]{1,4}', paste0('data-for=\"', prefix, counter), line)
			counter= counter+1
		}
		file_by_lines[i] = line
	}
	writeLines(file_by_lines, outfile)
}


render_clean_html= function(rscript){
	render(rscript, output_format = 'html_document')
	enumerate_htmlwidgets(sub('.R','.html',rscript))
}
