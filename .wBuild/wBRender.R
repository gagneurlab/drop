
#require(knitrBootstrap)
require(knitr)
opts_knit$set(root.dir = getwd())
require(rmarkdown)
opts_chunk$set(echo=TRUE,message=FALSE,error=FALSE,warning=FALSE, cache=F)
options(width=100)

source(".wBuild/render_child.R")
source(".wBuild/rmarkdown_show_hide_table.R")

#copy dependency to intermediate dir
intermediates_dir = tempfile()
i = dir.create(file.path(dirname(intermediates_dir), basename(intermediates_dir)), showWarnings = FALSE)
i = file.copy(".wBuild/rmarkdown_show_hide_function.html",intermediates_dir)

file_input = snakemake@input[['RScript']]
file_output = snakemake@output[['wBhtml']]


sFile = tempfile()
write(spin(text=readChar(file_input, file.info(file_input)$size),knit=FALSE),sFile)
kProcessor = default_output_format(sFile)
if (kProcessor$name=="html_document")
{
	libPath = paste0(paste0(rep('../',length(strsplit(file_input,'/')[[1]])-1),collapse=''),'Output/html/libR')
	format = html_document(toc = TRUE,toc_float = TRUE,fig_retina = NULL,code_folding="show",self_contained=FALSE,lib_dir = libPath,css=c('lib/add_content_table.css','lib/leo_style.css'),df_print ='tibble')
}else{
	require(knitrBootstrap)
	format = do.call(eval(parse(text=kProcessor$name)),kProcessor$options)
}

render(file_input,output_dir = dirname(file_output),clean=TRUE,intermediates_dir = intermediates_dir, output_file = basename(file_output),output_format = format)
#render(file_input,output_dir = dirname(file_output),clean=TRUE,intermediates_dir = intermediates_dir, output_file = basename(file_output),output_format = format)
