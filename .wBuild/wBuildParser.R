parseWBHeader = function(filename)
{
	if(exists('snakemake') && snakemake@rule != 'WB')
	{
		filename = snakemake@input[['RScript']]
	}
	
	library(stringr)
	wbPD = 'Output/ProcessedData'
	wildcards = list(
			`wbPD` = 'Output/ProcessedData',
			`wbP` = str_split(filename,'/')[[1]][2],
			`wbPP` = str_split(filename,'/')[[1]][3],
			`wbPD_P` = paste(wbPD,str_split(filename,'/')[[1]][2],sep ="/"),
			`wbPD_PP` = paste(wbPD,str_split(filename,'/')[[1]][2],str_split(filename,'/')[[1]][3], sep='/'))
	
	if(!exists('snakemake') || snakemake@rule == 'WB')
	{
		library(methods)
		library(yaml)
		library(magrittr)
		library(gsubfn)
		
		wb = readLines(filename)%>%str_subset("^#'")%>%str_replace("^#'",'')%>%paste0(collapse='\n')%>%
				str_replace("\\n---(.|\\n)*",'')%>%yaml.load%>%`[[`('wb')
		
		repList = wildcards
		names(repList) = paste0('{',names(repList),'}')
		
		if(is.null(wb$output))
		{output = list()}
		else
		{output = as.list(gsubfn('\\{[^\\}]+\\}', repList,unlist(wb$output)))}
		
		if(is.null(wb$input))
		{input = list()}
		else
		{input = as.list(gsubfn('\\{[^\\}]+\\}', repList,unlist(wb$input)))}
		
		
		Snakemake <- setClass(
				"Snakemake",
				slots = c(
						input = "list",
						output = "list",
						params = "list",
						wildcards = "list",
						threads = "numeric",
						log = "list",
						resources = "list",
						config = "list",
						rule = "character"
				)
		)
		
		snakemake <<- Snakemake(
				input = c(input,list(filename, "RScript" = c(filename))),
				output = output,
				params = list(),
				wildcards = wildcards,
				threads = 1,
				log = list(),
				resources = list(),
				config = list(),
				rule = 'WB'
		)
		
	}
	if(snakemake@rule != 'WB')
	{
		snakemake@wildcards <<- c(wildcards,snakemake@wildcards)
	}
	
	for (w in names(wildcards))
	{
		message(w,': ', wildcards[[w]])
	}
	return(invisible(NULL))
	
}
wbReadRDS = function(name)
{
	message('Read ', snakemake@input[[name]],' ...')
	readRDS(snakemake@input[[name]])
}

wbSaveRDS = function(obj, name)
{
	message('Save ', snakemake@output[[name]],' ...')
	saveRDS(obj, snakemake@output[[name]])
}