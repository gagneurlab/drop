
#' rmarkdown_show_hide_table
#' 
#' function to hide/show a container within a rmarkdown html
#' needs to be used with the code chunk option ==> results='asis'
#' 
#' @example within R
#' 		divname='raw'; button='show/hide new table'
#' 		#' <a id="displayText" href="javascript:toggle('`r divname`');">`r button`</a>
#' 		#' <div id="`r divname`" style="display: none" href="javascript:toggle();">
#' 		data_object
#' 		#' </div>
#' 		#' <a id="displayText" href="javascript:toggle('`r divname`');">`r button`</a>
#'
rmarkdown_show_hide_table <- function(name, data_view_cmd, view_function='eval'){

	# generate arbitraty table name; avoid duplicates by chance
	options(digits.secs=6)
	rand_table_name <- paste0("table_", sub(' ','_', Sys.time()))
	rand_id <- paste0("displayText_", sub(' ','_', Sys.time()))
	
	# def html table button
	table_button <- paste0("<p><a id=\"",rand_id,"\" href=\"javascript:toggle('", rand_table_name, 
			"');\">show/hide ", name, "</a></p>")
	
	# put button above container
	cat(table_button)
	
	# print container
	cat("\n<div id=\"", rand_table_name, "\" style=\"display:none\" href=\"javascript:toggle();\">\n", sep="")
	if(is.null(view_function)){
		writeLines('\n')
#		writeLines(data_view_cmd)
		data_view_cmd
		writeLines('\n')
	}else{
		cat("\n<pre>\n<code>\n\n")
		get(view_function)(data_view_cmd)
		cat("\n\n</code>\n</pre>\n\n")
	}
	cat("\n</div>\n")
	
	# put button below container
	cat(table_button)
}

#' knitr_show_hide_container
#' 
#' Shortcut for rmarkdown_show_hide_table, e.g. if data_view_cmd returns directly HTML code
#' 
knitr_show_hide_container <- function(name, data_view_cmd){
	rmarkdown_show_hide_table(name, data_view_cmd, view_function =NULL)
}

rmarkdown_show_hide_table2 <- function(name, data_view_cmd, view_function='eval'){
	
	# generate arbitraty table name; avoid duplicates by chance
	options(digits.secs=6)
	rand_table_name <- paste0("table_", sub(' ','_', Sys.time()))
	
	# put button above
	table_button <- paste("<a id=\"displayText\" href=\"javascript:toggle('", rand_table_name, 
			"');\">show/hide ", name, "</a>", sep="")
	
	# create output
	cat(table_button)
	
	# print container
	cat("\n<div id=\"", rand_table_name, "\" style=\"display: none\" href=\"javascript:toggle();\">\n", sep="")
	if(is.null(view_function)){
		writeLines('')
		writeLines(data_view_cmd)
	}else{
		#cat("\n<pre>\n<code>\n\n")
		get(view_function)(data_view_cmd)
		#cat("\n\n</code>\n</pre>\n\n")
	}
	cat("\n")
	cat("</div>\n")
	
	# put button below container
	cat(paste0("<p>",table_button, "</p>"))
}


