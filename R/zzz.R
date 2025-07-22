.onLoad <- function(libname, pkgname){

 #   library.dynam("PoweR", pkgname, libname)
    
    ns <- asNamespace(pkgname)

    reg_syms <- getDLLRegisteredRoutines(pkgname)$`.C`
    
    stat_names <- grep("^stat\\d+$", names(reg_syms), value = TRUE)
    law_names  <- grep("^law\\d+$",  names(reg_syms), value = TRUE)
    
    my_ff_call_db_stats <- lapply(stat_names, function(name) {
	function(...) {
            do.call(".C", c(list(name), list(...)))
        }
    })
    names(my_ff_call_db_stats) <- stat_names
    
    my_ff_call_db_laws <- lapply(law_names, function(name) {
        function(...) {
            do.call(".C", c(list(name), list(...)))
        }
        
    })
    names(my_ff_call_db_laws) <- law_names
    
    assign(".PoweR_stat_dispatch", my_ff_call_db_stats, envir = ns)
    assign(".PoweR_law_dispatch", my_ff_call_db_laws, envir = ns)
       
}

.onAttach <- function(lib, pkg){

    env <- as.environment("package:PoweR")
  
    source(system.file(package = "PoweR", "laws", "densities.R"), env)
    source(system.file(package = "PoweR", "laws", "moments.R"), env)
    source(system.file(package = "PoweR", "statsPureR", "statsPureR.R"), env)
    source(system.file(package = "PoweR", "computationsPureR", "pureR.R"), env)
#  for (file in list.files(system.file(package="PoweR","printTemplates"),pattern="^.*\\.R$")) source(system.file(package="PoweR","printTemplates",file),env)

}
