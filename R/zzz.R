# .onAttach <- function(...){
    # packageStartupMessage("Full description, Bug report, Suggestion and the latest codes:")
    # packageStartupMessage("https://github.com/YinLiLin/R-KAML")
# }

.onLoad <- function(libname, pkgname) {
    # Limit number of threads in veclib (MacOS MRO)
    if (Sys.info()["sysname"] == "Darwin") {
        Sys.setenv("VECLIB_MAXIMUM_THREADS" = "1")
    }
    return(invisible())
}
