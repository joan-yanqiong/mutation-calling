vscode_init_path <- file.path(Sys.getenv(
    if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"
), ".vscode-R", "init.R")
if (file.exists(vscode_init_path)) {
    source(vscode_init_path)
}

# Sys.setenv(RENV_PATHS_CACHE = paste0(getwd(), "/local_renv_cache:/opt/local_renv_cache:/cluster/projects/gaitigroup/Users/Joan/local_renv_cache:/Users/joankant/Desktop/gaitigroup/Users/Joan/local_renv_cache"))renv

# if (.Platform$OS.type == )
if (Sys.info()[["sysname"]] == "Darwin") {
    # If for mac store on local machine
    Sys.setenv(RENV_PATHS_CACHE = "/Users/joankant/renv_cache")
    Sys.setenv(CLANG = "/Users/joankant/miniforge-pypy3/envs/gbmcomm/bin/clang")
    Sys.setenv("CLANG++" = "/Users/joankant/miniforge-pypy3/bin/clang++")
    Sys.setenv("pkg-config" = "/Users/joankant/miniforge-pypy3/envs/gbmcomm/bin/pkg-config")
    Sys.setenv(GCC = "/Users/joankant/miniforge-pypy3/envs/gbmcomm/bin/gcc")
    Sys.setenv("libxext" = "/Users/joankant/miniforge-pypy3/envs/gbmcomm/lib/libXext.dylib")
} else if (Sys.info()[["sysname"]] == "Linux") {
    # If for linux store on cluster
    Sys.setenv(RENV_PATHS_CACHE = paste0("/Users/joankant/Desktop/gaitigroup/Users/Joan/renv_cache;/cluster/projects/gaitigroup/Users/Joan/renv_cache"), glue(getwd(), "/renv_cache"))
}

options(
    renv.config.sandbox.enabled = FALSE,
    renv.config.cache.enabled = TRUE,
    renv.settings.use.cache = TRUE,
    renv.consent = TRUE,
    renv.config.auto.snapshot = FALSE,
    renv.config.pak.enabled = TRUE
)
source("renv/activate.R")
