vscode_init_path <- file.path(Sys.getenv(
    if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"
), ".vscode-R", "init.R")
if (file.exists(vscode_init_path)) {
    source(vscode_init_path)
}

if (Sys.info()[["sysname"]] == "Darwin") {
    # If for mac store on local machine
    Sys.setenv(RENV_PATHS_CACHE = paste0("/Users/joankant/OneDrive - UHN/renv_cache"))
    Sys.setenv(CLANG = "/Users/joankant/miniforge-pypy3/envs/standard_env/bin/clang")
    Sys.setenv("CLANG++" = "/Users/joankant/miniforge-pypy3/bin/clang++")
    Sys.setenv("pkg-config" = "/Users/joankant/miniforge-pypy3/envs/standard_env/bin/pkg-config")
    Sys.setenv(GCC = "/Users/joankant/miniforge-pypy3/envs/standard_env/bin/gcc")
    Sys.setenv("libxext" = "/Users/joankant/miniforge-pypy3/envs/standard_env/lib/libXext.dylib")
} else if (Sys.info()[["sysname"]] == "Linux") {
    # If for linux store on cluster
    Sys.setenv(RENV_PATHS_CACHE = paste0("cluster/home/t119972uhn/renv_cache"))
}

options(
    renv.config.sandbox.enabled = FALSE,
    renv.config.cache.enabled = TRUE,
    renv.settings.use.cache = TRUE,
    renv.consent = TRUE,
    renv.config.auto.snapshot = FALSE,
    renv.config.pak.enabled = TRUE,
    renv.config.snapshot.type = "implicit"
)
source("renv/activate.R")
