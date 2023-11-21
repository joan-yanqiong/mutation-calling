vscode_init_path <- file.path(Sys.getenv(
    if (.Platform$OS.type == "windows") "USERPROFILE" else "HOME"
), ".vscode-R", "init.R")
if (file.exists(vscode_init_path)) {
    source(vscode_init_path)
}

if (Sys.info()[["sysname"]] == "Darwin") {
    # If for mac store on local machine
    Sys.setenv(RENV_PATHS_CACHE = "~/Desktop/gaitigroup/Users/Joan/renv_cache")
    Sys.setenv(CLANG = "~/miniforge3/envs/standard_env/bin/clang")
    Sys.setenv("CLANG++" = "~/miniforge3/bin/clang++")
    Sys.setenv("pkg-config" = "~/miniforge3/envs/standard_env/bin/pkg-config")
    Sys.setenv(GCC = "~/miniforge3/envs/standard_env/bin/gcc")
    Sys.setenv("libxext" = "~/miniforge3/envs/standard_env/lib/libXext.dylib")
} else if (Sys.info()[["sysname"]] == "Linux") {
    # If for linux store on cluster
    Sys.setenv(RENV_PATHS_CACHE = paste0("~/renv_cache"))
}

.libPaths(c(.libPaths(), "~/Desktop/gaitigroup/Users/Joan/h4h-mutation-calling/renv/library/R-4.2/x86_64-apple-darwin13.4.0"))

options(
    renv.config.sandbox.enabled = FALSE,
    renv.config.cache.enabled = TRUE,
    renv.settings.use.cache = TRUE,
    renv.consent = TRUE,
    renv.config.auto.snapshot = FALSE,
    renv.config.pak.enabled = FALSE,
    renv.config.snapshot.type = "explicit"
)
# source("renv/activate.R")
