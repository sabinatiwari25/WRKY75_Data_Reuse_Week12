# Restore the usual Windows user library when a sandbox omits R_LIBS_USER.
user_library <- file.path(Sys.getenv('LOCALAPPDATA'), 'R', 'win-library',
                         paste(R.version$major, strsplit(R.version$minor, '.', fixed=TRUE)[[1]][1], sep='.'))
fallback_library <- file.path(Sys.getenv('USERPROFILE'),'AppData','Local','R','win-library',basename(user_library))
if (.Platform$OS.type=='windows') .libPaths(c(user_library,fallback_library,.libPaths()))
if (dir.exists('.Rlib')) .libPaths(c(normalizePath('.Rlib'), .libPaths()))
