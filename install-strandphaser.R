#!/usr/bin/env Rscript
Sys.setenv(TAR = '/bin/tar')
library(devtools)
withr::with_libpaths(new = "R-packages", install_git("git://github.com/daewoooo/StrandPhaseR.git", branch = "master"))

