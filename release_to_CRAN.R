##
## 1. CHECK FOR REVISIONS TO THE PROCESS documented in 
##    the chapter on "Releasing to CRAN" in 
##    Wickham and Bryan, R Packages: 
##    https://r-pkgs.org/release.html

##
## 2. Open Project = the package file (that includes includes DESCRIPTION, etc.)
##    ... because the "process recommended "Release to CRAN" process
##    recommended by Wicham and Bryan assumes that. 

##
## 3. Pick a version number:  The version number 
##    in DESCRIPTION: must be later than any version on CRAN
##

## 
## 4. The submission process
##    Store submission comments in cran-comments.md
##

## 
## 5.  Test environments
##
# devtools::check_win_*()

devtools::check_win_devel()
devtools::check_win_release()
devtools::check_win_oldrelease()

devtools::check_rhub()

##
## 6. Reverse dependencies
##

# install.packages("revdepcheck")
#Warning in install.packages :
#  package ‘revdepcheck’ is not available for this version of R
# https://revdepcheck.r-lib.org/
# install using: 
# devtools::install_github('r-lib/revdepcheck')

revdepcheck::revdep_reset()

revdepcheck::revdep_check(num_workers = 4)

##
## 7. Update README.md and NEWS.md 
##

## 
## 8. spell check 
## 

devtools::spell_check()

##
## 9. Release to CRAN
##

release(Pkg)

## 
## 10. Update the version number for the future 
##

## 
## 11. Publicise
##
