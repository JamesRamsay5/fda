## fda_5.1.7,  11 November 2020
This release and the last release in April of this year were prepared using maxOS Catalina version 10.15.7
This release was spell checked with spelling::spell_check_package()
Checking for reverse dependencies was via command devtools::install_github("r-lib/revdepcheck")

Two packages were found in error:

Actigraphy:  Contact was established with the maintainer and he identified an error in his code for a call to 
fda function fRegress.stderr() which he fixed.

FDRreg:  This is a tiny package.  We alerted the maintainer of the check, and then downloaded the tar file.
There were only two lines of code calling a fda functions create.bspline.basis() and eval.basis(), and their calls 
were correct.  
We will wait a week or so before resubmission in case the maintainer responds to our email.

Jim Ramsay
