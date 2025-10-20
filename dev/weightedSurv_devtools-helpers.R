library(usethis)
library(devtools)

# Step 1: Create the package structure
# This will create the template directory
# with Rproject setup
#usethis::create_package(".. /GitHub/weightedSurv")

# Copy functions into R directory
# Step 3: Add README and MIT license
usethis::use_readme_rmd(open = FALSE)
usethis::use_mit_license("Larry Leon")

# Step 4: Add dependencies to DESCRIPTION
desc::desc_set_dep("survival", file = "DESCRIPTION")

#usethis::use_package("ggplot2")

# Step 5: Generate documentation
# Also, run this if revising R files such as @importFrom

# clean up old documentation
unlink("man/*.Rd")

devtools::document()

devtools::load_all()

devtools::check()

devtools::clean_dll()

# In terminal
# rm -r man/
usethis::use_git()

#usethis::use_github()

# If functions are in namespace but not directly loaded
# devtools::load_all()

# Or  access hidden files:  mypackage::my_function

# Notes
# Every \item in a \describe{} block must have both a label and a description.
# Do not leave a lone \item{...} without {...} after it.
# Do not next \item inside another \item
# Use plain text, dashes, or a single paragraph for subpoints.





