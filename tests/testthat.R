library("base")
requireNamespace("testthat")
library("w4mclstrpeakpics")

# The following line causes a warning in devtools::check() -
#   WARNING '::' or ':::' import not declared from: ‘testthat’
# but I have found no way to suppress it without importing testthat into the package.
# The problem is that 'testthat' is not needed by end users of the package.
testthat::test_check("w4mclstrpeakpics")
