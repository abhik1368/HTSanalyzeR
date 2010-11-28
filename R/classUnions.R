###############################################################################
# define class unions here
# Xin Wang <xw264@cam.ac.uk>
# Advisor: Florian Markowetz <florian.markowetz@cancer.org.uk> 
# University of Cambridge Deparment of Oncology
# Cancer Research UK - Cambridge Research Institute
# At 08:59:59, on 20 Nov 2010
###############################################################################
setClassUnion("graphNEL_Or_NULL",c("graphNEL","NULL"))
setClassUnion("numeric_Or_integer_Or_NULL",c("numeric","integer","NULL"))
setClassUnion("numeric_Or_integer",c("numeric","integer"))
setClassUnion("character_Or_NULL", c("character", "NULL"))
setClassUnion("numeric_Or_NULL", c("numeric", "NULL"))
