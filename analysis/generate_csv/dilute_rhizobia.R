#' This script takes the rhizobia input and calculate volume for making OD=0.1
#' The formula is adapted from Corlett's excel sheet entitled "Diluting rhizobia cultures to inoculate plants"


# Blank the spec using the uninoculated media

OD_undiluted <- c(0.234, 0.223, 0.213, 0.211, 0.145, 0.522)

OD_half_diluted <- c(0.119, 0.107, 0.106, 0.106, 0.08, 0.282)

# Dilution factor to get OD=0.1
dilution_factor <- 0.1 / ((OD_undiluted - OD_half_diluted) / 0.5)

# Total volume of diluted culture. Set it to the required mL
v <- 50

v * dilution_factor
v * (1-dilution_factor)

#v * dilution_factor / (dilution_factor + 1)
#v * 1 / (dilution_factor + 1)
