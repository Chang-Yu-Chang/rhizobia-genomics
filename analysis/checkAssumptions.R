# AUTHOR:   Corlett Wolfe Wood
# DATE:     7 April 2016

# PURPOSE:  Create plots to check three major regression assumptions:

# 1) Normality: Histogram of the residuals and a qqplot
#    Are the residuals normally distributed?
#    -- Plot a histogram of the residuals
#    -- Create a normal quantile-quantile plot

# 2) Linearity: Fitted vs residuals plot
#    Is there a linear relationship between the dependent and independent variables?
#    -- Plot the residuals versus the fitted values with a loess line
#    -- Slope should be 0

# 3) Homoscedasticity: Scale-location plot
#    Is the absolute value of the residuals constant across the range of x?
#    -- Plot the sq. root of the abs. value of residuals against fitted values
#    -- A slope != 0 shows a change in spread

#------------------------------------------------------------------------------#
# SOURCES
#   https://data.library.virginia.edu/diagnostic-plots/
#   http://glmm.wdfiles.com/local--files/examples/Banta_2011_part1.pdf
#------------------------------------------------------------------------------#

# FUNCTION: check.assumptions
check.assumptions <- function(mod){
  # Multiple panels
  par(mfrow=c(2,2))

  # Normality
  hist(resid(mod))
  qqnorm(resid(mod))
  qqline(resid(mod))

  # Linearity
  # --> Plot residual versus fitted values with loess line
  plot(fitted(mod), resid(mod), main = "Residual vs fitted values")
  fr.fit <- loess(resid(mod)~fitted(mod))
  fr.vec <- seq(min(fitted(mod)), max(fitted(mod)), length.out=100)
  lines(fr.vec,predict(fr.fit,fr.vec), col="red")

  # Homoscedasticity
  # --> Plot the absolute value of the residuals against the fitted values
  plot(fitted(mod), sqrt(abs(resid(mod))), main="Scale-location")
  sl.fit <- loess(sqrt(abs(resid(mod)))~fitted(mod))
  sl.vec <- seq(min(fitted(mod)), max(fitted(mod)), length.out=100)
  lines(sl.vec,predict(sl.fit,sl.vec), col="red")
}


# Function to calculate variance inflation factors
#SOURCE: https://jonlefcheck.net/2012/12/28/dealing-with-multicollinearity-using-variance-inflation-factors/
vif.lme <- function (fit) {
    ## adapted from rms::vif
    v <- vcov(fit)
    nam <- names(fixef(fit))
    ## exclude intercepts
    ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
    if (ns > 0) {
        v <- v[-(1:ns), -(1:ns), drop = FALSE]
        nam <- nam[-(1:ns)] }
    d <- diag(v)^0.5
    v <- diag(solve(v/(d %o% d)))
    names(v) <- nam
    v}
