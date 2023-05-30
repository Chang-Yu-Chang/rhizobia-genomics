#' This script fit the to

require(tidyverse)
library(janitor)
require(viridis)
require(data.table)
require(mgcv)
require(egg)
source(here::here("analysis/00-metadata.R"))

compute.gam <- function(x)
{
    # Remove "blank"
    x <- x[order(t)]
    x[, abs := abs-0.034]
    blk <- x$abs[1]

    if (any(x$abs>.05, na.rm=TRUE)) {

        x[, abs := ifelse(abs<0.006, 0.006, abs)]
        x[, lOD := log(abs)]

        # First smooth the raw data, which is needed for the lag phase
        gam0 <- gam(abs ~ s(t, bs = 'ad'), data=x)

        t <- unique(x$t)
        newd <- data.frame(t = t)
        pred0 <- predict(gam0, newd,
                         se.fit = TRUE) %>%
            as_tibble() %>%
            mutate(t = t)

        X0 <- predict(gam0, newd, type="lpmatrix")

        eps <- 1e-7 ## finite difference interval
        newd <- data.frame(t = t+eps)
        X1 <- predict(gam0, newd, type="lpmatrix")

        Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives

        df <- Xp%*%coef(gam0)              ## ith smooth derivative
        df.sd <- rowSums(Xp%*%gam0$Vp*Xp)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5

        pred0 <- pred0 %>%
            mutate(deriv.fit = df[,1], deriv.sd = df.sd,
                   csource = x$csource[1],
                   isolate = x$seq[1])

        # Now model the log OD for getting the growth rate
        gam1 <- try(gam(lOD ~ s(t, bs = 'ad'), data=x), silent=TRUE)

        t <- unique(x$t)
        newd <- data.frame(t = t)
        pred <- predict(gam1, newd,
                        se.fit = TRUE) %>%
            as_tibble() %>%
            mutate(t = t)

        X0 <- predict(gam1, newd, type="lpmatrix")

        eps <- 1e-5 ## finite difference interval
        newd <- data.frame(t = t+eps)
        X1 <- predict(gam1, newd, type="lpmatrix")

        Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives

        df <- Xp%*%coef(gam1)              ## ith smooth derivative
        df.sd <- rowSums(Xp%*%gam1$Vp*Xp)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5

        # GET THE DATA
        pred <- pred %>%
            mutate(deriv.fit = df[,1], deriv.sd = df.sd,
                   csource = x$csource[1],
                   seq = x$seq[1],
                   OD.frac = pred0$fit/max(pred0$fit),
                   lOD = x$lOD,
                   well=x$well
                   # date=as.character(x$date)
            )

        # GET GROWTH PARAMETERS

        # discard all the first X points with OD lower than 0.015 or t<1h
        aux <- rle(pred$lOD<log(0.015))
        aux <- ifelse(aux$values[1], pred$t[aux$lengths[1]], 0)
        aux <- ifelse(aux<1, 1, aux)

        ## compute lag
        maxgr <- pred %>%
            filter(t>aux) %>%
            filter(deriv.fit==max(deriv.fit))

        t.maxGr <- maxgr$t[1]

        OD.maxGr <- exp(maxgr$fit[1])

        p0 <- setDT(pred0)[t==t.maxGr]
        slope.maxGr <- p0$deriv.fit

        ## y = y_0 + slope * t, so:
        y_0 <- OD.maxGr - (slope.maxGr*t.maxGr)

        # I take as lag the intersection with y = min (OD)
        t_lag <- (min(pred0$fit)-y_0)/slope.maxGr # lag

        prm <- tibble(seq = x$seq[1],
                      well = x$well[1],
                      csource = x$csource[1],
                      strain = x$strain[1],
                      # date = as.character(x$date[1]),
                      r = maxgr$deriv.fit[1],
                      t.r = t.maxGr,
                      lag = t_lag,
                      maxOD.fit = max(pred0$fit),
                      startOD = blk )
        return(list('data' = pred,
                    'params' = prm))
    } else {
        prm <- tibble(seq = x$seq[1],
                      well = x$well[1],
                      csource = x$csource[1],
                      strain = x$strain[1],
                      # date = as.character(x$date[1]),
                      r = 0,
                      t.r = NA,
                      lag = NA,
                      maxOD.fit = NA,
                      startOD = blk)
        return(list('data' = NA,
                    'params' = prm))
    }
}


gc <- read_csv(paste0(folder_data, "raw/growth_curve/rhizobia_growth_curve.csv"))
gc_plate <- read_csv(paste0(folder_data, "raw/growth_curve/gc_plate.csv"))
clean_well_names <- function (x) {
    y <- paste0(str_sub(x, 1, 1), str_sub(x, 2, 3) %>% as.numeric %>% sprintf("%02d", .))
    return(y)
}

# Tidy up
gc <- gc %>%
    # Calculate the time interval in minutes
    mutate(t = as.numeric(difftime(Time, Time[1], units = "hours"))) %>%
    #mutate(Time = hour(gc$Time) + minute(gc$Time)/60) %>%
    select(-Time, -`T° 600`) %>%
    pivot_longer(cols = -t, names_to = "well", values_to = "od600") %>%
    mutate(well = clean_well_names(well)) %>%
    left_join(gc_plate) %>%
    filter(strain != "blank") %>%
    rename(abs = od600)

gc <- as.data.table(gc)
gc <- split(gc, by = c('well', 'strain'))
gc.fits <- lapply(gc, compute.gam)
gc.prm <- do.call(rbind, lapply(gc.fits, function(x) x[[2]]))
gc <- do.call(rbind, lapply(gc.fits, function(x) x[[1]]))

## get statistics over replicates
gc.prm.stat <- gc.prm %>%
    group_by(strain) %>%
    summarize(r.sem = sd(r)/sqrt(n()), r = mean(r),
              t.r.sem = sd(t.r)/sqrt(n()), t.r = mean(t.r),
              lag.sem = sd(lag)/sqrt(n()), lag = mean(lag),
              maxOD.sem = sd(maxOD.fit)/sqrt(n()), maxOD = mean(maxOD.fit))

setDT(gc.prm)
setDT(gc.prm.stat)
fwrite(gc.prm, file=paste0(folder_data, 'temp/04c-gc_prm.csv'))
fwrite(gc.prm.stat, file=paste0(folder_data, 'temp/04c-gc_prm_summ.csv'))

#d <- fread(paste0(folder_data, 'raw/growth_curve/test_fit/raw_gcurves_all.csv'))
#d <- split(d, by = c(well', 'strain', 'seq'))
#d.fits <- lapply(d, compute.gam)
# d.prm <- do.call(rbind, lapply(d.fits, function(x) x[[2]]))
# d <- do.call(rbind, lapply(d.fits, function(x) x[[1]]))

## get statistics over replicates
# d.prm.stat <- d.prm %>%
#     group_by(seq, csource) %>%
#     summarize(r.sem = sd(r)/sqrt(n()), r = mean(r),
#               t.r.sem = sd(t.r)/sqrt(n()), t.r = mean(t.r),
#               lag.sem = sd(lag)/sqrt(n()), lag = mean(lag),
#               maxOD.sem = sd(maxOD.fit)/sqrt(n()), maxOD = mean(maxOD.fit))
#
# setDT(d.prm.stat)
