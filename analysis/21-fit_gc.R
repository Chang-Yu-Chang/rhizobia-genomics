#' This script fits the growth curve data

renv::load()
suppressPackageStartupMessages({
    library(tidyverse)
    library(janitor)
    library(data.table)
    library(mgcv)
    source(here::here("analysis/00-metadata.R"))
})

gc_plate <- read_csv(paste0(folder_data, "raw/growth_curve2/gc_plate.csv"), show_col_types = F)

list_gcs <- rep(list(NA), 3)
list_gcs[[1]] <- read_csv(paste0(folder_data, "raw/growth_curve2/rhizobia_growth_curve.csv"), show_col_types = F) # 30C
list_gcs[[2]] <- read_csv(paste0(folder_data, "raw/growth_curve3/rhizobia_growth_curve.csv"), show_col_types = F) # 35C
list_gcs[[3]] <- read_csv(paste0(folder_data, "raw/growth_curve4/rhizobia_growth_curve.csv"), show_col_types = F) # 25C
names(list_gcs) <- c("30c", "35c", "25c")

# Tidy up time series
clean_well_names <- function (x) paste0(str_sub(x, 1, 1), str_sub(x, 2, 3) %>% as.numeric %>% sprintf("%02d", .))
wide_to_long <- function (gc) {
    #' This function clean up the variable names and pivot the wide form to long
    gc %>%
        # Calculate the time interval in minutes
        mutate(t = as.numeric(difftime(Time, Time[1], units = "hours"))) %>%
        select(-Time, -`T° 600`) %>%
        pivot_longer(cols = -t, names_to = "well", values_to = "od600") %>%
        mutate(well = clean_well_names(well)) %>%
        left_join(gc_plate) %>%
        rename(abs = od600)

}
extract_blank <- function (gc) {
    #' This function takes the long format of gc to calculate the average abs of blank at each time point
    gc %>%
    filter(exp_id == "blank") %>%
        group_by(t) %>%
        summarize(abs_blank = mean(abs))
}
subtract_blank <- function (gc, gc_blank) {
    #' This function takes the blank and subtract it
    gc %>%
        left_join(gc_blank) %>%
        mutate(abs = abs - abs_blank) %>%
        mutate(abs = ifelse(abs < 0, 0, abs)) %>%
        filter(exp_id != "blank")

}
summarize_gc <- function (gc) {
    #' This function takes the gc (after blank) and average across replicates
    gc %>%
    group_by(t, exp_id) %>%
        summarize(mean_abs = mean(abs), sd_abs = sd(abs))
}
# Smooth and fit GAM to growth curves
compute.gam <- function(x) {
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
                      exp_id = x$exp_id[1],
                      # date = as.character(x$date[1]),
                      r = maxgr$deriv.fit[1],
                      t.r = t.maxGr,
                      lag = t_lag,
                      maxOD = max(pred0$fit),
                      startOD = blk )
        return(list('data' = pred,
                    'params' = prm))
    } else {
        prm <- tibble(seq = x$seq[1],
                      well = x$well[1],
                      csource = x$csource[1],
                      exp_id = x$exp_id[1],
                      # date = as.character(x$date[1]),
                      r = 0,
                      t.r = NA,
                      lag = NA,
                      maxOD = NA,
                      startOD = blk)
        return(list('data' = NA,
                    'params' = prm))
    }
}
calculate_prm <- function (gc) {
    #' This function takes the gc time series and calculate the gc parameters
    gc <- as.data.table(gc)
    gc <- split(gc, by = c('well', 'exp_id'))
    gc_fits <- lapply(gc, compute.gam)
    gc_prm <- do.call(rbind, lapply(gc_fits, function(x) x[[2]]))
    return(gc_prm)
}
summarize_prm <- function (gc_prm) {
    #' This function gets the statistics over replicates
    gc_prm %>%
        group_by(exp_id) %>%
        summarize(r.sem = sd(r)/sqrt(n()), r = mean(r),
                  t.r.sem = sd(t.r)/sqrt(n()), t.r = mean(t.r),
                  lag.sem = sd(lag)/sqrt(n()), lag = mean(lag),
                  maxOD.sem = sd(maxOD)/sqrt(n()), maxOD = mean(maxOD))

}


list_gcs <- list_gcs %>% lapply(function(gc) {
        gc <- wide_to_long(gc)
        gc_blank <- extract_blank(gc)
        gc <- subtract_blank(gc, gc_blank)
        gc_summ <- summarize_gc(gc)
        gc_prm <- calculate_prm(gc)
        gc_prm_summ <- summarize_prm(gc_prm)
        #
        gc <- gc %>% filter(t <= 48)
        gc_summ <- gc_summ %>% filter(t <= 48)
        return(list(gc = gc, gc_summ = gc_summ, gc_prm = gc_prm, gc_prm_summ = gc_prm_summ))
    })
gcs <- list_gcs %>% lapply(function(x) `[[`(x, "gc")) %>% bind_rows(.id = "temperature")
gc_summs <- list_gcs %>% lapply(function(x) `[[`(x, "gc_summ")) %>% bind_rows(.id = "temperature")
gc_prms <- list_gcs %>% lapply(function(x) `[[`(x, "gc_prm")) %>% bind_rows(.id = "temperature")
gc_prm_summs <- list_gcs %>% lapply(function(x) `[[`(x, "gc_prm_summ")) %>% bind_rows(.id = "temperature")


write_csv(gcs, paste0(folder_data, 'temp/21-gcs.csv'))
write_csv(gc_summs, paste0(folder_data, 'temp/21-gc_summs.csv'))
write_csv(gc_prms, paste0(folder_data, 'temp/21-gc_prms.csv'))
write_csv(gc_prm_summs, paste0(folder_data, 'temp/21-gc_prm_summs.csv'))
