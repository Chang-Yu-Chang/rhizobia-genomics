#' This script fits the growth curve data

library(tidyverse)
library(data.table)
library(mgcv) # For smoothing
source(here::here("metadata.R"))

gc_plate <- read_csv(paste0(folder_data, "raw/growth_curves/growth_curve2/gc_plate.csv")) # plate layout
list_gcs <- rep(list(NA), 4)
list_gcs[[1]] <- read_csv(paste0(folder_data, "raw/growth_curves/growth_curve2/rhizobia_growth_curve.csv")) # 30C
list_gcs[[2]] <- read_csv(paste0(folder_data, "raw/growth_curves/growth_curve3/rhizobia_growth_curve.csv")) # 35C
list_gcs[[3]] <- read_csv(paste0(folder_data, "raw/growth_curves/growth_curve4/rhizobia_growth_curve.csv")) # 25C
list_gcs[[4]] <- read_csv(paste0(folder_data, "raw/growth_curves/growth_curve5/rhizobia_growth_curve.csv")) # 40C
names(list_gcs) <- c("30c", "35c", "25c", "40c")

wide_to_long <- function (gc) {
    #' This function clean up the variable names and pivot the wide form to long
    clean_well_names <- function (x) paste0(str_sub(x, 1, 1), str_sub(x, 2, 3) %>% as.numeric %>% sprintf("%02d", .))

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
# gc <- list_gcs[[2]]
# x <- gcl %>% filter(well == "G03")
compute_gam <- function(x) {
    #' Smooth and fit GAM to a growth curve
    x <- x %>%
        mutate(abs = ifelse(abs<0.01, 0.01, abs)) %>% # Avoid 0 value
        mutate(lod = log(abs))

    if (any(x$abs > 0.01, na.rm = T)) {

        # First smooth the raw data, which is needed for the lag phase
        gam0 <- gam(abs ~ s(t, bs = "ad"), data = x)

        ut <- unique(x$t)
        newd <- tibble(t = ut)
        X0 <- predict(gam0, newd, type="lpmatrix")
        eps <- 1e-7 ## finite difference interval
        newd <- tibble(t = ut + eps)
        X1 <- predict(gam0, newd, type="lpmatrix")
        Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives

        df <- Xp%*%coef(gam0)              ## ith smooth derivative
        df.sd <- rowSums(Xp%*%gam0$Vp*Xp)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5

        pred0 <- predict(gam0, newd, se.fit = T) %>%
            as_tibble() %>% mutate(t = ut) %>%
            mutate(deriv.fit = df[,1], deriv.sd = df.sd)

        # Now model the log OD for getting the growth rate
        gam1 <- gam(lod ~ s(t, bs = 'ad'), data = x)

        newd <- tibble(t = ut)
        pred <- predict(gam1, newd, se.fit = TRUE) %>% as_tibble() %>% mutate(t = ut)
        X0 <- predict(gam1, newd, type="lpmatrix")
        eps <- 1e-5 ## finite difference interval
        newd <- tibble(t = ut+eps)
        X1 <- predict(gam1, newd, type="lpmatrix")
        Xp <- (X1-X0)/eps ## maps coefficients to (fd approx.) derivatives

        df <- Xp%*%coef(gam1)              ## ith smooth derivative
        df.sd <- rowSums(Xp%*%gam1$Vp*Xp)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5

        # GET THE DATA
        pred <- pred %>%
            mutate(
                deriv.fit = df[,1],
                deriv.sd = df.sd,
                OD.frac = pred0$fit/max(pred0$fit),
                lod = x$lod
            )

        # GET GROWTH PARAMETERS
        # discard all the first X points with OD lower than 0.001 or t<1h
        aux <- rle(pred$lod<log(0.001))
        aux <- ifelse(aux$values[1], pred$t[aux$lengths[1]], 0)
        aux <- ifelse(aux<1, 1, aux)

        ## compute lag
        # Max growth rate
        rmax <- pred %>%
            filter(t > aux) %>%
            filter(deriv.fit == max(deriv.fit))
        t_rmax <- rmax$t[1] # time at which max gorwth rate
        od_rmax <- exp(rmax$fit[1]) # od at which max growth rate
        p0 <- pred0 %>% filter(t == t_rmax)
        slope <- p0$deriv.fit # slope at the absolute od scale

        # At time t1 when r is max: y1 = y0 + slope * t1, where y1 is the od at time t1, the slope is at absolute od not log od
        # The intercept is then
        y0 <- od_rmax - (slope*rmax$t[1])

        # I take lag as the intersection with y = min (OD)
        # Find the intersection with y = slope * t + y0
        t_lag <- (min(pred0$fit)-y0)/slope # lag

        prm <- tibble(
            r = rmax$deriv.fit,
            t_rmax = t_rmax, # time at which read the maximum growth rate
            lag = t_lag,
            maxOD = max(pred0$fit)
        )
    } else {
        pred <- NA
        prm <- tibble(
            r = NA,
            t_rmax = NA,
            lag = NA,
            maxOD = NA
        )
    }
    return(list(data = pred, params= prm))
}
summarize_gc <- function (gc) {
    #' This function takes the  (after blank) and average across replicates
    gc %>%
        group_by(t, exp_id) %>%
        summarize(mean_abs = mean(abs), sd_abs = sd(abs)) %>%
        ungroup()
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
    # Format
    gcl <- wide_to_long(gc) # to long format
    gcl <- gcl %>% filter(t <= 48)

    # Blank
    gcl_blank <- extract_blank(gcl)
    gcl <- subtract_blank(gcl, gcl_blank)

    # Smooth growth curves
    gc_smooth <- gcl %>%
        nest(data = c(-well, -exp_id)) %>%
        mutate(prm = map(data, compute_gam))

    # Smoothed curve
    gcl_smooth <- gc_smooth %>%
        unnest(prm) %>%
        slice(seq(1, n(), 2)) %>%
        select(well, exp_id, prm) %>%
        unnest(prm) %>%
        mutate(abs_fit = exp(fit)) %>%
        select(t, well, exp_id, abs_fit)

    gcl_smooth$abs_fit <- as.vector(gcl_smooth$abs_fit)

    # Parameter
    gtw <- gc_smooth %>%
        unnest(prm) %>%
        slice(seq(2, n(), 2)) %>%
        select(well, exp_id, prm) %>%
        unnest(prm)


    gts <- gtw %>%
        group_by(exp_id) %>%
        summarize(r.sem = sd(r)/sqrt(n()), r = mean(r),
                  lag.sem = sd(lag)/sqrt(n()), lag = mean(lag),
                  maxOD.sem = sd(maxOD)/sqrt(n()), maxOD = mean(maxOD))

    return(list(gcl = gcl, gcl_smooth = gcl_smooth, gcl_blank = gcl_blank, gtw = gtw, gts = gts))
})

gcl <- list_gcs %>% lapply(function(x) `[[`(x, "gcl")) %>% bind_rows(.id = "temperature")
gcl_smooth <- list_gcs %>% lapply(function(x) `[[`(x, "gcl_smooth")) %>% bind_rows(.id = "temperature")
gtw <- list_gcs %>% lapply(function(x) `[[`(x, "gtw")) %>% bind_rows(.id = "temperature")
gts <- list_gcs %>% lapply(function(x) `[[`(x, "gts")) %>% bind_rows(.id = "temperature")

write_csv(gcl, paste0(folder_phenotypes, 'growth/gcl.csv')) # Raw growth curves
write_csv(gcl_smooth, paste0(folder_phenotypes, 'growth/gcl_smooth.csv')) # Smooth growth curves
write_csv(gtw, paste0(folder_phenotypes, 'growth/gtw.csv')) # Growth traits per well
write_csv(gts, paste0(folder_phenotypes, 'growth/gts.csv')) # Growth traits per isolate
