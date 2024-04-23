#' This script makes the table of statistics

renv::load()
library(tidyverse)
library(lme4)
library(sjPlot)
library(webshot)
source(here::here("metadata.R"))

data("sleepstudy")
data("efc")
efc$cluster <- as.factor(efc$e15relat)
m1 <- lmer(neg_c_7 ~ c160age + c161sex + e42dep + (1 | cluster), data = efc)
m2 <- lmer(Reaction ~ Days + (1 + Days | Subject), data = sleepstudy)

tab_model(m1, m2, file = here::here("plots/Tab1.html"))
webshot(here::here("plots/Tab1.html"), here::here("plots/Tab1.png"), vwidth = 15, vheight = 10, zoom = 5)
