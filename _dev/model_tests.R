# Model tests
data("nhanes3lead")
library(splines)

# lm
fit <- lm(Math ~ ns(logBLL, df = 5) * Male *
            (Age + Race + PIR + Enough_Food + Smoke_in_Home +
               Smoke_Pregnant + NICU),
          data = nhanes3lead)

a <- adrf(fit, treat = "logBLL")
plot(a, sim = F)

# glm
fit <- glm(Block >= 12 ~ ns(logBLL, df = 5) * Male *
             (Age + Race + PIR + Enough_Food + Smoke_in_Home +
                Smoke_Pregnant + NICU),
           data = nhanes3lead,
           family = "quasibinomial")

a <- adrf(fit, treat = "logBLL")
plot(a, sim = F)

# MASS::glm.nb()
fit <- MASS::glm.nb(Math ~ ns(logBLL, df = 5) * Male *
                      (Age + Race + PIR + Enough_Food + Smoke_in_Home +
                         Smoke_Pregnant + NICU),
                    data = nhanes3lead)

a <- adrf(fit, treat = "logBLL")

plot(a, sim = F)

# WeightIt::glm_weightit()
fit <- WeightIt::glm_weightit(Block >= 12 ~ ns(logBLL, df = 5) * Male *
                                (Age + Race + PIR + Enough_Food + Smoke_in_Home +
                                   Smoke_Pregnant + NICU),
                              data = nhanes3lead,
                              family = "quasibinomial")

a <- adrf(fit, treat = "logBLL")
plot(a, sim = F)

# WeightIt::glm_weightit() (with weightit)
w_fit <- WeightIt::weightit(logBLL ~ Male + Age + Race + PIR + Enough_Food +
                              Smoke_in_Home + Smoke_Pregnant + NICU,
                            data = nhanes3lead,
                            method = "glm")

fit <- WeightIt::glm_weightit(Block >= 12 ~ ns(logBLL, df = 5) * Male *
                                (Age + Race + PIR + Enough_Food + Smoke_in_Home +
                                   Smoke_Pregnant + NICU),
                              data = nhanes3lead,
                              family = "quasibinomial",
                              weightit = w_fit)

a <- adrf(fit, treat = "logBLL")
plot(a, sim = F)

# survey::svyglm()
des <- survey::svydesign(~1, weights = ~MEC_wt, data = nhanes3lead)

fit <- survey::svyglm(Block >= 12 ~ ns(logBLL, df = 5) *
                        (Age + Male + Race + PIR + Enough_Food + Smoke_in_Home +
                           Smoke_Pregnant + NICU),
                      design = des,
                      family = "quasibinomial")

a <- adrf(fit, treat = "logBLL")
plot(a, sim = F)

# dbarts::bart2() continuous
fit <- dbarts::bart2(Math ~ logBLL + Male +
                       Age + Race + PIR + Enough_Food + Smoke_in_Home +
                       Smoke_Pregnant + NICU,
                     data = nhanes3lead,
                     keepTrees = TRUE,
                     verbose = FALSE)

a <- adrf(fit, treat = "logBLL", n = 25)
plot(a, sim = F)

# dbarts::bart2() binary
fit <- dbarts::bart2(Block >= 12 ~ logBLL + Male +
                       Age + Race + PIR + Enough_Food + Smoke_in_Home +
                       Smoke_Pregnant + NICU,
                     data = nhanes3lead,
                     keepTrees = TRUE,
                     verbose = FALSE)

a <- adrf(fit, treat = "logBLL", n = 25)
plot(a, sim = F)

# estimatr::lm_robust()
fit <- estimatr::lm_robust(Math ~ ns(logBLL, df = 5) * Male *
                             (Age + Race + PIR + Enough_Food + Smoke_in_Home +
                                Smoke_Pregnant + NICU),
                           data = nhanes3lead)

a <- adrf(fit, treat = "logBLL")
plot(a, sim = F)

# fixest::feols()
fit <- fixest::feols(Math ~ ns(logBLL, df = 5) * Male *
                       (Age + Race + PIR + Enough_Food + Smoke_in_Home +
                          Smoke_Pregnant + NICU),
                     data = nhanes3lead)

a <- adrf(fit, treat = "logBLL")
plot(a, sim = F)

# fixest::feglm
fit <- fixest::feglm(Block >= 12 ~ ns(logBLL, df = 5) * Male *
                       (Age + Race + PIR + Enough_Food + Smoke_in_Home +
                          Smoke_Pregnant + NICU),
                     data = nhanes3lead,
                     family = "binomial")

a <- adrf(fit, treat = "logBLL")
plot(a, sim = F)

# fixest::fenegbin
fit <- fixest::fepois(Math ~ ns(logBLL, df = 5) * Male *
                        (Age + Race + PIR + Enough_Food + Smoke_in_Home +
                           Smoke_Pregnant + NICU),
                      data = nhanes3lead)

a <- adrf(fit, treat = "logBLL")
plot(a, sim = F)

# fixest::fenegbin
fit <- fixest::fenegbin(Math ~ ns(logBLL, df = 5) * Male *
                          (Age + Race + PIR + Enough_Food + Smoke_in_Home +
                             Smoke_Pregnant + NICU),
                        data = nhanes3lead)

a <- adrf(fit, treat = "logBLL")
plot(a, sim = F)

# rstanarm::stan_lm
fit <- rstanarm::stan_lm(Math ~ ns(logBLL, df = 5) * Male *
                           (Age + Race + PIR + Enough_Food + Smoke_in_Home +
                              Smoke_Pregnant + NICU),
                         data = nhanes3lead,
                         refresh = 0)

a <- adrf(fit, treat = "logBLL")
plot(a, sim = F)

# rstanarm::stan_glm
fit <- rstanarm::stan_glm(Block >= 12 ~ ns(logBLL, df = 5) * Male *
                            (Age + Race + PIR + Enough_Food + Smoke_in_Home +
                               Smoke_Pregnant + NICU),
                          data = nhanes3lead,
                          family = "binomial",
                          refresh = 0)

a <- adrf(fit, treat = "logBLL")
plot(a, sim = F)

# brms::brm (linear)
fit <- brms::brm(Math ~ ns(logBLL, df = 5) * Male *
                   (Age + Race + PIR + Enough_Food + Smoke_in_Home +
                      Smoke_Pregnant + NICU),
                 data = nhanes3lead,
                 refresh = 0)

a <- adrf(fit, treat = "logBLL")
plot(a, sim = F)

# brms::brm (probit)
fit <- brms::brm(Block >= 12 ~ ns(logBLL, df = 5) * Male *
                   (Age + Race + PIR + Enough_Food + Smoke_in_Home +
                      Smoke_Pregnant + NICU),
                 data = nhanes3lead,
                 family = binomial("probit"),
                 refresh = 0)

a <- adrf(fit, treat = "logBLL")
plot(a, sim = F)
