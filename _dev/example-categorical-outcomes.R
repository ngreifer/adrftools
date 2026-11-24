library(adrftools)
data("nhanes3lead")

# An ordered outcome: Math score in three groups
nhanes3lead$Math3 <- cut(nhanes3lead$Math, breaks = c(-Inf, 6, 12, Inf),
                         labels = c("low", "middle", "high"),
                         ordered_result = TRUE)
table(nhanes3lead$Math3)

# Ordinal logistic regression
fit <- MASS::polr(Math3 ~ splines::ns(logBLL, df = 3) * Male +
                    Age + Race + PIR + Enough_Food,
                  data = nhanes3lead, Hess = TRUE)

# One ADRF per category (probabilities by default)
adrf1 <- adrf(fit, treat = "logBLL")
adrf1
adrf1(logBLL = c(0, 1, 2))
adrf1(logBLL = c(0, 1, 2), outcome = "high")
summary(adrf1)
plot(adrf1)
plot(adrf1, outcome = c("low", "high"))
reference_curve(adrf1, reference = 0)(logBLL = 2, outcome = "low") |> summary()

# Subgroups: contrasts are within each category
adrf2 <- adrf(fit, treat = "logBLL", by = ~Male)
adrf2(logBLL = 1, outcome = "low", subset = Male == 1)
curve_contrast(adrf2) |> summary(outcome = "low")
plot(curve_contrast(adrf2))

# Only some categories
adrf3 <- adrf(fit, treat = "logBLL", outcome = "low")
adrf3
