info <- read.delim("C:/Users/2026mr001/Desktop/Stage/Identification cif/cif-like/Check conta/Tentative/info.txt")
info <- info[info$Clade != "Insert", ]
library(betareg)
library(emmeans)

# GC ~ Clade
par(mfrow = c(1,1))
if (max(info$GC, na.rm = TRUE) > 1) {
  info$GC <- info$GC / 100
}
range(info$GC)

mod <- betareg(GC~ Clade, data = info)
summary(mod)
emm <- emmeans(mod, ~ Clade, type = "response")
plot(emm)

# coding_density  ~ Clade

if (max(info$coding_density, na.rm = TRUE) > 1) {
  info$coding_density <- info$coding_density / 100
}
mod <- betareg(coding_density ~ Clade, data = info)
summary(mod)
emm <- emmeans(mod, ~ Clade, type = "response")
plot(emm)

res <- aov(coding_density ~ Clade, data = info)

# hypotheticals vs annotated ~ Clade

library(emmeans)
mod <- glm(cbind(CDSs - hypotheticals, hypotheticals) ~ Clade,
           data = info,
           family = quasibinomial)

emm <- emmeans(mod, ~ Clade, type = "response")
plot(emm)

sum(residuals(mod, type = "pearson")^2) / df.residual(mod)
Anova(mod)
summary(mod)


