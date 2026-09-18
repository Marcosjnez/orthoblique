library(latent)
library(plot.matrix)

# https://journals.sagepub.com/doi/abs/10.1177/21677026211055170
# https://osf.io/hvskz/

#### Extract the data ####

R <- mtmm$ABCD$R # Extract the correlation matrix
R[, c("ARGUc", "ARGUp")] # duplicated items
remove <- -match(c("ARGUc", "ARGUp"), rownames(R))
R <- R[remove, ][, remove]
items <- rownames(R)

#### Separate traits and methods ####

# Traits:
Attention_Problems <- c("IMATc", "FFINc", "CONCc", "RESTc", "IMPLc", "DISTc",
                        "IMATp", "FFINp", "CONCp", "RESTp", "IMPLp", "DISTp",
                        "IMATt", "FFINt", "CONCt", "RESTt", "IMPLt", "DISTt")
Externalizing <- c("DESTc", "DOHOc", "DISCc", "IRRIc", "TEMPc", "THRTc",
                   "DESTp", "DOHOp", "DISCp", "IRRIp", "TEMPp", "THRTp",
                   "ARGUt", "DESTt", "DISCt", "IRRIt", "TEMPt", "THRTt")
Internalizing <- c("FEARc", "GUILc", "EMBAc", "WRTHc", "DEPRc", "WORRc",
                   "FEARp", "GUILp", "EMBAp", "WRTHp", "DEPRp", "WORRp",
                   "FEARt", "GUILt", "EMBAt", "WRTHt", "DEPRt", "WORRt")

# Methods:
Children <- c("IMATc", "FFINc", "CONCc", "RESTc", "IMPLc", "DISTc",
              "DESTc", "DOHOc", "DISCc", "IRRIc", "TEMPc", "THRTc",
              "FEARc", "GUILc", "EMBAc", "WRTHc", "DEPRc", "WORRc")
Parents <- c("IMATp", "FFINp", "CONCp", "RESTp", "IMPLp", "DISTp",
             "DESTp", "DOHOp", "DISCp", "IRRIp", "TEMPp", "THRTp",
             "FEARp", "GUILp", "EMBAp", "WRTHp", "DEPRp", "WORRp")
Teachers <- c("IMATt", "FFINt", "CONCt", "RESTt", "IMPLt", "DISTt",
              "ARGUt", "DESTt", "DISCt", "IRRIt", "TEMPt", "THRTt",
              "FEARt", "GUILt", "EMBAt", "WRTHt", "DEPRt", "WORRt")

#### Create target matrices ####

Target <- matrix(0, nrow = length(items), ncol = 6)
rownames(Target) <- items
Target[Attention_Problems, 1] <- 1
Target[Externalizing, 2] <- 1
Target[Internalizing, 3] <- 1
Target[Children, 4] <- 1
Target[Parents, 5] <- 1
Target[Teachers, 6] <- 1
PsiTarget <- diag(6); PsiTarget[1:3, 1:3] <- 1; PsiTarget[3:6, 3:6] <- 1
diag(PsiTarget) <- 0

plot(Target, breaks = c(0, 0.5, 1),
     main = "Target on factor loadings",
     col = c("black", "yellow"))
plot(PsiTarget, breaks = c(0, 0.5, 1),
     main = "Target on factor correlations",
     col = c("black", "yellow"))

#### MTMM analysis POB (Target; 3 Traits + 3 Methods) ####

set.seed(2026)
efa <- lefa(sample.cov = R, sample.nobs = 2119,
            nfactors = 6, estimator = "uls",
            projection = "poblq", oblique = c(3, 3),
            rotation = "target", target = Target,
            se = FALSE)

#### Plots ####

# Align the results to the theoretical Target:
lambda <- efa@transformed_pars$lambda_rotated
psi <- efa@transformed_pars$psi_rotated; diag(psi) <- 1
x <- fungible::faAlign(F1 = Target, F2 = abs(lambda), Phi2 = abs(psi))
lambda <- lambda[, c(x$FactorMap[2, ])]
psi <- psi[c(x$FactorMap[2, ]), c(x$FactorMap[2, ])]

# Dimnames:
rownames(lambda) <- NULL
colnames(lambda) <- c("Attention Problems", "Externalizing", "Internalizing",
                      "Children", "Parents", "Teachers")
rownames(psi) <- colnames(psi) <- c("Attention Problems", "Externalizing", "Internalizing",
                                    "Children", "Parents", "Teachers")

model_name <- "target"

# Plot absolute loadings:
# file_dir1 <- paste("examples/multitrait_multimethod/",
#                    model_name, " (loadings).pdf", sep = "")
# pdf(file_dir1, width = 8, height = 6)
par(mar = c(5.1, 4.1, 4.1, 4.1))
plot(abs(lambda),
     breaks = c(0, 0.20, 0.50, 1),
     main = "Absolute values of factor loadings",
     xlab = "Factors", ylab = "Items",
     col = c("black", "orange", "yellow"))
# dev.off()

# Plot absolute factor correlations:
# file_dir2 <- paste("examples/multitrait_multimethod/",
#                    model_name, " (correlations).pdf", sep = "")
# pdf(file_dir2, width = 8, height = 6)
par(mar = c(5.1, 4.1, 4.1, 4.1))
plot(abs(psi),
     breaks = c(0, 0.15, 0.30, 1),
     main = "Absolute value of factor correlations",
     xlab = "Factors", ylab = "Factors",
     col = c("black", "orange", "yellow"))
# dev.off()

#### Save results in an excel file ####

# create_xlsx(lambda, psi, digits = 2L,
#             file = paste("examples/multitrait_multimethod/", 
#                          model_name, ".xlsx", sep = ""))
