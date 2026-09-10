library(latent)
library(fungible)
library(plot.matrix)

# Data from https://osf.io/72zp3/

#### Extract the data ####

# Check samples and sample sizes
samples <- unique(icases$sample) # industry mooc fire student dutch
Ns <- sapply(samples, FUN = function(x) sum(icases$sample == x))
names(Ns) <- samples
Ns

# Subset of items pertaining to the HEXACO-100:
selection <- 5:104
full <- icases[, selection]
# Reorder the dataframe:
reorder <- c("hexemfea146" , "hexemfea170" , "hexemfea74" , "hexemfea2" ,
             "hexemanx128" , "hexemanx8" , "hexemanx80" , "hexemanx176" ,
             "hexemdep62" , "hexemdep182" , "hexemdep134" , "hexemdep158" ,
             "hexemsen44" , "hexemsen164" , "hexemsen20" , "hexemsen68" ,

             "hexexsse3" , "hexexsse99" , "hexexsse51" , "hexexsse171" ,
             "hexexsbo33" , "hexexsbo57" , "hexexsbo177" , "hexexsbo81" ,
             "hexexscb111" , "hexexscb135" , "hexexscb63" , "hexexscb39" ,
             "hexexliv165" , "hexexliv93" , "hexexliv69" , "hexexliv45" ,

             "hexcoorg125" , "hexcoorg173" , "hexcoorg29" , "hexcoorg53" ,
             "hexcodil35" , "hexcodil107" , "hexcodil11" , "hexcodil131" ,
             "hexcoper17" , "hexcoper41" , "hexcoper65" , "hexcoper89" ,
             "hexcopru191" , "hexcopru95" , "hexcopru47" , "hexcopru71" ,

             "hexopaes6" , "hexopaes54" , "hexopaes78" , "hexopaes150" ,
             "hexopinq60" , "hexopinq180" , "hexopinq132" , "hexopinq12" ,
             "hexopcre162" , "hexopcre138" , "hexopcre186" , "hexopcre42" ,
             "hexopunc48" , "hexopunc144" , "hexopunc120" , "hexopunc192" ,

             "hexagfor148" , "hexagfor172" , "hexagfor4" , "hexagfor52" ,
             "hexaggen154" , "hexaggen106" , "hexaggen130" , "hexaggen82" ,
             "hexagfle112" , "hexagfle16" , "hexagfle184" , "hexagfle160" ,
             "hexagpat46" , "hexagpat142" , "hexagpat94" , "hexagpat70" ,

             "hexhosin1" , "hexhosin121" , "hexhosin73" , "hexhosin97" ,
             "hexhofai127" , "hexhofai7" , "hexhofai175" , "hexhofai79" ,
             "hexhogre109" , "hexhogre85" , "hexhogre37" , "hexhogre157" ,
             "hexhomod43" , "hexhomod139" , "hexhomod187" , "hexhomod67")
full <- full[, reorder] # 96 items

# Extract the MOOC sample:
mooc <- full[icases$sample == samples[2], ]
dim(mooc) # 4286   96

#### Create target matrices ####

Target_g <- rbind(diag(6) %x% rep(1, 16))
Target_s <- diag(24) %x% rep(1, 4)
Target <- cbind(Target_g, Target_s)
PsiTarget <- diag(30); PsiTarget[1:6, 1:6] <- 1; diag(PsiTarget) <- 0

plot(Target, breaks = c(0, 0.5, Inf),
     main = "Target on factor loadings",
     col = c("black", "yellow"))
plot(PsiTarget, breaks = c(0, 0.5, Inf),
     main = "Target on factor correlations",
     col = c("black", "yellow"))

#### Bi-Factor Analysis POB (Target; 6 general + 24 specific factors) ####

rotation_model <- list(
  list(geomin = list(epsilon = 0.01, factors = 1:6),
       target = list(target = Target[, 7:30], factors = 7:30)),
  list(oblimin = list(gamma = 0, factors = 1:6),
       target = list(target = Target[, 7:30], factors = 7:30)),
  list(geomin = list(epsilon = 0.01, factors = 1:6),
       geomin = list(epsilon = 0.01, factors = 7:30)),
  list(oblimin = list(gamma = 0, factors = 1:6),
       oblimin = list(gamma = 0, factors = 7:30))
)

names(rotation_model) <- c("geomin + target", "oblimin + target",
                           "geomin + geomin", "oblimin + oblimin")
nmodels <- length(rotation_model)

set.seed(2026)
for(i in seq_len(nmodels)) {
  
  efa <- lefa(data = mooc, nfactors = 30,
              ordered = TRUE, estimator = "dwls",
              projection = "poblq", oblique = 6,
              rotation = rotation_model[[i]],
              control.efa = list(maxit = 5000L, rstarts = 3L, cores = 3L),
              control.rotation = list(rstarts = 10L, cores = 10L),
              se = FALSE)
  
  # Align the results to the theoretical Target:
  lambda <- efa@transformed_pars$lambda_rotated
  psi <- efa@transformed_pars$psi_rotated; diag(psi) <- 1
  x <- fungible::faAlign(F1 = Target, F2 = abs(lambda), Phi2 = abs(psi))
  lambda <- lambda[, c(x$FactorMap[2, ])]
  psi <- psi[c(x$FactorMap[2, ]), c(x$FactorMap[2, ])]
  
  # Dimnames:
  rownames(lambda) <- NULL
  colnames(lambda) <- c(paste("G", 1:6, sep = ""), paste("S", 1:24, sep = ""))
  rownames(psi) <- colnames(psi) <-
    c(paste("G", 1:6, sep = ""), paste("S", 1:24, sep = ""))
  
  # Save results in an excel file:
  model_name <- names(rotation_model)[i]
  create_xlsx(lambda, psi, digits = 2L,
              file = paste("examples/correlated_general_factors/", 
                           model_name, ".xlsx", sep = ""))
  
  # Plot absolute loadings:
  subfix <- paste("(", model_name, ")", sep = "")
  file_dir1 <- paste("examples/correlated_general_factors/",
                     model_name, " (loadings).pdf", sep = "")
  pdf(file_dir1, width = 8, height = 6)
  par(mar = c(5.1, 4.1, 4.1, 4.1))
  plot(abs(lambda),
       breaks = c(0, 0.20, 0.50, 1),
       main = paste("Absolute value of factor loadings", 
                    subfix, collapse = " "),
       xlab = "Factors", ylab = "Items",
       col = c("black", "orange", "yellow"))
  dev.off()
  
  # Plot absolute factor correlations:
  file_dir2 <- paste("examples/correlated_general_factors/",
                     model_name, " (correlations).pdf", sep = "")
  pdf(file_dir2, width = 7, height = 6)
  par(mar = c(5.1, 4.1, 4.1, 4.1))
  plot(abs(psi),
       breaks = c(0, 0.15, 0.30, 1),
       main = paste("Absolute value of factor correlations", 
                    subfix, collapse = " "),
       xlab = "Factors", ylab = "Factors",
       col = c("black", "orange", "yellow"))
  dev.off()
  
}
