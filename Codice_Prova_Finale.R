################################################################################
# Codice prova finale Statistica
# Studente: Fonzino Adriano
################################################################################

library("tibble")
library("dplyr")
source("FunzioniR.R")
library(pbapply)
library(multtest)
library(glmnet)

Prostata <- read.table(file = "Prostata.txt", header = TRUE)
dim(Prostata)


# estrazione variabile di risposta
y <- as.factor(Prostata$Status)

# estrazione variabili esplicative
X <- as.matrix(Prostata[, -1])

rm(Prostata)

alpha <- 0.05

################################################################################
######################### SCREENING PROCEDURES (MTP) ###########################
################################################################################

######################
# THEORETICAL METHOD #
######################

out.ttest.theoretical <- marginTTest(X = X, y = y, NullDistr = "theoretical")
out.ttest.theoretical
hist(out.ttest.theoretical)


out.ttest.theoretical$Summary %>% arrange(p.val) %>% pull(p.val) %>% 
  mt.rawp2adjp(proc = c("Bonferroni", "SidakSS", "Holm", "SidakSD", "Hochberg", 
                        "BH", "BY")) -> adj_pval.theor


Decision.theor <- ifelse(adj_pval.theor$adjp <= alpha, "reject Ho", "do not reject H0")

out.ttest.theoretical$Summary %>% arrange(p.val) -> MTP.Results.theor
for (i in c("Bonferroni", "SidakSS", "Holm", "SidakSD", "Hochberg", 
            "BH", "BY")) {
  MTP.Results.theor[i] = as.tibble(Decision.theor)[i]
}
MTP.Results.theor

MTP.Results.theor %>% select(Decision:BY) %>% 
  apply(2, table) -> MTP.Summary.05.theor
MTP.Summary.05.theor


tibble(Procedure = colnames(MTP.Summary.05.theor),
       NoRej_H0 = MTP.Summary.05.theor[1L, ],
       Rej_H0 = MTP.Summary.05.theor[2L, ]) -> MTP.Summary.05.theor
MTP.Summary.05.theor
write.table(MTP.Summary.05.theor , file = "immagini_tabelle\\1_Summary_theor.csv")


######################
# PERMUTATION METHOD #
######################

source("FunzioniR.R")

out.ttest.permutation <- marginTTest(X = X, y = y, NullDistr = "permutation", 
                                     B = 5000)
out.ttest.permutation
hist(out.ttest.permutation)


out.ttest.permutation$Summary %>% arrange(p.val) %>% pull(p.val) %>% 
  mt.rawp2adjp(proc = c("Bonferroni", "SidakSS", "Holm", "SidakSD", "Hochberg", 
                        "BH", "BY")) -> adj_pval.perm


Decision.perm <- ifelse(adj_pval.perm$adjp <= alpha, "reject Ho", "do not reject H0")

out.ttest.permutation$Summary %>% arrange(p.val) -> MTP.Results.perm
for (i in c("Bonferroni", "SidakSS", "Holm", "SidakSD", "Hochberg", 
            "BH", "BY")) {
  MTP.Results.perm[i] = as.tibble(Decision.perm)[i]
}
MTP.Results.perm

MTP.Results.perm %>% select(Decision:BY) %>% 
  apply(2, table) -> MTP.Summary.05.perm
MTP.Summary.05.perm


tibble(Procedure = colnames(MTP.Summary.05.perm),
       NoRej_H0 = MTP.Summary.05.perm[1L, ],
       Rej_H0 = MTP.Summary.05.perm[2L, ]) -> MTP.Summary.05.perm
MTP.Summary.05.perm
write.table(MTP.Summary.05.perm , file = "immagini_tabelle\\2_Summary_perm.csv")

# Standardize data
X.scaled <- scale(X)

# divide into train and test sets
set.seed(1234)
id1 <- sample(which(y == "control"), size = 10)
id2 <- sample(which(y == "cancer"), size = 10) 

Xtrain <- X.scaled[-c(id1,id2), ]
ytrain <- y[-c(id1, id2)]
dim(Xtrain)
length(ytrain)

Xtest <- X.scaled[c(id1,id2),]
ytest <- y[c(id1, id2)]
dim(Xtest)  
length(ytest)


######################################
# REGRESSION MODEL WITHOUT SCREENING #
######################################

# selects significant genes without screening techniques

MTP.Results.perm %>% select(Var, Decision) %>% 
  filter(Decision == "reject H0") %>% pull(Var) -> genes.sign05.perm
genes.sign05.perm


Xtrain.sign05.perm <- Xtrain[, genes.sign05.perm]
colnames(Xtrain.sign05.perm) == genes.sign05.perm

# p = 481 ==> p >> n ==> lasso selection!
obj.lasso <- glmnet(Xtrain.sign05.perm, ytrain, family = "binomial", alpha = 1)
names(obj.lasso)

plot(obj.lasso, xvar = "lambda", xlab = expression(log(lambda)),
     main = "Logistic Regression with lasso Penalization",
     lwd = 2, cex.lab = 1.6)

# select best lambda value by Cross Validation
out.lasso.cv <- cv.glmnet(x = Xtrain.sign05.perm, y = ytrain, 
                          family = "binomial", type.measure = "mse")
out.lasso.cv

plot(out.lasso.cv, xvar = "lambda", xlab = expression(log(lambda)), 
     main = "High-Dimesional Lasso Regression Model", lwd = 2, cex.lab = 1.6)

names(out.lasso.cv)
out.lasso.cv$nzero

# get lambda with the min. MSE
out.lasso.cv$lambda.min

# extract informations about the non zero coefficients selected from CV
idNonZero <- predict(out.lasso.cv, newx = XnoNA, type = "nonzero", s = "lambda.min")
lasso.coeff <- predict(out.lasso.cv, newx = XnoNA, type = "coef", s = "lambda.min")
lasso.genes <- names(lasso.coeff[lasso.coeff[,1] != 0, ])[-1]
lasso.genes

# select genes obtained from lasso feature selection
Xtrain.sign05.perm.lasso <- Xtrain.sign05.perm[, lasso.genes]
dim(Xtrain.sign05.perm.lasso)



# Ridge Regression for Predictive Model

out.ridge.cv <- cv.glmnet(x = Xtrain.sign05.perm.lasso, y = ytrain, alpha = 0, 
                          family = "binomial", type.measure = "mse", keep=TRUE)
plot(out.ridge.cv, main = "Linear Regression Model with Ridge Penalty")
out.ridge.cv$lambda.min

coef(out.ridge.cv, s = "lambda.min")

final.ridge.model <- glmnet(x = Xtrain.sign05.perm.lasso, y = ytrain, 
                            family = "binomial", 
                            lambda = out.ridge.cv$lambda.min, alpha = 0)

yhat.ridge.cv <- predict(final.ridge.model, newx = Xtest[, lasso.genes], 
                         type= "response", s = "lambda.min")
yhat.ridge.cv

yhat.ridge.cv.discrete <- ifelse(yhat.ridge.cv > 0.5, "control", "cancer")
#View(cbind(as.factor(yhat.ridge.cv.discrete), ytest))


confusion.glmnet(final.ridge.model, 
                 newx = Xtest[, lasso.genes], newy = ytest, 
                 family = "binomial") -> confusion_05
confusion_05

write.table(confusion_05 , file = "immagini_tabelle\\3d_confusionM_05.csv")

#######################################
# REGRESSION MODEL WITH MTP SCREENING #
#######################################

MTP.Summary.05.perm

# Selecting genes from BH MTP screening
MTP.Results.perm %>% select(Var, BH) %>% 
  filter(BH == "reject Ho") %>% pull(Var) -> genes.BH.perm
genes.BH.perm

Xtrain.BH.perm <- Xtrain[, genes.BH.perm]
colnames(Xtrain.BH.perm) == genes.BH.perm

# p = 35 ==> lasso selection!
obj.lasso.bh <- glmnet(Xtrain.BH.perm, ytrain, family = "binomial", alpha = 1)
names(obj.lasso.bh)

plot(obj.lasso.bh, xvar = "lambda", xlab = expression(log(lambda)),
     main = "Logistic Regression with lasso Penalization on BH genes",
     lwd = 2, cex.lab = 1.6)

# select best lambda value by Cross Validation
out.lasso.cv.bh <- cv.glmnet(x = Xtrain.BH.perm, y = ytrain, 
                          family = "binomial", type.measure = "mse")
out.lasso.cv.bh

plot(out.lasso.cv.bh, xvar = "lambda", xlab = expression(log(lambda)), 
     main = "High-Dimesional Lasso Regression Model", lwd = 2, cex.lab = 1.6)

names(out.lasso.cv.bh)
out.lasso.cv.bh$nzero

# get lambda with the min. MSE
out.lasso.cv.bh$lambda.min

# extract information about the non zero coefficients selected from CV
idNonZero.bh <- predict(out.lasso.cv.bh, newx = XnoNA, type = "nonzero", s = "lambda.min")
lasso.coeff.bh <- predict(out.lasso.cv.bh, newx = XnoNA, type = "coef", s = "lambda.min")
lasso.genes.bh <- names(lasso.coeff.bh[lasso.coeff.bh[,1] != 0, ])[-1]
lasso.genes.bh
length(lasso.genes.bh)

#We selected 27 genes from 35.
# select genes obtained from lasso feature selection
Xtrain.bh.perm.lasso <- Xtrain.BH.perm[, lasso.genes.bh]
dim(Xtrain.bh.perm.lasso)



# Ridge Regression for Predictive Model
out.ridge.cv.bh <- cv.glmnet(x = Xtrain.bh.perm.lasso, y = ytrain, alpha = 0, 
                          family = "binomial", type.measure = "mse", keep=TRUE)
plot(out.ridge.cv.bh, main = "Linear Regression Model with Ridge Penalty")
out.ridge.cv.bh$lambda.min

coef(out.ridge.cv.bh, s = "lambda.min")

final.ridge.model.bh <- glmnet(x = Xtrain.bh.perm.lasso, y = ytrain, 
                            family = "binomial", 
                            lambda = out.ridge.cv.bh$lambda.min, alpha = 0)

yhat.ridge.cv.bh <- predict(final.ridge.model.bh, newx = Xtest[, lasso.genes.bh], 
                         type= "response", s = "lambda.min")
yhat.ridge.cv.bh

yhat.ridge.cv.bh.discrete <- ifelse(yhat.ridge.cv.bh > 0.5, "control", "cancer")
View(cbind(as.factor(yhat.ridge.cv.bh.discrete), ytest))


confusion.glmnet(final.ridge.model.bh, 
                 newx = Xtest[, lasso.genes.bh], newy = ytest, 
                 family = "binomial") -> confusion_bh
confusion_bh


write.table(confusion_bh , file = "immagini_tabelle\\4d_confusionM_bh.csv")


# The CV-Ridge Logistic Regression model built through the use of the genes 
# screened firstly by MTP with BH's correction (35 significant expression levels)
# and after by Cross Validated Lasso Logistic Regression's features selection, 
# gave back an accuracy of 100%. 
# The same model, built without the MTP screening of the genes (481 significant 
# genes' expression levels) has been less accurate, with an overall accuracy of 
# 75%, misclassifying 3 control and 2 cancer samples, showing thus, an overfitting 
# (higher variance) with respect to the less complicated model, trained on a 
# minor number of features.

as.matrix(coef(final.ridge.model.bh)) %>% 
  row.names()


tibble(Coeff = as.matrix(coef(final.ridge.model.bh))[,1]) -> coeff
coeff["Gene Name"] = as.matrix(coef(final.ridge.model.bh)) %>% 
  row.names()

write.csv(coeff, file="immagini_tabelle\\5_final_model_coeff.csv")

coeff <- coeff[-1,]
coeff %>% arrange(desc(Coeff)) -> coeff
as.matrix(coeff$Coeff)
barplot(t(as.matrix(coeff$Coeff)), las=2, ylab = "abs(Coefficent)")

# single box plots
for (gene in lasso.genes.bh) {
  png(filename = paste("immagini_tabelle\\boxplots\\", gene, ".png", sep = ""))
  boxplot(X[, gene] ~ y, main= gene, ylab = "Gene Expression Level", 
          xlab = "")
  dev.off()
}

# first graph
par(mfrow=c(2,7))
for (gene in lasso.genes.bh[1:14]) {
  boxplot(X[, gene] ~ y, main= gene, ylab = "Gene Expression Level", 
          xlab = "")
}

# second graph
par(mfrow=c(2,7))
for (gene in lasso.genes.bh[15:27]) {
  boxplot(X[, gene] ~ y, main= gene, ylab = "Gene Expression Level", 
          xlab = "")
}



