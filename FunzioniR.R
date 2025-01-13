marginTTest <- function(X, y, NullDistr = c("theoretical", "permutation"), 
                        alpha = 0.05, B = 500) {
    this.call <- match.call()
    NullDistr <- match.arg(NullDistr)
    if (NullDistr == "theoretical") {
        print("Testing using Theoretical Null Distribution")
        marginTTest.out <- pbapply(X, 2L, marginTTest.theodistr, y = y, 
                                   alpha = alpha)
        }
    else {
        print("Testing using Permutation Null Distribution")
        marginTTest.out <- pbapply(X, 2L, marginTTest.permdistr, y = y, 
                                 alpha = alpha, B = B)
        
        }
    marginTTest.out <- data.frame(Var = colnames(X), t(marginTTest.out))
    Decision <- rep("reject H0", times = ncol(X))
    Decision[marginTTest.out$p >= alpha] <- "do not reject H0"
    marginTTest.out$Decision <- Decision
    marginTTest.out <- as_tibble(marginTTest.out)
    out <- list(call = this.call, Summary = marginTTest.out)
    class(out) <- "marginTTest"
    out
}

marginTTest.theodistr <- function(x, y, alpha){
    out <- vector(mode = "numeric", length = 4L)
    names(out) <- c("TTest", "cv1", "cv2", "p-val")
    n <- length(y)
    X0 <- x[y == "control"]
    X1 <- x[y == "cancer"]
    out.ttest <- t.test(X0, X1)
    out["TTest"] <- out.ttest$statistic
    out["cv1"] <- qt(p = alpha / 2, df = n - 1)
    out["cv2"] <- qt(p = 1- alpha / 2, df = n - 1)
    out["p-val"] <- out.ttest$p.value
    out
}

marginTTest.permdistr <- function(x, y, alpha, B){
    out <- vector(mode = "numeric", length = 4L)
    names(out) <- c("TTest", "cv1", "cv2", "p-val")
    n <- length(y)
    X0 <- x[y =="control"]
    X1 <- x[y =="cancer"]
    out.ttest <- t.test(X0, X1)
    out["TTest"] <- out.ttest$statistic
    TStat.perm <- vector(mode = "numeric", length = B)
    for(b in seq_len(B)) {
        id <- sample(n)
        y.perm <- y[id]
        X0.perm <- x[y.perm == "control"]
        X1.perm <- x[y.perm == "cancer"]
        out.ttest.perm <- t.test(X0.perm, X1.perm)
        TStat.perm[b] <- out.ttest.perm$statistic
    }
    out[c("cv1", "cv2")] <- quantile(TStat.perm, probs = c(alpha / 2, 1 - alpha / 2))
    out["p-val"] <- mean(abs(TStat.perm) >= abs(out["TTest"]))
    out
}

########################################
# Definizione funzioni metodo

print.marginTTest <- function(x, ...) {
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
    cat("\nOuput della Procedura di Screening\n")
    cat("Numero Test Eseguiti = ", nrow(x$Summary), "\n\n")
    print(x$Summary, ...)
}

hist.marginTTest <- function(x, ...) {
    x$Summary %>% pull(p.val) -> p.val
    hist(p.val, ...)
}














