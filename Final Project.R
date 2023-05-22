prior <- function(betas){
    return(prod(sapply(betas, dnorm, mean=0, sd=10)))
}

likelihood <- function(X, y, betas){
    e <- exp(X %*% betas)
    px <- e/(1+e)
    return(px^y * ((1-px)^(1-y)))
}

log_posterior <- function(X, y, betas){
    log_prior <- log(prior(betas))
    log_likelihoods <- log(likelihood(X, y, betas))
    return(log_prior + sum(log_likelihoods))
}

MCMC <- function(X, y, n, beta_start=c(0,0,0), jump_dist_sd=0.1){
       B <- ncol(X)   # number of betas (coefficients)
       beta <- matrix(nrow=n, ncol=B)
       beta[1,] <- beta_start
       for(i in 2:n){
             current_betas <- beta[i-1,]
             new_betas <- current_betas + rnorm(B, mean=0, sd=jump_dist_sd)
             for(j in 1:B){
                   test_betas <- current_betas
                   test_betas[j] <- new_betas[j]
                   rr <- log_posterior(X, y, test_betas) - log_posterior(X, y, current_betas)
                   if(log(runif(1)) < rr){
                         beta[i,j] <- new_betas[j]
                   } 
                   else {
                           beta[i,j] <- current_betas[j]
                   }
              }
       }
     return(beta)
}

# data <- read.csv("ai4i2020.csv")

burn_in_length <- 10000
N <- 50000   # total number of points
y <- data$Machine.failure
X <- cbind(1, as.matrix(data[,c(7,8)]))   # pad with a column of ones for intercept

beta <- MCMC(X, y, n=N, jump_dist_sd=0.05)


beta2 <- beta[(burn_in_length+1):N,]

MCMC_results <- apply(beta2, MARGIN=2, function(X){c(mean(X), quantile(X, c(0.025, 0.5, 0.975)))})
colnames(MCMC_results) <- c("Intercept", "Torque..Nm.", "Tool.wear..min.")

glm_model <- glm(Machine.failure ~ Torque..Nm. + Tool.wear..min., data=data, family=binomial(link="logit"))

sc <- summary(glm_model)$coefficients
glm_95 <- rbind(sc[,1] - 1.96*sc[,2], sc[,1] + 1.96*sc[,2])

g <- ggplot(data=NULL, aes(x=1:N)) + xlab("") + geom_vline(xintercept=burn_in_length, color="red")
p_int <- g + geom_line(aes(y=beta[,1])) + ylab("Intercept")
p_tor <- g + geom_line(aes(y=beta[,2])) + ylab("Torque..Nm.")
p_twe <- g + geom_line(aes(y=beta[,3])) + ylab("Tool.wear..min.")

ggarrange(p_int, p_tor, p_twe, ncol = 2, nrow = 2)

g2 <- ggplot(data=NULL) + ylab("") + theme_minimal()
p2_int <- g2 + geom_histogram(aes(x=beta2[,1])) + geom_vline(xintercept=c(MCMC_results[2,1], MCMC_results[4,1])) + geom_vline(xintercept=glm_95[,1], linetype=2, color = "blue") + xlab("Value") + ggtitle("Intercept")
p2_tor <- g2 + geom_histogram(aes(x=beta2[,2])) + geom_vline(xintercept=c(MCMC_results[2,2], MCMC_results[4,2])) + geom_vline(xintercept=glm_95[,2], linetype=2, color = "blue") + xlab("Coefficient") + ggtitle("Torque")
p2_twe <- g2 + geom_histogram(aes(x=beta2[,3])) + geom_vline(xintercept=c(MCMC_results[2,3], MCMC_results[4,3])) + geom_vline(xintercept=glm_95[,3], linetype=2, color = "blue") + xlab("Coefficient") + ggtitle("Tool wear")
ggarrange(p2_int, p2_tor, p2_twe, ncol = 2, nrow = 2)
