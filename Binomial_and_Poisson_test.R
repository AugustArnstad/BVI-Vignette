library(lme4) # for simulating data
library(INLA)

# If not already installed, install the 'devtools' package
if(!require(devtools)) install.packages("devtools")
# Install VariableImportanceINLA
devtools::install_github("AugustArnstad/VariableImportanceINLA")
library(VariableImportanceINLA)
## BINOMIAL

# Parameters
n = 1000
n_groups <- 50
group_effect <- rnorm(n_groups, mean = 0, sd = 1)
group <- rep(1:n_groups, each = n / n_groups)

# Simulate fixed effect
x1 <- rnorm(n, mean = 0, sd = 1)
x2 <- rnorm(n, mean = 0, sd = 1)
fixed_effect <- 1 # effect size for x

# Simulate latent variable z using the probit model approach
eta <- fixed_effect * x1 + x2 + group_effect[group]

probit_p <- pnorm(eta)

# Simulate binomial response based on these probabilities
y_bin <- rbinom(n, size = 1, prob = probit_p)

# Update the data frame
data_binomial <- data.frame(y = y_bin, x1 = x1, x2=x2, group = as.factor(group))

# Probability of success
#logit_p <- fixed_effect * x + group_effect[group]
#p <- exp(logit_p) / (1 + exp(logit_p))
# Simulate binomial response
#y <- rbinom(n, size = 1, prob = p) # assuming a single trial per observation


# Convert z to binary outcome y using threshold (e.g., 0)
# This simulates the binary outcome based on the probit model assumption

formula_binomial <- y ~ x1 + x2 + f(group, model = "iid")

# Fit the model using INLA
result_binomial <- inla(formula_binomial,
                        family = "binomial",
                        data = data_binomial,
                        control.family = list(link = "probit"))

# Summary of the results
summary(result_binomial)


result_VII_binomial <- VariableImportanceINLA::perform_inla_analysis(data_binomial, formula_binomial, family = "binomial", "probit")


sample_posterior_count <- function(model, formula, data, n_samp=10, additive_param=NULL, param_of_interest=NULL) {
  
  # Make sure its correct with distributional variance
  # Make the param_of_interest a general input object for all functions. That was a nice way to put it
  # Scaling not complete
  
  
  response <- all.vars(formula)[1]
  scaled_response <- scale(data[, response])
  scale_const <- attr(scaled_response, 'scaled:scale')
  
  #fixed <- model$names.fixed[2:length(model$names.fixed)]
  effects <- extract_effects(formula)
  fixed <- effects$fixed_effects
  
  fam <- model$.args$family
  
  link <- model$.args$control.family[[1]]$link
  
  distribution = paste0("Distributional variance: ", link)
  
  response <- all.vars(formula)[1]
  
  variance_marginals_list <- lapply(model$marginals.hyperpar, function(x) inla.tmarginal(function(t) 1/t, x))
  random <- names(variance_marginals_list)
  random <- gsub("Precision for the ", "", random)
  random <- gsub("Precision for ", "", random)
  
  random <- c(distribution, random)
  
  beta_mat <- matrix(NA, nrow=n_samp, ncol=length(fixed))
  scaled_beta_mat <- matrix(NA, nrow=n_samp, ncol=length(fixed))
  importance_mat <- matrix(NA, nrow=n_samp, ncol=length(fixed))
  scaled_importance_mat <- matrix(NA, nrow=n_samp, ncol=length(fixed))
  names(importance_mat) <- fixed
  R2_mat <- matrix(NA, nrow=n_samp, ncol=1)
  R2_cond_mat <- matrix(NA, nrow=n_samp, ncol=1)
  var_pred_mat <- matrix(NA, nrow=n_samp, ncol=1)
  h2_mat <- matrix(NA, nrow=n_samp, ncol=1)
  random_mat <- matrix(NA, nrow=n_samp, ncol=length(random))
  scaled_random_mat <- matrix(NA, nrow=n_samp, ncol=length(random))
  
  if(!is.null(fixed)){
    SVD <- VariableImportanceINLA::SVD_decomp(data[fixed])
    
    lambda <- SVD$lambda
  }else {
    lambda <- diag(length(fixed))
  }
  
  samps_Z <- inla.posterior.sample(model, n = n_samp)
  
  latent_row_names <- rownames(samps_Z[[1]]$latent)
  
  output_length=length(samps_Z[[1]]$latent)
  
  not_na = which(!is.na(data[response]))
  
  
  for (i in 1:n_samp){
    # Extract all sampled values, separate them by covariate/predictor, and assign them to the sampled matrix
    samples_length <- 0
    predictor <- paste0("^Predictor:")
    predictor_names <- grep(predictor, latent_row_names, value = TRUE)
    predictor_samples <- samps_Z[[i]]$latent[predictor_names, , drop = FALSE]
    samples_tot <- length(predictor_samples)
    
    if (fam == "binomial"){
      if (link == "probit"){
        distribution_var <- 1
      } else if (link == "logit"){
        distribution_var <- pi^2/3
      }
    }else if (fam == "poisson"){
      if (link == "log"){
        intercept <- samps_Z[[i]]$latent[output_length-length(fixed)]
        distribution_var <- log((1/exp(intercept)) + 1)
      }else if (link == "root"){
        distribution_var <- 0.25
      }
    }
    
    
    random_mat[i, 1] <- distribution_var
    
    
    if (length(random)>1){
      for (j in 2:length(random)){
        pattern <- paste0("^", random[j], ":")
        random_names <- grep(pattern, latent_row_names, value = TRUE)
        random_samples <- samps_Z[[i]]$latent[random_names, , drop = FALSE]
        samples_tot <- samples_tot + length(random_samples)
        random_mat[i, j] <- var(random_samples)
      }
    }else{
      if (i==n_samp){
        print("No random effects, only resiudal variance")
      }
      
    }
    
    beta <- samps_Z[[i]]$latent[(samples_tot+2):output_length]  #Skip intercept
    #print(beta)
    beta_mat[i, ] <- beta
    #scaled_beta <- beta/scale_const
    #scaled_beta_mat[i, ] <- scaled_beta
    importance_mat[i, ] <- lambda^2 %*% beta^2
    #scaled_importance_mat[i, ] <- lambda^2 %*% scaled_beta^2
    
    
  }
  
  rowsum <- rowSums(random_mat) + rowSums(importance_mat) + distribution_var
  scaled_random_mat <- random_mat/rowsum
  scaled_beta_mat <- beta_mat/rowsum
  scaled_importance_mat <- importance_mat/rowsum
  
  # Do not think these are correct!!!!! They do not contain the guassian observations for example.
  # Could possibly use the variance of the predictor as the measure
  
  R2_mat <- rowSums(importance_mat) / (rowSums(importance_mat) + rowSums(random_mat) + distribution_var)
  
  if (length(random)>1){
    R2_cond <- (rowSums(importance_mat) + rowSums(random_mat)) / (rowSums(importance_mat) + rowSums(random_mat) + distribution_var)
  }else if (length(random)==1){
    R2_cond <- (rowSums(importance_mat) + random_mat ) / (rowSums(importance_mat) + rowSums(random_mat) + distribution_var)
  }else{
    R2_cond <- rowSums(importance_mat)  / (rowSums(importance_mat) + rowSums(random_mat) + distribution_var)
  }
  
  
  beta_mat <- as.data.frame(beta_mat)
  names(beta_mat) <- fixed
  importance_mat <- as.data.frame(importance_mat)
  names(importance_mat) <- fixed
  scaled_beta_mat <- as.data.frame(scaled_beta_mat)
  names(scaled_beta_mat) <- fixed
  scaled_importance_mat <- as.data.frame(scaled_importance_mat)
  names(scaled_importance_mat) <- fixed
  
  random_mat <- as.data.frame(random_mat)
  names(random_mat) <- random
  scaled_random_mat <- as.data.frame(scaled_random_mat)
  names(scaled_random_mat) <- random
  
  R2_mat <- as.data.frame(R2_mat)
  names(R2_mat) <- "Marginal R2"
  R2_cond <- as.data.frame(R2_cond)
  names(R2_cond) <- "Conditional R2"
  
  if (!is.null(additive_param)){
    h2_mat <- random_mat[, additive_param]/(rowSums(importance_mat) + rowSums(random_mat) + distribution_var)
    h2_mat <- as.data.frame(h2_mat)
    names(h2_mat) <- paste0("Heritability of: ", additive_param)
  }
  
  
  return(list(beta_samples = beta_mat,
              importance_samples = importance_mat,
              scaled_beta_samples = scaled_beta_mat,
              scaled_importance_samples = scaled_importance_mat,
              random_samples = random_mat,
              scaled_random_samples = scaled_random_mat,
              R2_marginal = R2_mat,
              R2_conditional = R2_cond,
              var_y = var_pred_mat,
              heritability = h2_mat))
}



VariableImportanceINLA::sample_posterior_count(result_VII_binomial, formula_binomial, data_binomial, n_samp=5, additive_param=NULL, param_of_interest=NULL)#, distribution_var=1)

summary(result_VII_binomial)


abc <- inla.posterior.sample(result_VII_binomial, n=5)
length(abc[[1]]$latent)
abc[[1]]$latent[153-2]

## POISSON

# Simulate random effects
group_effect_pois <- rnorm(n_groups, mean = 0, sd = 1)

# Simulate fixed effect
fixed_effect_pois <- 1 # effect size for x

# Linear predictor and expected count
lambda <- exp(fixed_effect_pois * x1 + x2 + group_effect_pois[group])

# Simulate Poisson response
y_pois <- rpois(n, lambda = lambda)
# Create data frame
data_poisson <- data.frame(y = y_pois, x1 = x1, x2=x2, group = as.factor(group))

formula_poisson <- y ~ x1 + x2 + f(group, model = "iid")

# Fit the model using INLA
result_poisson <- inla(formula_poisson,
                        family = "poisson",
                        data = data_poisson,
                        control.family = list(link = "log"))

# Summary of the results
summary(result_poisson)


result_VII_poisson <- VariableImportanceINLA::perform_inla_analysis(data_poisson, formula_poisson, family = "poisson", "log")

beta_0 <- result_VII_poisson$summary.fixed$mean[1]
sigma_d <- (1/exp(beta_0)) + 1
sigma_d

VariableImportanceINLA::sample_posterior_count(result_VII_poisson, formula_poisson, data_poisson, n_samp=5, additive_param=NULL, param_of_interest=NULL)#, distribution_var=sigma_d)


summary(result_VII_poisson)


beta_test <- c(1, 2)

as.matrix(data_poisson[, c("x1", "x2")]) %*% beta_test

data_poisson


inla.list.models()

