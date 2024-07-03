# sink("bias_29_04_2024.txt")

Évaluation du biais de l'imputation sous MAR avec le package mice dans le contexte


```{r}
ti <- Sys.time()

library(parallel)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(mice)
```


# Simulation



```{r}
simulation_no_hierachical <- function(n,length_x1, length_x2, formula_x1, formula_y, para_x1, betas) {
  x2 <- factor(sample(1:length_x2, n, replace = TRUE, prob = rep(0.25,4)))
  eta_x1 <- model.matrix(as.formula(formula_x1), data.frame(x2 = x2))[,-1] %*% as.matrix(para_x1)
  prob_x1 <-  t(apply(eta_x1, 1, function(x) exp(x - max(x)) / sum(exp(x - max(x)))))
  x1 <- as.factor(apply(prob_x1, 1, function(p) sample(1:length(p), 1, prob = p)))
  eta_y <- model.matrix(as.formula(formula_y), data = data.frame(x1 = x1, x2 = x2)) %*% betas
  prob_y <- plogis(eta_y)
  y <- as.factor(rbinom(n, 1, prob_y))
  
  return(data.frame(id = seq_len(n), y = y, x1 = x1, x2 = x2))
}
```




```{r}
set.seed(10000)
nsim <- 1000
n <- 2000

length_x1 <- 5
length_x2 <- 4

formula_x1 <- "~ x2"
para_x1 <- matrix(c(rep(2,3), rep(1,9),rep(2.5,3) ),ncol = 5)


formula_y <- "~ x1 + x2  "
betas <- c(-1.5, 1, -2, 1.5, 2, 2,  1, 2)
```





```{r}
datas <- lapply(1:nsim, function(simulation) {
  simulation_no_hierachical(n,length_x1, length_x2, formula_x1, formula_y, para_x1, betas)
})
```










## Données manquantes sur la variable ordinale (MNAR)


```{r}
simulate.mnar <- function(data, repY, id, ord_var, seed, A, B) {
  set.seed(seed)
  id.A <- data[[id]][data[[repY]] %in% 1 & data[[ord_var]] %in% A]
  sample.A <- (sample(id.A, size = round(length(id.A) * 0.3), replace = FALSE))
  data[[ord_var]][sample.A] <- NA
  id.B <- data[[id]][data[[repY]] %in% 0 & data[[ord_var]] %in% B]
   sample.B <- (sample(id.B, size = round(length(id.B) * 0.3), replace = FALSE))
  data[[ord_var]][sample.B] <- NA
  return(data)
}
```




```{r}
A = 1
B = 5
seed = 10000

datas <- lapply(datas, function(data) {
  data$x1.mis <- data$x1
  data <- simulate.mnar(data = data, A = A, B = B,id = "id", repY = "y", ord_var = "x1.mis", seed = seed)
  data$x1.imp <- data$x1.mis
  return(data)
})
```



# Imputation mutliple avec le module mice


```{r}
M = 10 # nombre d'imputation

cl <- makeCluster(round(detectCores()-1)) 
clusterExport(cl, c("datas", "M")) 
mar_data <- parLapply(cl, 1:nsim, function(s) {
  mice::mice(
    data= datas[[s]][, c("y", "x1.imp", "x2")],  
    m = M,         
    maxit = 10,      
    method = c("logreg", "polyreg",  "polyreg"),  # Méthodes d'imputation 
    print = FALSE,
    seed = 10000
  )
})

stopCluster(cl)
```

```{r}
mar_data <- lapply(1:nsim, function(s) {
  lapply(1:M, function(m) {
    cbind(
      mice::complete(mar_data[[s]], m),
      datas[[s]][, c("id", "x1", "x1.mis")]
    )
  })
})
```





## Coefficients beta et intercepts zeta du modele de regression ordinale full data without missing data

```{r}
ordinal_full <- function(data) {
  ord <- MASS::polr(x1 ~ y + x2, method = "probit", Hess = TRUE, data = data)
  beta <- coef(ord)
  zeta <- ord$zeta
  sd_beta_zeta <- summary(ord)$coefficients[, "Std. Error"]
  return(list(beta = beta, zeta = zeta, sd_beta_zeta = sd_beta_zeta))
}
```


```{r}
results <- lapply(datas, ordinal_full)
```



# Extraire les moyennes des coefficients beta et zeta


```{r}
beta_full <- rowMeans(sapply(results, function(res) res$beta))
zeta_full <- rowMeans(sapply(results, function(res) res$zeta))
sd_full <- rowMeans(sapply(results, function(res) res$sd_beta_zeta))
```





## Regression ordinale des données imputées sous MAR

```{r}
ordinal_mar <- function(data,M, length_x1, seed = NULL) {
  set.seed(seed)
  ord.regression <- lapply(data, function(imp_data) {
    MASS::polr(x1.imp ~ y + x2, method = "probit", Hess = TRUE, data = imp_data)
  })
  
   beta <- sapply(ord.regression, function(ord_reg) {
    return(as.numeric(coef(ord_reg)))
  })
  
  zeta <- sapply(ord.regression, function(ord_reg) {
    as.numeric(ord_reg$zeta)
  })
  
  sd_beta_zeta <- sapply(ord.regression, function(ord_reg) {
    summary(ord_reg)$coefficients[, "Std. Error"]
  })
  
  return(list( beta = beta, zeta = zeta, sd_beta_zeta = sd_beta_zeta))
}
```





```{r}
cl <- makeCluster(round(detectCores()-1)) 
clusterExport(cl, c("ordinal_mar", "mar_data", "M", "length_x1")) 
out_ordinal_mar <- parLapply(cl, 1:nsim, function(s) {
  ordinal_mar(
    data = mar_data[[s]], 
    M = M,
    length_x1 = length_x1,
    seed = 10000)
})

stopCluster(cl) # Arret du cluster
```








### Coeficients issus du modele de regression ordinale ou X1 est MAR


```{r}
beta_mar <- matrix(NA,nrow =nsim, ncol = length(beta_full))
for(s in 1:nsim){
    beta_mar[s,] <- rowMeans(as.data.frame(out_ordinal_mar[[s]]$beta))
}
beta_mar <-colMeans(beta_mar)
names(beta_mar) <- names(beta_full)
```



```{r}
zeta_mar <- matrix(NA,nrow =nsim, ncol = length(zeta_full))
for(s in 1:nsim){
    zeta_mar[s,] <- rowMeans(as.data.frame(out_ordinal_mar[[s]]$zeta))
}
zeta_mar <-colMeans(zeta_mar)
names(zeta_mar) <-  names( zeta_full)
```


```{r}
sd_mar <- matrix(NA,nrow =nsim, ncol = length(sd_full))
for(s in 1:nsim){
    sd_mar[s,] <- rowMeans(as.data.frame(out_ordinal_mar[[s]]$sd_beta_zeta))
}
sd_mar <-colMeans(sd_mar)
names(sd_mar) <-  names( sd_full)
```




```{r}
Beta <- as.data.frame(cbind(beta_full, sd_full = sd_full[c("y1",  "x22", "x23", "x24")],beta_mar,sd_mar = sd_mar[1:4]))
Beta$Bias <- ((Beta$beta_mar - Beta$beta_full)/Beta$beta_full)*100
Beta
```

```{r}
Zeta <- as.data.frame(cbind(zeta_full, sd_full = sd_full[c("1|2", "2|3" ,"3|4", "4|5")], zeta_mar, sd_mar = sd_mar[5:8]))
Zeta$Bias <- ((Zeta$zeta_mar - Zeta$zeta_full )/Zeta$zeta_full)*100
Zeta
```





```{r}
tf <- Sys.time()
code_time = tf - ti
code_time
```


```{r}
couleur <- c("red", "black",  "blue",  "orange","green", "magenta","darkgreen", "purple", "cyan","#CCFF00")
# Fonction pour les graphes
mar_simule_graph <- function(data, title, labels, colors,xlab, ylab) {
  p <- ggplot() +
    geom_point(data = as.data.frame(data), aes(x = rownames(data), y = data[,1], color = labels[1], shape = labels[1]), size = 2) +
    geom_point(data = as.data.frame(data), aes(x = rownames(data), y = data[,2], color = labels[2], shape = labels[2]), size = 2) +
    scale_color_manual(name = "Legend", values = setNames(colors, labels)) +
    scale_shape_manual(name = "Legend", values = setNames(c(17, 19), labels)) +
    theme_clean(base_size = 7) +
    ylim(range(data)[1], range(data)[2])+
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5)) +
    ggtitle(title) +
    xlab(xlab) +      
    ylab(ylab) 
  return(p)
}
```

```{r}
beta_graph <- mar_simule_graph(data = Beta[, c("beta_full", "beta_mar")], title = "Simulated Coefficients vs. Coefficients Under MAR", 
                               labels = c( "Simule", "Mar"),colors = couleur[1:2], xlab = "Coefficients", 
                               ylab = "Estimate")
zeta_graph <- mar_simule_graph(data = Zeta[, c("zeta_full", "zeta_mar")], title = "Simulated Intercepts vs Intercepts under Mar", 
                               labels = c( "Simule", "Mar"),colors = couleur[1:2], xlab = "Intercepts", 
                               ylab = "Estimate")

grid.arrange(grobs = list(beta_graph, zeta_graph), ncol = 2)
```

```{r}
#sink()
```
