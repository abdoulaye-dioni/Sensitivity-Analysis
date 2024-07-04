
**Bias Evaluation of Imputation under MAR Using the jomo Package in hiechical Context**


```{r}
#sink("bias_hierachique_29_04_2024.txt")
```

```{r}
ti <- Sys.time()

library(mitml)
library(parallel)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(jomo)
```



# Simulation



```{r}
simulation_hierachical <- function(length_x1, length_x2, formula_x1, formula_y, para_x1,  betas,n_clus,n_obs,sd_U) {
  x2 <- factor(sample(1:length_x2, n_clus*n_obs, replace = TRUE, prob = rep(0.25,4)))
  eta_x1 <- model.matrix(as.formula(formula_x1), data.frame( x2 = x2))[,-1] %*% as.matrix(para_x1)
  prob_x1 <-  t(apply(eta_x1, 1, function(x) exp(x - max(x)) / sum(exp(x - max(x)))))
  x1 <- ordered(apply(prob_x1, 1, function(p) sample(1:length(p), 1, prob = p)))
  clus <- rep(1:n_clus, each=n_obs)
  U <- rnorm(n_clus, mean=0, sd = sd_U)
  ZU <- U[clus]
  eta_y <- model.matrix(as.formula(formula_y), data = data.frame(x1 = x1, x2 = x2)) %*%betas + ZU
  prob_y <- plogis(eta_y)
  y <- as.factor(rbinom(n_clus*n_obs, 1, prob_y))
  
  return(data.frame(id = seq_len(n_clus*n_obs), y = y, x1 = x1,  x2 = x2, clus = clus))
}
```


## Parameters 

```{r}
set.seed(100)
nsim = 500
n_clus = 10
n_obs = 200
length_x1 <- 3
length_x2 = 4
sd_U = 0.45
A = 1
B = 3
```

```{r}
formula_x1 <- "~ x2"
para_x1 <- matrix(c(rep(2,3), rep(1,3),rep(2.5,3) ),ncol = 3)


formula_y <- "~ x1 + x2"
betas <- c(-1, 1, -2, 2,  1, 2)
```





```{r}
datas <- lapply(1:nsim, function(simulation) {
    data <- simulation_hierachical(length_x1, length_x2, formula_x1, formula_y, para_x1,  betas,n_clus,n_obs,sd_U)
  data$x1.mis <- data$x1
  return(data)
  })
```








## Missing Not At Random on the variable ordinale


```{r}
missing_prob <- data.frame(matrix(c(0.2,0.3,0.1,0.4,0.4,0.1,0.1,0.3), nrow = 2, byrow = FALSE,dimnames = list(c("A","B"),paste0("proba",1:length_x2))))
missing_prob 
```


```{r}
simulate_mnar <-  function(data, proba, cat_var) {
  
  mnar_by_strata <- function(data, repY, id, ord_var,probA, probB, A, B, seed) {
  set.seed(seed)
  id.A <- data[[id]][data[[repY]] %in% 1 & data[[ord_var]] %in% A]
  sample.A <- sample(id.A, size = round(length(id.A) * probA), replace = FALSE)
  data[[ord_var]][data[[id]]%in%sample.A] <- NA
  id.B <- data[[id]][data[[repY]] %in% 0 & data[[ord_var]] %in% B]
  sample.B <- sample(id.B, size = round(length(id.B) * probB), replace = FALSE)
  data[[ord_var]][data[[id]]%in%sample.B] <- NA
  return(data)
}
  result <- lapply(unique(data[[cat_var]]), function(i) {
    mnar_by_strata(data = subset(data, subset = data[[cat_var]] == i), 
                  repY = "y", id = "id", ord_var = "x1.mis", 
                  probA = proba[[i]][1], probB = proba[[i]][2], 
                  A = A, B = B, seed = 100)
  })
  return(do.call(rbind, result))
}
```


```{r}
for(s in 1:nsim){
  resultat_s <- simulate_mnar(data = datas[[s]], cat_var = "x2", proba = missing_prob)
  datas[[s]] <- resultat_s[order(resultat_s$id),]
  datas[[s]]$x1.imp <- datas[[s]]$x1.mis
}
```




```{r}
PACKAGE <- function(package){
  for(p in package){
    library(p, character.only = TRUE)
  }
}
```


```{r}
M <- 10
nburn <- 100
nbetween <- 100
cl <- makeCluster(round(detectCores()/2))
formula = y ~ x1.imp + x2 + (1|clus) 
clusterExport(cl, c("datas","M", "nburn", "nbetween", "formula","nsim"))
invisible(clusterCall(cl, PACKAGE, c("jomo","mitml")))

# Calcule en parallele
mar_data <- parLapply(cl, 1:nsim, function(s) {
  jomo.glmer(
    formula = as.formula(formula),
    data = datas[[s]][, c("y", "x1.imp", "x2", "clus")],
    nburn = nburn,
    nbetween = nbetween,
    nimp = M,
    output = 0
  )
})

mar_data <- parLapply(cl, mar_data, function(jomo_imp) {
  jomo2mitml.list(jomo_imp)
})

stopCluster(cl)
```


```{r}
mar_data <- lapply(1:nsim, function(s) {
  lapply(1:M, function(m) {
    cbind(
      mar_data[[s]][[m]],
      datas[[s]][, c("x1", "x1.mis")]
    )
  })
})
```








## Beta Coefficients and Zeta Intercepts of the Ordinal Regression Model with Full Data without Missing Data

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




# Extracting the Means of Beta Coefficients and Zeta Intercepts


```{r}
beta_full <- rowMeans(sapply(results, function(res) res$beta))
zeta_full <- rowMeans(sapply(results, function(res) res$zeta))
sd_full <- rowMeans(sapply(results, function(res) res$sd_beta_zeta))
```




## Ordinal Regression of Data Imputed under MAR

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







### Coefficients from the Ordinal Regression Model where X1 is MAR


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
Beta <- as.data.frame(cbind(beta_full, sd_full = sd_full[c("y1",  "x22", "x23", "x24")],beta_mar,sd_mar = sd_mar[c("y1","x22" ,"x23", "x24")]))
Beta$Bias <- ((Beta$beta_mar - Beta$beta_full)/Beta$beta_full)*100
Beta
```

```{r}
Zeta <- as.data.frame(cbind(zeta_full, sd_full = sd_full[c("1|2", "2|3")], zeta_mar, sd_mar = sd_mar[c("1|2", "2|3")]))
Zeta$Bias <- ((Zeta$zeta_mar - Zeta$zeta_full )/Zeta$zeta_full)*100
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
    theme_clean(base_size = 9) +
    ylim(range(data)[1], range(data)[2])+
    theme(axis.text.x = element_text(angle = 45, hjust = 0.5)) +
    ggtitle(title) +
    xlab(xlab) +      
    ylab(ylab) 
  return(p)
}
```

```{r}
beta_graph <- mar_simule_graph(data = Beta[, c("beta_full", "beta_mar")], title = "Coefficents of Simulated vs Mar", 
                               labels = c( "Simule", "Mar"),colors = couleur[1:2], xlab = "Coefficients", 
                               ylab = "Estimate")
zeta_graph <- mar_simule_graph(data = Zeta[, c("zeta_full", "zeta_mar")], title = "Intercepts of Simulated vs Mar", 
                               labels = c( "Simule", "Mar"),colors = couleur[1:2], xlab = "Intercepts", 
                               ylab = "Estimate")

grid.arrange(grobs = list(beta_graph, zeta_graph), ncol = 2)
```

```{r}
sink()
```

