# =====================================================================
# IMPORTATION DES LIBRAIRIES
# =====================================================================

library(ggplot2)  # Pour les visualisations graphiques


# =====================================================================
# FONCTIONS DE BASE POUR LA RÉGRESSION LINÉAIRE
# =====================================================================

# Fonction pour ajouter une constante (intercept) à la matrice X
add_constant <- function(X) {
  cbind(rep(1, nrow(X)), X)
}

# Fonction pour ajuster une régression linéaire
linear_regression <- function(X, y, add_intercept = TRUE) {
  if (!is.matrix(X)) X <- as.matrix(X)         # Convertir en matrice si nécessaire
  if (!is.vector(y)) y <- as.vector(y)         # Convertir y en vecteur si besoin
  
  if (add_intercept) {
    X <- add_constant(X)                       # Ajouter constante si nécessaire
  }
  
  XtX <- t(X) %*% X                             # Produit matriciel XᵀX
  XtX_inv <- solve(XtX)                         # Inverse de XᵀX
  Xty <- t(X) %*% y                             # Produit Xᵀy
  
  beta <- XtX_inv %*% Xty                       # Calcul des coefficients
  return(beta)
}

# Fonction pour prédire les valeurs avec les coefficients
predict_linear <- function(X, beta, add_intercept = TRUE) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (add_intercept) X <- add_constant(X)
  predictions <- X %*% beta
  return(predictions)
}

# Fonction pour calculer le coefficient de détermination R²
r_squared <- function(y_true, y_pred) {
  y_mean <- mean(y_true)
  ss_total <- sum((y_true - y_mean)^2)         # Somme totale des carrés
  ss_residual <- sum((y_true - y_pred)^2)      # Somme des carrés des résidus
  r2 <- 1 - (ss_residual / ss_total)           # Formule du R²
  return(r2)
}


# =====================================================================
# STATISTIQUES DE LA RÉGRESSION
# =====================================================================

# Fonction pour calculer les statistiques : erreur type, t-stats, p-values, etc.
regression_stats <- function(X, y, beta, add_intercept = TRUE) {
  if (!is.matrix(X)) X <- as.matrix(X)
  if (add_intercept) X <- add_constant(X)
  
  n <- nrow(X)
  p <- ncol(X)
  
  y_pred <- X %*% beta
  residuals <- y - y_pred
  SCR <- sum(residuals^2)                      # Somme des carrés des résidus
  mse <- SCR / (n - p)                          # Erreur quadratique moyenne
  
  vcov <- mse * solve(t(X) %*% X)               # Matrice de variance-covariance
  std_errors <- sqrt(diag(vcov))               # Erreurs types
  t_stats <- beta / std_errors                 # Statistiques t
  p_values <- 2 * (1 - pt(abs(t_stats), df = n - p))  # p-values bilatérales
  
  SCT <- sum((y - mean(y))^2)
  r2 <- 1 - (SCR / SCT)                         # R²
  r2_adj <- 1 - (1 - r2) * ((n - 1) / (n - p - 1))  # R² ajusté
  
  f_stat <- ((SCT - SCR) / (p - 1)) / (SCR / (n - p))  # Stat F
  f_p_value <- 1 - pf(f_stat, df1 = p - 1, df2 = n - p)
  
  return(list(
    coefficients = beta,
    std_errors = std_errors,
    t_stats = t_stats,
    p_values = p_values,
    residuals = residuals,
    fitted_values = y_pred,
    r_squared = r2,
    adj_r_squared = r2_adj,
    f_statistic = f_stat,
    f_p_value = f_p_value,
    df = c(p, n - p),
    mse = mse,
    vcov = vcov
  ))
}


# =====================================================================
# DÉTECTION DE MULTICOLINÉARITÉ (VIF)
# =====================================================================

# Fonction pour calculer les VIF à partir de la matrice de corrélation
vif <- function(X) {
  # 1. S'assurer que X est une matrice
  X <- as.matrix(X)
  
  # 2. Supprimer les colonnes constantes (sinon corr = NaN)
  constant_cols <- apply(X, 2, function(col) var(col) == 0)
  if (any(constant_cols)) {
    warning("Colonnes constantes supprimées : ", paste(colnames(X)[constant_cols], collapse = ", "))
    X <- X[, !constant_cols]
  }

  # 3. Calculer la matrice de corrélation
  corr_matrix <- cor(X)

  # 4. Inverser la matrice de corrélation
  inv_corr <- solve(corr_matrix)

  # 5. Le VIF est la diagonale de l’inverse de la matrice de corrélation
  vif_values <- diag(inv_corr)

  # 6. Nommer les VIF avec les noms des colonnes
  names(vif_values) <- colnames(X)
  
  return(vif_values)
}


# =====================================================================
# IMPORTATION DES DONNÉES ET UTILISATION DES FONCTIONS
# =====================================================================

# Charger les données
data <- read.csv("/Users/salim/Desktop/Insea/S2/PARTIE 1/Regression Lineaire/R/DataSetMacroeconomy.csv", sep = ",", header = TRUE)
summary(data)

# Séparer les variables explicatives (X) et la variable cible (y)
X <- data[, -ncol(data)]
y <- data[, ncol(data)]

# Estimation des coefficients
beta <- linear_regression(X, y)

# Prédiction
y_pred <- predict_linear(X, beta)

# Affichage du R²
r2 <- r_squared(y, y_pred)
print(paste("Coefficient de détermination R² :", r2))

# Affichage des statistiques de la régression
stats <- regression_stats(X, y, beta)
print(stats)

# Affichage des valeurs de VIF
vif_values <- vif(X)
print("Valeurs VIF :")
print(vif_values)


# =====================================================================
# DIAGNOSTICS VISUELS
# =====================================================================

# Graphique des résidus vs prédictions (homoscédasticité)
y_pred <- stats$fitted_values
residuals = stats$residuals
ggplot(data, aes(x = y_pred, y = stats$residuals)) +
  geom_point() +
  geom_hline(yintercept = 0, color = "red") +
  labs(title = "Résidus vs Valeurs Prédites",
       x = "Valeurs Prédites",
       y = "Résidus") +
  theme_minimal()

# QQ plot (normalité des résidus)
ggplot(data, aes(sample = stats$residuals)) +
  stat_qq() +
  stat_qq_line(color = "red") +
  labs(title = "QQ Plot des Résidus",
       x = "Quantiles Théoriques",
       y = "Quantiles Observés") +
  theme_minimal()
