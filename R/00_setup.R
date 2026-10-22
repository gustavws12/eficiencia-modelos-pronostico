## =========================================================================
## 00_setup.R  —  Configuración global del estudio de simulación
## Proyecto: Eficiencia y sensibilidad de modelos de pronóstico
## =========================================================================

suppressPackageStartupMessages({
  library(stats)
  library(forecast)
})

## --- Rutas ---------------------------------------------------------------
if (!exists("RUTA_BASE")) RUTA_BASE <- getwd()
RUTA_R          <- file.path(RUTA_BASE, "R")
RUTA_RESULTADOS <- file.path(RUTA_BASE, "resultados")
RUTA_FIGURAS    <- file.path(RUTA_BASE, "figuras")
RUTA_DATOS      <- file.path(RUTA_BASE, "datos")
for (d in c(RUTA_RESULTADOS, RUTA_FIGURAS, RUTA_DATOS)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

## --- Opciones globales ---------------------------------------------------
options(stringsAsFactors = FALSE, warn = 1)

## Semilla maestra del experimento
SEMILLA_MAESTRA <- 20260828L

## Horizontes de pronóstico evaluados
HORIZONTES <- c(1L, 5L, 10L)
H_MAX      <- max(HORIZONTES)

## Niveles de los intervalos de predicción
NIVELES <- c(80, 95)

## Número de trayectorias simuladas para densidades predictivas
B_PRED <- 500L

## Burn-in de los DGP
BURN_IN <- 200L

## --- Disponibilidad de paquetes opcionales -------------------------------
## El núcleo del estudio no depende de ellos: TAR/SETAR, MS-AR y el modelo
## bayesiano están implementados en este proyecto (ver 02_modelos.R).
## Estos paquetes se usan, si están disponibles, como verificación cruzada
## o como modelos adicionales.
PKG_OPCIONALES <- c("tsDyn", "MSwM", "bsts", "xgboost", "gbm", "MCS", "torch")
PKG_OK <- setNames(
  vapply(PKG_OPCIONALES, function(p) requireNamespace(p, quietly = TRUE), logical(1)),
  PKG_OPCIONALES
)

## --- Filtro de Hamilton compilado (opcional pero muy recomendable) -------
## Reduce el costo del MS-AR en dos órdenes de magnitud. Si falla, el
## proyecto sigue funcionando con la implementación en R puro.
HAMILTON_CPP <- FALSE
if (requireNamespace("Rcpp", quietly = TRUE)) {
  .cpp <- file.path(RUTA_R, "hamilton.cpp")
  if (file.exists(.cpp)) {
    HAMILTON_CPP <- tryCatch({
      Rcpp::sourceCpp(.cpp, cacheDir = file.path(tempdir(), "rcpp_cache"),
                      verbose = FALSE)
      TRUE
    }, error = function(e) {
      message("Rcpp no pudo compilar hamilton.cpp: ", conditionMessage(e))
      FALSE
    })
  }
}

mensaje_entorno <- function() {
  cat("--- Entorno ---------------------------------------------------\n")
  cat("R:", R.version.string, "\n")
  cat("Filtro de Hamilton compilado (Rcpp):", HAMILTON_CPP, "\n")
  cat("Paquetes opcionales disponibles:\n")
  print(PKG_OK)
  cat("---------------------------------------------------------------\n")
}
