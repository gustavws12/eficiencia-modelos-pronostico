## =========================================================================
## 05_piloto.R  —  Corrida piloto
##
## Propósito: (i) validar que todo el pipeline corre de extremo a extremo,
## (ii) medir tiempos por modelo para dimensionar la corrida definitiva,
## (iii) verificar que los resultados van en la dirección teórica esperada.
##
## Uso:  Rscript R/05_piloto.R  [n_rep]
## =========================================================================

RUTA_BASE <- if (nzchar(Sys.getenv("PROYECTO_BASE"))) Sys.getenv("PROYECTO_BASE") else getwd()
source(file.path(RUTA_BASE, "R", "00_setup.R"))
source(file.path(RUTA_BASE, "R", "01_dgp.R"))
source(file.path(RUTA_BASE, "R", "02_modelos.R"))
source(file.path(RUTA_BASE, "R", "03_evaluacion.R"))
source(file.path(RUTA_BASE, "R", "04_monte_carlo.R"))

mensaje_entorno()

args <- commandArgs(trailingOnly = TRUE)
N_REP      <- if (length(args) >= 1) as.integer(args[1]) else 5L
N_ORIGENES <- if (length(args) >= 2) as.integer(args[2]) else 8L
N_CORES    <- if (length(args) >= 3) as.integer(args[3]) else
                max(1L, min(parallel::detectCores(), 8L))

cat(sprintf("\nPILOTO: %d réplicas por celda, %d orígenes, %d núcleos\n\n",
            N_REP, N_ORIGENES, N_CORES))

## Diseño reducido del piloto: todos los DGP, un solo T y un solo sigma
diseno_piloto <- construir_diseno(
  dgps   = names(DGP_GRID),
  Ts     = 250L,
  sigmas = 1.0,
  n_rep  = N_REP
)

t0 <- Sys.time()
mc <- correr_montecarlo(diseno_piloto,
                        n_origenes = N_ORIGENES,
                        n_cores = N_CORES,
                        etiqueta = "piloto",
                        verbose = TRUE)
cat(sprintf("\nDuración total: %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))

ag <- agregar_montecarlo(mc$resumen)
write.csv(ag, file.path(RUTA_RESULTADOS, "piloto_agregado.csv"), row.names = FALSE)

agdm <- agregar_dm(mc$dm)
if (!is.null(agdm))
  write.csv(agdm, file.path(RUTA_RESULTADOS, "piloto_dm_agregado.csv"), row.names = FALSE)

## --- Reportes rápidos en consola ----------------------------------------
cat("\n=== RMSE relativo a ARIMA (h = 1) ===\n")
sub <- ag[ag$h == 1, c("dgp", "etiqueta", "rmse_med", "rmse_rel_arima",
                       "crps_med", "picp95", "tiempo_med")]
for (d in unique(sub$dgp)) {
  cat("\n--", DGP_GRID[[d]]$etiqueta, "--\n")
  s <- sub[sub$dgp == d, ]
  s <- s[order(s$rmse_rel_arima), ]
  print(format(s[, -1], digits = 3), row.names = FALSE)
}

cat("\n=== Tiempo mediano de ajuste+pronóstico por modelo (segundos) ===\n")
tt <- tapply(ag$tiempo_med, ag$etiqueta, median, na.rm = TRUE)
print(round(sort(tt), 3))

cat("\n=== Tasa de éxito por modelo ===\n")
print(round(tapply(ag$tasa_ok, ag$etiqueta, mean, na.rm = TRUE), 3))

if (is.null(agdm))
  cat("\nAVISO: no se calculó Diebold-Mariano. La prueba exige al menos 8\n",
      "orígenes comunes; esta corrida usó ", N_ORIGENES,
      ". La corrida definitiva debe usar >= 15.\n", sep = "")

## --- Extrapolación de costo de la corrida definitiva ---------------------
## Basada en el tiempo de reloj efectivamente medido, no en una suma de
## medianas (que ignora la sobrecarga del bucle de orígenes).
dur_min      <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
n_tareas     <- nrow(diseno_piloto) * N_REP
seg_nucleo   <- dur_min * 60 * N_CORES / n_tareas      # costo por réplica
n_celdas_def <- length(DGP_GRID) * 2 * 2               # DGP x T x sigma
N_ORIG_DEF   <- 15L
FACTOR_T1000 <- 3.5   # supuesto: costo ~ lineal-superlineal en T; verificar

for (n_rep_def in c(200L, 500L)) {
  costo <- seg_nucleo * (N_ORIG_DEF / N_ORIGENES) *
    (n_celdas_def / 2) * (1 + FACTOR_T1000) / 2 * n_rep_def * 2 / 3600
  cat(sprintf("\nCorrida definitiva con %d réplicas, %d orígenes, %d celdas:\n",
              n_rep_def, N_ORIG_DEF, n_celdas_def))
  cat(sprintf("  ~%.0f horas-núcleo | ~%.0f h con %d núcleos | ~%.1f h con 32 núcleos\n",
              costo, costo / N_CORES, N_CORES, costo / 32))
}
cat(sprintf("\n(Costo medido por réplica en este piloto: %.1f s-núcleo con %d orígenes)\n",
            seg_nucleo, N_ORIGENES))

cat("\nResultados en:", RUTA_RESULTADOS, "\n")
