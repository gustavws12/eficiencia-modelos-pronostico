## =========================================================================
## 06_datos.R  —  Descarga de las series del anexo empírico
##
## IMPORTANTE: este script requiere acceso a internet. El entorno de cómputo
## del estudio de simulación no lo tiene, así que se ejecuta por separado y
## los CSV resultantes quedan en datos/.
##
##   Rscript R/06_datos.R
##
## Cada función devuelve un data.frame con columnas: fecha (Date), valor.
## Ninguna serie se transforma aquí: la limpieza y las transformaciones
## (log-retornos, diferencias, alineación de frecuencias) van en 07_aplicacion.R,
## para que el dato crudo quede auditable.
## =========================================================================

RUTA_BASE  <- if (nzchar(Sys.getenv("PROYECTO_BASE"))) Sys.getenv("PROYECTO_BASE") else getwd()
RUTA_DATOS <- file.path(RUTA_BASE, "datos")
if (!dir.exists(RUTA_DATOS)) dir.create(RUTA_DATOS, recursive = TRUE)

FECHA_INICIO <- as.Date("2005-01-01")

guardar <- function(df, nombre) {
  if (is.null(df) || !nrow(df)) { message("SIN DATOS: ", nombre); return(invisible(NULL)) }
  df <- df[order(df$fecha), ]
  df <- df[df$fecha >= FECHA_INICIO & is.finite(df$valor), ]
  f <- file.path(RUTA_DATOS, paste0(nombre, ".csv"))
  utils::write.csv(df, f, row.names = FALSE)
  message(sprintf("%-14s %5d obs  %s a %s", nombre, nrow(df),
                  min(df$fecha), max(df$fecha)))
  invisible(df)
}

## -------------------------------------------------------------------------
## FRED (Reserva Federal de St. Louis) — sin llave de API
## -------------------------------------------------------------------------
fred <- function(id) {
  url <- sprintf("https://fred.stlouisfed.org/graph/fredgraph.csv?id=%s", id)
  d <- tryCatch(utils::read.csv(url, stringsAsFactors = FALSE),
                error = function(e) { message("FRED ", id, ": ", conditionMessage(e)); NULL })
  if (is.null(d)) return(NULL)
  names(d)[1:2] <- c("fecha", "valor")
  d$fecha <- as.Date(d$fecha)
  d$valor <- suppressWarnings(as.numeric(as.character(d$valor)))
  d[, c("fecha", "valor")]
}

## -------------------------------------------------------------------------
## datos.gov.co (Socrata) — series del Banco de la República
## -------------------------------------------------------------------------
socrata <- function(dataset, col_fecha, col_valor, limite = 200000) {
  url <- sprintf("https://www.datos.gov.co/resource/%s.csv?$limit=%d&$order=%s",
                 dataset, limite, col_fecha)
  d <- tryCatch(utils::read.csv(url, stringsAsFactors = FALSE),
                error = function(e) { message("Socrata ", dataset, ": ", conditionMessage(e)); NULL })
  if (is.null(d) || !all(c(col_fecha, col_valor) %in% names(d))) {
    if (!is.null(d)) message("Columnas disponibles en ", dataset, ": ",
                             paste(names(d), collapse = ", "))
    return(NULL)
  }
  data.frame(fecha = as.Date(substr(d[[col_fecha]], 1, 10)),
             valor = suppressWarnings(as.numeric(d[[col_valor]])))
}

## -------------------------------------------------------------------------
## Series
## -------------------------------------------------------------------------

## 1. TRM COP/USD (diaria) — dataset histórico de la Superintendencia Financiera
descargar_trm <- function() {
  d <- socrata("32sa-8pi3", "vigenciadesde", "valor")
  guardar(d, "trm_cop_usd")
}

## 2. USD/EUR (diaria) — FRED DEXUSEU (dólares por euro)
descargar_usdeur <- function() guardar(fred("DEXUSEU"), "usd_eur")

## 3. Brent (diaria) — FRED DCOILBRENTEU;  WTI — DCOILWTICO
descargar_brent <- function() guardar(fred("DCOILBRENTEU"), "brent")
descargar_wti   <- function() guardar(fred("DCOILWTICO"),   "wti")

## 4. Inflación Colombia (mensual, variación anual del IPC) — FRED
##    COLCPIALLMINMEI es el índice; la variación anual se calcula en 07.
descargar_ipc <- function() guardar(fred("COLCPIALLMINMEI"), "ipc_colombia")

## 5. PIB Colombia (trimestral, volumen encadenado) — FRED
descargar_pib <- function() guardar(fred("NGDPRSAXDCCOQ"), "pib_colombia")

## 6. Tasa de política monetaria Colombia — FRED (tasa de interés de política)
descargar_tasa_politica <- function() guardar(fred("IRSTCI01COM156N"), "tasa_politica")

## -------------------------------------------------------------------------
## Series que requieren descarga manual (no hay API abierta estable)
## -------------------------------------------------------------------------
## 7. COLCAP  — Banco de la República, "Índice de la Bolsa de Valores"
##    https://www.banrep.gov.co/es/estadisticas/indice-bursatil-colombia
##    Descargar el Excel, exportar a datos/colcap.csv con columnas fecha,valor.
##
## 8. IBR overnight — Banco de la República
##    https://www.banrep.gov.co/es/estadisticas/indicador-bancario-referencia-ibr
##    Guardar como datos/ibr.csv
##
## 9. EMBI Colombia — J.P. Morgan, difundido por el Banco de la República
##    https://www.banrep.gov.co/es/estadisticas/indicadores-riesgo-pais
##    Guardar como datos/embi_colombia.csv
##
## Si alguna de estas fuentes cambia de formato, la alternativa es
## investing.com (COLCAP) o Bloomberg/Refinitiv si hay licencia.

verificar_manuales <- function() {
  faltan <- c()
  for (f in c("colcap", "ibr", "embi_colombia")) {
    if (!file.exists(file.path(RUTA_DATOS, paste0(f, ".csv")))) faltan <- c(faltan, f)
  }
  if (length(faltan))
    message("\nPendientes de descarga manual: ", paste(faltan, collapse = ", "),
            "\nVer las URL en los comentarios de este script.")
}

## -------------------------------------------------------------------------
if (sys.nframe() == 0L) {
  message("Descargando series del anexo empírico...\n")
  descargar_trm()
  descargar_usdeur()
  descargar_brent()
  descargar_wti()
  descargar_ipc()
  descargar_pib()
  descargar_tasa_politica()
  verificar_manuales()
  message("\nListo. Archivos en: ", RUTA_DATOS)
}
