# Protocolo de diseño Monte Carlo
### Sensibilidad de modelos de pronóstico: no lineales clásicos, deep learning y bayesianos

**Proyecto:** *Modelos de pronóstico temporal en finanzas: una comparación entre enfoques no lineales clásicos, deep learning y métodos bayesianos*
**Versión:** 0.1 — 2026-08-28
**Entorno:** R (>= 4.4), manuscrito en RMarkdown

---

## 1. Objetivo del experimento

El estudio de simulación responde a la primera pregunta de investigación del artículo:

> ¿En qué contextos (tipo de proceso generador, tamaño muestral, nivel de ruido, presencia y forma de los regímenes) cada familia de modelos supera a las demás en precisión predictiva?

La estrategia es un experimento factorial controlado en el que **se conoce el proceso generador de datos (DGP)**. Esto permite separar tres efectos que en datos reales quedan confundidos:

1. **Sesgo de especificación**: qué pierde un modelo flexible (LSTM, XGBoost) cuando el DGP es lineal y simple.
2. **Ganancia por no linealidad**: qué gana un modelo de umbral o de cambio de régimen cuando el DGP efectivamente tiene regímenes, y cuánta señal se necesita para detectarlos.
3. **Costo de la incertidumbre paramétrica**: cómo se degrada cada familia cuando T es pequeño o la razón señal-ruido es baja.

El resultado esperado es un **mapa de dominancia condicional** (qué familia gana, bajo qué condiciones), que luego se contrasta con la aplicación empírica del Anexo.

---

## 2. Procesos generadores de datos (DGP)

Notación: `e_t ~ iid N(0, sigma^2)`; `T` es la longitud de la serie simulada; se descartan 200 observaciones de calentamiento (*burn-in*) en todos los casos.

### 2.1 Bloque lineal (referencia)

| ID | Proceso | Especificación | Parámetros base |
|----|---------|----------------|-----------------|
| **D1** | AR(1) | `y_t = c + phi*y_{t-1} + e_t` | `c = 0`, `phi in {0.5, 0.9}` |
| **D2** | ARMA(1,1) | `y_t = c + phi*y_{t-1} + e_t + theta*e_{t-1}` | `phi = 0.6`, `theta = 0.4` |
| **D3** | ARIMA(1,1,1) | `Dy_t = mu + phi*Dy_{t-1} + e_t + theta*e_{t-1}` | `mu = 0.02`, `phi = 0.4`, `theta = 0.3` |
| **D4** | ARIMA(1,1,0) | `Dy_t = mu + phi*Dy_{t-1} + e_t` | `mu = 0.02`, `phi = 0.4` |
| **D5** | ARIMA(0,0,1) = MA(1) | `y_t = c + e_t + theta*e_{t-1}` | `c = 0`, `theta = 0.6` |

Función de estos cinco DGP: son los casos en que **el benchmark lineal es el modelo correcto**. Miden el *costo de la flexibilidad*. La hipótesis de trabajo es que ARIMA domina, que BSTS iguala a ARIMA con mejor calibración de intervalos, y que el deep learning pierde de forma sistemática en T pequeño y converge hacia ARIMA cuando T crece.

`phi = 0.9` en D1 introduce un segundo régimen de persistencia (raíz cercana a la unidad) que replica el comportamiento de series macro como la inflación o la tasa de política.

### 2.2 Bloque de umbrales

Aquí está el núcleo del pedido: **modelos de umbral con variable de referencia continua y con variable de referencia categórica**. Se distinguen tres mecanismos.

**D6 — SETAR(2;1,1): umbral autoexcitado (variable de umbral continua endógena)**

```
y_t = c1 + phi1*y_{t-d} ... :  y_t = c1 + phi1*y_{t-1} + e_t   si  y_{t-d} <= r
                               y_t = c2 + phi2*y_{t-1} + e_t   si  y_{t-d} >  r
```

Parámetros base: `d = 1`, `r = 0`, `c1 = -0.3`, `phi1 = 0.85`, `c2 = 0.3`, `phi2 = 0.30`.
Se verifica estacionariedad ergódica por la condición suficiente de Petruccelli–Woolford (`phi1 < 1`, `phi2 < 1`, `phi1*phi2 < 1`).

**Factor de intensidad de no linealidad:** `delta = |phi1 - phi2| in {0.20, 0.55}`. Es la palanca central del experimento: mide cuánta no linealidad hace falta para que TAR/MS/DL le ganen a ARIMA.

**D7 — TAR con variable de umbral continua exógena**

```
y_t = c1 + phi1*y_{t-1} + e_t   si  z_{t-1} <= r
y_t = c2 + phi2*y_{t-1} + e_t   si  z_{t-1} >  r
```

con `z_t = 0.7*z_{t-1} + u_t`, `u_t ~ N(0,1)`, independiente de `e_t`, y `r = quantil 0.5 de z`.

Este caso separa dos habilidades distintas: **conocer la variable de transición** (ventaja del TAR bien especificado) frente a **inferirla de los datos** (lo que hacen el MS-AR, el LSTM y el XGBoost al recibir `z` como covariable). Se corren dos versiones de cada modelo: con `z` observable y sin ella.

**D8 — TAR con variable de umbral categórica exógena (regímenes observables)**

```
s_t in {1,2,3}   (cadena de Markov de primer orden, matriz P)
y_t = c_{s_t} + phi_{s_t}*y_{t-1} + sigma_{s_t}*e_t
```

Parámetros base (interpretación macro: expansión / neutral / estrés):

| Régimen | `c` | `phi` | `sigma` |
|---------|-----|-------|---------|
| 1 (expansión) |  0.40 | 0.30 | 0.8 |
| 2 (neutral)   |  0.00 | 0.70 | 1.0 |
| 3 (estrés)    | -0.60 | 0.85 | 2.2 |

Matriz de transición base (regímenes persistentes):

```
P = [0.95 0.04 0.01
     0.03 0.94 0.03
     0.02 0.08 0.90]
```

Aquí `s_t` es **observable** para el modelo TAR categórico y **latente** para el MS-AR. La comparación D8 mide directamente el valor de información del régimen observado, que es la pregunta práctica relevante en gestión de riesgo (¿sirve condicionar a un indicador de régimen publicado, o conviene estimarlo?).

**Variante D8b (robustez):** `s_t` observable con error de clasificación del 10 %, para evaluar la sensibilidad del TAR categórico a la calidad del indicador de régimen.

---

## 3. Factores experimentales

| Factor | Niveles | Justificación |
|--------|---------|---------------|
| DGP | D1(phi=.5), D1(phi=.9), D2, D3, D4, D5, D6(delta baja), D6(delta alta), D7, D8 → **10 celdas** | Cubre lineal / integrado / umbral continuo / umbral categórico |
| Tamaño muestral `T` | 250, 1000 | Frecuencia baja (macro trimestral/mensual) vs. alta (financiera diaria) |
| Razón señal-ruido | `sigma in {0.5, 2.0}` | Aproxima el contraste macro (poco ruido) vs. financiero (mucho ruido) |
| Horizonte `h` | 1, 5, 10 | Táctico vs. estratégico |
| Réplicas `R` | 500 (definitivo) / 20 (piloto) | Error de Monte Carlo del RMSE relativo < 1.5 % |

Total de celdas de diseño: `10 x 2 x 2 = 40`; series simuladas en la corrida definitiva: `40 x 500 = 20.000`.

Semilla reproducible: `set.seed(20260828 + id_replica)`, con flujo L'Ecuyer-CMRG para paralelización.

---

## 4. Modelos evaluados

### 4.1 Núcleo del Monte Carlo (un representante por familia + benchmarks)

| ID | Modelo | Familia | Implementación en R |
|----|--------|---------|---------------------|
| **M0** | Random walk / naive | Benchmark trivial | `forecast::naive` |
| **M1** | ARIMA automático | Benchmark lineal | `forecast::auto.arima` (BIC) |
| **M2** | SETAR / TAR | No lineal clásico | `tsDyn::setar`, umbral por rejilla + BIC |
| **M3** | MS-AR (2 y 3 regímenes) | No lineal clásico | `MSwM::msmFit` (EM) |
| **M4** | BSTS (nivel local + AR) | Bayesiano | `bsts::bsts`, 2.000 draws MCMC |
| **M5** | XGBoost sobre rezagos | ML clásico | `xgboost`, `p = 12` rezagos |
| **M6** | Red neuronal autorregresiva | Deep learning | `forecast::nnetar` (NNAR) como base; `torch` LSTM si el entorno lo permite |

> **Nota de implementación (M6):** la corrida definitiva usa LSTM en `torch` (R). El piloto valida primero con `nnetar`, que es liviano y sirve de cota inferior del desempeño de una red. Si `torch` no queda operativo en el entorno de cómputo, la LSTM se reserva para la aplicación empírica y el Monte Carlo se reporta con NNAR, declarándolo explícitamente en la sección de metodología.

### 4.2 Zoológico completo (solo en la aplicación empírica del Anexo)

GRU, CNN-LSTM, CNN-LSTM-Attention, Transformer tipo EMAT, LightGBM, Random Forest, BVAR, MS-GARCH bayesiano.

---

## 5. Diseño de evaluación

1. **Partición**: para cada serie simulada de longitud `T`, los primeros 70 % son entrenamiento, el 15 % siguiente validación (selección de hiperparámetros), el 15 % final prueba estrictamente fuera de muestra.
2. **Origen móvil**: ventana expansiva (*expanding window*) con re-estimación cada `k` pasos (`k = 10` en T=250, `k = 25` en T=1000) para acotar el costo computacional; se documenta el efecto de `k` en un chequeo de robustez.
3. **Pronósticos**: directos para `h > 1` en los modelos de ML/DL; iterados en los modelos econométricos, que es su forma natural.

### 5.1 Métricas de precisión puntual

- **RMSE** y **MAE** (escala del DGP).
- **MASE** en lugar de MAPE. *Razón:* varios DGP cruzan el cero y el MAPE se vuelve inestable o infinito; MASE es escala-libre y comparable entre celdas. Se reporta MAPE solo en la aplicación empírica sobre series estrictamente positivas (IPC, TRM, COLCAP, Brent).
- **U de Theil** relativo al naive.

### 5.2 Métricas de incertidumbre

- **PICP** (cobertura empírica del intervalo al 80 % y 95 %) y **MPIW** (ancho medio normalizado).
- **CRPS** para los modelos con densidad predictiva completa (BSTS, MS-AR, ARIMA); es la métrica que permite argumentar la ventaja bayesiana más allá del error puntual.
- **Pinball loss** en los cuantiles 5 % y 95 % (puente hacia el backtesting de VaR de la aplicación).

### 5.3 Pruebas de superioridad predictiva

- **Diebold–Mariano** (con corrección de muestra pequeña de Harvey–Leybourne–Newbold) por réplica, contra ARIMA como pivote; se reporta la **proporción de rechazos** al 5 % en cada celda.
- **Giacomini–White** para comparaciones con parámetros estimados y ventanas móviles (el caso correcto aquí, dado que se re-estima).
- **Model Confidence Set** (Hansen, Lunde y Nason) agregando las réplicas de cada celda: entrega el conjunto de modelos estadísticamente indistinguibles del mejor, que es la forma honesta de reportar el "ganador".

### 5.4 Costo computacional

Se registra tiempo de estimación y de pronóstico por modelo y celda (mediana y P90). El artículo reporta el **frente de Pareto precisión–costo**, que es una de las contribuciones prometidas en la introducción.

---

## 6. Hipótesis contrastables

| # | Hipótesis | Celdas que la prueban |
|---|-----------|-----------------------|
| H1 | Con DGP lineal, ARIMA y BSTS son indistinguibles en RMSE, pero BSTS domina en CRPS y cobertura | D1–D5 |
| H2 | El deep learning solo alcanza a ARIMA con DGP lineal cuando `T >= 1000` y ruido bajo | D1–D5 × T × sigma |
| H3 | La ventaja de TAR/MS-AR aparece a partir de un umbral de intensidad `delta` y se desvanece con ruido alto | D6 (delta baja vs. alta) × sigma |
| H4 | Conocer la variable de transición (continua u observable) domina a inferirla; la brecha se cierra cuando T crece | D7, D8 (versiones con y sin `z`/`s`) |
| H5 | Con regímenes categóricos observables y ruido heterocedástico por régimen, la familia bayesiana domina en cuantificación de incertidumbre aunque no en RMSE | D8 |
| H6 | El ML de árboles (XGBoost) es el más robusto a especificación errónea pero el peor calibrado en intervalos | Todas |

---

## 7. Anexo empírico: series de aplicación

| Serie | Frecuencia | Fuente | Rol en el diseño |
|-------|-----------|--------|------------------|
| Inflación Colombia (IPC var. anual) | Mensual | DANE / BanRep | Macro persistente, con quiebre 2021–2023 |
| PIB Colombia (var. real) | Trimestral | DANE | Muestra corta, régimen recesivo observable |
| TRM COP/USD | Diaria | BanRep | Financiera, alta volatilidad, ligada a Brent |
| USD/EUR | Diaria | BCE / FRED | Financiera desarrollada, benchmark internacional |
| Brent / WTI | Diaria | EIA / FRED | Commodity con quiebres estructurales fuertes |
| COLCAP | Diaria | BVC | Índice bursátil emergente, regímenes de volatilidad |
| Tasa de política / IBR | Diaria–mensual | BanRep | Muy persistente; **variable de umbral categórica natural** (ciclo alcista / estable / bajista) |
| EMBI Colombia | Diaria | JP Morgan / BanRep | **Variable de umbral continua natural**; insumo para el backtesting de VaR |

**Diseño del anexo.** Las series no se tratan como ocho ejercicios independientes, sino como tres bloques que replican los factores del Monte Carlo:

- **Bloque macro (baja frecuencia, poco ruido):** inflación, PIB, tasa de política. Contrasta H1 y H2 en la región donde el deep learning debería perder.
- **Bloque financiero (alta frecuencia, mucho ruido):** TRM, USD/EUR, Brent, COLCAP. Contrasta H2 y H6.
- **Bloque de régimen:** EMBI como variable de umbral continua y la fase del ciclo de política monetaria como variable de umbral categórica, aplicadas sobre TRM y COLCAP. Contrasta H3, H4 y H5 con datos reales, y conecta con el backtesting de VaR (Kupiec, Christoffersen) y el Sharpe de una estrategia condicionada al régimen.

Horizonte objetivo: 2005–2025 (>= 20 años) para las series diarias; desde 2000 para inflación y PIB.

---

## 8. Reproducibilidad

- Semillas fijas y registro de `sessionInfo()` en cada corrida.
- Resultados crudos por réplica en `resultados/mc_raw_*.parquet`; agregados en `resultados/mc_agg_*.csv`.
- El manuscrito (`manuscrito/articulo.Rmd`) lee únicamente los agregados, de modo que compila sin re-ejecutar la simulación.
- Control de versiones de paquetes vía `renv` (pendiente de activar al fijar la corrida definitiva).

---

## 9. Decisiones pendientes

1. Confirmar `torch`/LSTM en el entorno de cómputo o dejar NNAR como representante de DL en el Monte Carlo.
2. Fijar `R = 500` o reducir a `R = 300` según los tiempos medidos en el piloto.
3. Definir la fuente exacta y el método de descarga de EMBI Colombia (no siempre disponible en abierto).
4. Revista objetivo: define el límite de palabras y el estilo de referencias del `.Rmd`.
