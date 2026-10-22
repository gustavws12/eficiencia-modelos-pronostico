# Hallazgos de la corrida piloto

**Fecha:** 2026-08-28
**Configuración:** 10 DGP × T = 250 × σ = 1 × 5 réplicas × 15 orígenes × 9 modelos
**Duración:** 15,1 min en 2 núcleos · 0 réplicas fallidas · tasa de éxito 100 % en los 9 modelos

---

## 1. El pipeline funciona de extremo a extremo

Los nueve modelos ajustaron y pronosticaron en las 10 celdas sin una sola falla.
Eso incluye los tres mecanismos de umbral, que era el punto de mayor riesgo del
diseño. Las métricas de densidad (CRPS, PICP, pinball) se calculan para todos los
modelos porque todos entregan trayectorias simuladas, no solo la media.

## 2. Los resultados van en la dirección teórica esperada

**DGP lineales (D1–D5).** ARIMA y el AR(p) bayesiano ocupan los primeros
lugares en las cinco celdas; el ML de árboles y la red neuronal quedan
consistentemente últimos, con RMSE entre 5 % y 19 % peor que ARIMA. Es
exactamente H1 y H2: la flexibilidad se paga cuando el DGP es simple y T = 250.

**Umbral categórico observable (D8).** El TAR con régimen observable es el mejor
modelo de la celda (RMSE 5 % por debajo de ARIMA) y supera al MS-AR, que tiene
que inferir el mismo régimen. Es la evidencia directa de H4 —el valor de
información del indicador de régimen— y es el resultado más limpio del piloto.

**Calibración de intervalos.** La red neuronal sub-cubre de forma sistemática:
PICP al 95 % entre 0,75 y 0,92 en siete de diez celdas, frente a 0,93–0,97 de
los modelos econométricos y bayesianos. Es el germen de H6.

## 3. Una celda hay que recalibrar: D7 (umbral continuo exógeno)

Con σ = 1, el TAR con variable de umbral continua no le gana a ARIMA
(RMSE relativo 1,03), aunque sí le gana al SETAR autoexcitado mal especificado
(1,11). El mecanismo se detecta, pero la señal no alcanza para dominar.

Diagnóstico: con `φ1 = 0,85`, `φ2 = 0,30`, `c1 = −0,30`, `c2 = 0,30` y umbral en
la mediana de `z`, la diferencia entre regímenes en la media condicional es
pequeña frente a σ = 1. Tres respuestas posibles, en orden de preferencia:

1. **Esperar al diseño completo.** La celda σ = 0,5 del diseño factorial debería
   revelar el mecanismo sin tocar los parámetros. Es la opción metodológicamente
   más limpia: la sensibilidad al ruido es justamente lo que el experimento mide.
2. **Aumentar la separación de interceptos** a `c1 = −0,6`, `c2 = 0,6`,
   documentándolo como calibración de diseño previa a la corrida definitiva.
3. **Mover el umbral fuera de la mediana** (p. ej. al percentil 70) para generar
   un régimen minoritario más distintivo, que es lo que se observa en datos
   reales de estrés financiero.

Recomendación: (1) más (3). No tocar `φ` ni `σ`, que son los factores del diseño.

## 4. Dos errores encontrados y corregidos

Quedan documentados porque son del tipo que no falla ruidosamente:

- **`optim` y el emparejamiento parcial de argumentos.** Pasar `p = p` a `optim()`
  hace que R lo empareje con el formal `par`; el MS-AR estaba optimizando un
  vector de un solo elemento y devolviendo "no estimable" sin error. Corregido
  renombrando los argumentos internos (`orden`, `theta`) e invocando `optim` con
  `par =` y `fn =` explícitos.
- **Intervalos de los ensambles de árboles.** Calibrados con residuos in-sample
  daban intervalos al 95 % más angostos que el intervalo al 80 % de ARIMA, lo que
  habría inflado artificialmente el desempeño del ML en cobertura. Ahora la
  dispersión se estima con un holdout temporal del 20 % final del entrenamiento.
- **Alineación del TAR cuando `d > p`.** La variable de umbral `y_{t−d}` y la
  matriz de rezagos arrancaban en índices distintos, desalineando residuos y
  regímenes. Corregido anclando ambos en `max(p, d)`.

## 5. El MS-AR dejó de ser el cuello de botella

El filtro de Hamilton en R puro costaba 6,4 s por ajuste y representaba ~80 % del
tiempo total del experimento. Reescrito en C++ (`R/hamilton.cpp`, cargado con
`Rcpp::sourceCpp`), el ajuste completo —incluida la selección de K por BIC y la
simulación de 500 trayectorias— baja a **0,63 s**, un factor de 10. El proyecto
sigue funcionando sin Rcpp: `00_setup.R` detecta la disponibilidad y usa la
implementación en R puro como respaldo.

Tiempos medianos por ajuste + pronóstico (T = 250, h ≤ 10):

| Modelo | Segundos |
|--------|---------:|
| Naive (RW) | 0,005 |
| AR(p) bayesiano | 0,012 |
| TAR umbral categórico | 0,014 |
| SETAR autoexcitado | 0,038 |
| TAR umbral continuo | 0,043 |
| ARIMA (auto) | 0,095 |
| ML de árboles | 0,143 |
| MS-AR | 0,155 |
| Red neuronal AR | 1,850 |

El nuevo cuello de botella es la red neuronal (nnetar con 20 repeticiones y 200
trayectorias simuladas para los intervalos). Si hace falta recortar, ahí está el
margen: bajar a 100 trayectorias reduce el costo total del experimento a la mitad
con un efecto menor sobre la precisión de las métricas de cobertura.

## 6. Costo de la corrida definitiva

Costo medido: **36,2 segundos-núcleo por réplica** (T = 250, 15 orígenes, 9
modelos). Extrapolando al diseño completo (40 celdas: 10 DGP × 2 T × 2 σ), con el
supuesto de que T = 1000 cuesta 3,5× más que T = 250:

| Réplicas | Horas-núcleo | Con 8 núcleos | Con 32 núcleos |
|---------:|-------------:|--------------:|---------------:|
| 200 | ~181 | ~23 h | ~5,7 h |
| 500 | ~452 | ~57 h | ~14 h |

**Esta es la decisión que hay que tomar.** El entorno donde corrió el piloto tiene
2 núcleos, así que la corrida definitiva no cabe aquí. Tres caminos:

- **200 réplicas en una máquina de 8–16 núcleos.** El error de Monte Carlo del
  RMSE relativo queda cerca de 2 %, suficiente para el mapa de dominancia.
- **500 réplicas en 32 núcleos.** Es lo que pide un referee exigente para las
  celdas de umbral, donde la varianza entre réplicas es mayor.
- **Diseño escalonado:** 500 réplicas en las cuatro celdas de umbral (D6a, D6b,
  D7, D8), que son las que sostienen la contribución del artículo, y 200 en las
  seis celdas lineales, que solo cumplen función de referencia. Costo total
  ~270 horas-núcleo. Es la opción con mejor relación potencia/costo.

## 7. Ajuste pendiente en la prueba de Diebold–Mariano

Con 15 orígenes por réplica, el DM individual casi no rechaza (proporciones de
0 a 0,4 y en su mayoría 0). No es un error: es falta de potencia de una prueba
aplicada a 15 observaciones. Dos correcciones para la corrida definitiva:

1. Reportar la proporción de rechazos **entre réplicas**, que es lo que ya hace
   `agregar_dm()` y que sí tiene potencia con 200–500 réplicas.
2. Añadir un **DM agrupado**, que apile los errores de todas las réplicas de una
   celda antes de contrastar. Es la versión con potencia real y complementa al
   Model Confidence Set.

## 8. Qué sigue

| # | Tarea | Bloqueante |
|---|-------|-----------|
| 1 | Decidir réplicas y máquina para la corrida definitiva | **Sí** |
| 2 | Recalibrar D7 (opción 1 + 3 de la sección 3) | Sí |
| 3 | Añadir el DM agrupado y el Model Confidence Set a `03_evaluacion.R` | No |
| 4 | Descargar las ocho series del anexo (`R/06_datos.R`, requiere internet) | No |
| 5 | Verificar las referencias 2025–2026 del bloque B de `referencias.bib` | No |
| 6 | Evaluar `torch` para LSTM/GRU en la aplicación empírica | No |
