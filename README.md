# Eficiencia y sensibilidad de modelos de pronóstico

Proyecto de investigación: *Modelos de pronóstico temporal en finanzas: una comparación entre
enfoques no lineales clásicos, deep learning y métodos bayesianos.* Estudio de simulación Monte
Carlo sobre diez procesos generadores, evaluación con origen móvil y prueba de Diebold-Mariano,
aplicación empírica y manuscrito en RMarkdown.

```
protocolo/protocolo_simulacion.md   Diseño experimental completo (leer primero)
protocolo/hallazgos_piloto.md       Resultados de la corrida piloto
R/00_setup.R                        Configuración, rutas, semillas, chequeo de paquetes
R/01_dgp.R                          Procesos generadores de datos
R/02_modelos.R                      Modelos de pronóstico con interfaz común
R/03_evaluacion.R                   Métricas, Diebold-Mariano, origen móvil
R/04_monte_carlo.R                  Motor del experimento y agregación
R/05_piloto.R                       Corrida piloto y estimación de costo
R/06_datos.R                        Series de la aplicación empírica
R/hamilton.cpp                      Filtro de Hamilton en C++ (Rcpp)
resultados/                         Salidas del piloto (CSV)
manuscrito/articulo.Rmd · referencias.bib · articulo.pdf
```

## Cómo correr

```bash
PROYECTO_BASE=$(pwd) Rscript R/05_piloto.R 3 6 2     # piloto: réplicas, orígenes, núcleos
```

Autor: Gustavo A Romero
