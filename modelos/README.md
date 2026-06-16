# README — carpeta `modelos`

Esta carpeta contiene el pipeline principal de simulación N-cuerpos del proyecto.

## Objetivo

Estimar y refinar el encuentro Tierra–Apophis incrementando progresivamente la complejidad del sistema dinámico.

## Modelos incluidos

- **`modelo1_SEMAJ.ipynb`**  
  Modelo SEMAJ (5 cuerpos): Sol, Tierra, Luna, Apophis y Júpiter.

- **`modelo2_SEMAJV.ipynb`**  
  Extensión SEMAJV (6 cuerpos): añade Venus al modelo anterior.

- **`modelo3_completo.ipynb`**  
  Modelo final (9 cuerpos): incluye perturbadores adicionales y consolida el análisis principal del proyecto.

## Rol dentro del proyecto

`modelos/` produce las trayectorias y métricas principales (como la distancia Tierra–Apophis), mientras que `analisis/` profundiza la interpretación física y matemática de esos resultados.
