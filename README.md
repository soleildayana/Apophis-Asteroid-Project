# Apophis Asteroid Project

Proyecto de mecánica celeste para estudiar el encuentro cercano del asteroide **(99942) Apophis** con la Tierra en 2029 mediante simulación numérica y análisis orbital.

## Estructura actual del repositorio

```text
Apophis-Asteroid-Project/
├── README.md
├── agents.md
├── Orbit_Viewer_3D.ipynb
├── analisis/
│   ├── README.md
│   ├── nb01_TransfHeliocentricaObservaciones.ipynb
│   ├── nb02_Cuadraturas_2BRP.ipynb
│   ├── nb03_Anomalias.ipynb
│   ├── nb4_kepler_benchmark.ipynb
│   ├── nb5_hodografo.ipynb
│   ├── nb6_crtbp_Lpuntos_Cjacobi.ipynb
│   ├── nb7_perturbaciones_gauss.ipynb
│   └── observaciones_astrometricas.csv
└── modelos/
    ├── README.md
    ├── modelo1_SEMAJ.ipynb
    ├── modelo2_SEMAJV.ipynb
    └── modelo3_completo.ipynb
```

## Qué hay en cada carpeta

- **`modelos/`**: pipeline principal de integración N-cuerpos con complejidad creciente (SEMAJ → SEMAJV → modelo final de 9 cuerpos).  
  Ver: [`modelos/README.md`](modelos/README.md)
- **`analisis/`**: notebooks de análisis y fundamentos orbitales (determinación de órbita, cuadraturas, anomalías, Kepler, hodógrafo, CRTBP/Jacobi y Gauss).  
  Ver: [`analisis/README.md`](analisis/README.md)
- **`Orbit_Viewer_3D.ipynb`**: visualización 3D de órbitas con rotaciones de Euler y elementos orbitales.

## Flujo recomendado de lectura

1. Revisar **`modelos/`** para entender la construcción del resultado dinámico principal.
2. Profundizar en **`analisis/`** para validar e interpretar la física orbital.
3. Usar **`Orbit_Viewer_3D.ipynb`** para visualización cualitativa de geometría orbital.

## Dependencia principal

- [`pymcel`](https://github.com/seap-udea/pymcel): utilidades de mecánica celeste, efemérides e integración usadas en los notebooks.
