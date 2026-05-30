# Apophis Asteroid Project

<div align="center">
  <img
    src="https://media.giphy.com/media/v1.Y2lkPWVjZjA1ZTQ3YjhzdTg4ZXducGZpZXU0aWYzc3ZtZ2lsdDM4ZTNydTI1dXg3ZWlieSZlcD12MV9naWZzX3JlbGF0ZWQmY3Q9Zw/CFONrl19EZ0s0/giphy.gif"
    alt="Apophis GIF"
    width="520"
  />
</div>

Proyecto de **Mecánica Celeste** que estudia la aproximación a la Tierra del asteroide **(99942) Apophis** en **abril de 2029**, mediante integración numérica del **problema de N‑cuerpos** (gravedad newtoniana), con:
- comparación progresiva de modelos (3 → 4 → 5 → 6 → 13 cuerpos) para evaluar **convergencia** de la distancia mínima, y
- un modelo final de **13 cuerpos** con **refinamiento temporal** alrededor de la ventana del encuentro.

---

## Estructura del repositorio

```
Apophis-Asteroid-Project/
├── main.ipynb                        ← Reporte principal consolidado (académico/científico)
├── modelos/                         ← Integración N-cuerpos (pipeline principal)
│   ├── modelo1_SEMAJ.ipynb
│   ├── modelo2_SEMAJV.ipynb
│   └── modelo3_completo.ipynb
├── analisis/                        ← Análisis avanzados (mecánica celeste profunda)
│   ├── nb4_kepler_benchmark.ipynb
│   ├── nb5_hodografo.ipynb
│   ├── nb6_crtbp_jacobi_tisserand.ipynb
│   └── nb7_perturbaciones_gauss.ipynb
├── Orbit_Viewer.ipynb
├── cuadraturas_PRelativo2cuerpos.ipynb
├── README.md
└── agents.md
```

### Flujo narrativo recomendado

```
main.ipynb (reporte final consolidado)
       ↑
modelos + analisis (notebooks de soporte)
```

---

## Modelos N-cuerpos (`modelos/`)

Esta carpeta contiene la evolución progresiva del modelo dinámico, donde el **Modelo 3** es el resultado final:

- **Modelo 1 — SEMAJ (5 cuerpos):** Sol–Tierra–Luna–Apophis–Júpiter  
  Notebook: [`modelos/modelo1_SEMAJ.ipynb`](modelos/modelo1_SEMAJ.ipynb)

- **Modelo 2 — SEMAJV (6 cuerpos):** Sol–Tierra–Luna–Apophis–Júpiter–Venus  
  Notebook: [`modelos/modelo2_SEMAJV.ipynb`](modelos/modelo2_SEMAJV.ipynb)

- **Modelo 3 (final) — 9 cuerpos + convergencia + incertidumbre:**  
  Sol–Mercurio–Venus–Tierra–Luna–Marte–Júpiter–Saturno–Apophis  
  Notebook: [`modelos/modelo3_completo.ipynb`](modelos/modelo3_completo.ipynb)

---

## Análisis avanzados (`analisis/`)

Cuatro notebooks independientes que exploran la mecánica celeste de Apophis con mayor profundidad teórica:

### NB4 — Solver de Kepler: Benchmarking
**Archivo:** [`analisis/nb4_kepler_benchmark.ipynb`](analisis/nb4_kepler_benchmark.ipynb)

Implementa una clase `KeplerSolver` que compara tres métodos iterativos para resolver la ecuación trascendental de Kepler $M = E - e\sin E$:
- Newton-Raphson (convergencia cuadrática)
- Punto Fijo (convergencia lineal)
- Laguerre-Conway (convergencia de orden superior, robusto)

Evalúa el número de iteraciones y tiempo CPU a lo largo de la órbita de Apophis, con énfasis en la región del perigeo.

### NB5 — El Hodógrafo: Geometría de la Velocidad
**Archivo:** [`analisis/nb5_hodografo.ipynb`](analisis/nb5_hodografo.ipynb)

Visualiza el hodógrafo de Apophis: el vector velocidad describe un **círculo** en el espacio $(v_x, v_y)$ (Teorema de Hamilton). Muestra cómo el hodógrafo se desplaza antes y después del encuentro de 2029, codificando el cambio en excentricidad, y valida que las velocidades N-cuerpos en tramos keplerianos caen sobre el círculo teórico.

### NB6 — CRTBP, Constante de Jacobi y Parámetro de Tisserand
**Archivo:** [`analisis/nb6_crtbp_jacobi_tisserand.ipynb`](analisis/nb6_crtbp_jacobi_tisserand.ipynb)

Dos bloques en un solo notebook (misma física de tres cuerpos):

**Bloque A — CRTBP:** Calcula la constante de Jacobi $C_J = 2\Omega - v^2$ y grafica las curvas de velocidad cero del sistema Sol-Tierra. Determina si la "puerta" en $L_1$ o $L_2$ está abierta o cerrada para Apophis: ¿puede ser capturado gravitacionalmente?

**Bloque B — Tisserand:** Calcula $T_P = a_T/a + 2\cos i\,\sqrt{(a/a_T)(1-e^2)}$ en función del tiempo (2027–2031) y demuestra que, aunque $a$ y $e$ cambian drásticamente tras el encuentro, $T_P$ se conserva ($\Delta T_P/T_P < 1\%$). Es la "prueba de identidad" orbital de Apophis.

### NB7 — Ecuaciones de Gauss: Perturbaciones Orbitales
**Archivo:** [`analisis/nb7_perturbaciones_gauss.ipynb`](analisis/nb7_perturbaciones_gauss.ipynb)

Modela cómo la fuerza perturbadora de la Tierra cambia la inclinación $i$ y el nodo ascendente $\Omega$ de Apophis durante el encuentro, usando las ecuaciones de Gauss-Lagrange (VOP):

$$\frac{di}{dt} = \frac{r\cos(\omega+f)}{h}\,W \qquad \frac{d\Omega}{dt} = \frac{r\sin(\omega+f)}{h\sin i}\,W$$

Compara los resultados con la integración N-cuerpos directa y menciona las variables de Delaunay como el formalismo canónico equivalente.

> Nota: `cuadraturas_PRelativo2cuerpos.ipynb` corresponde a un **cuaderno de clase** sobre el problema relativo de dos cuerpos y determinación orbital, que sirve como base teórica para los análisis posteriores.

---

## Contexto teórico

El problema físico es la dinámica gravitacional de Apophis dentro del Sistema Solar, en particular durante la ventana del **encuentro cercano con la Tierra (abril de 2029)**.  
Aunque la atracción dominante es la del Sol, para estimar con fidelidad la fecha y la distancia del máximo acercamiento se incorporan perturbaciones de cuerpos relevantes (Tierra–Luna, Venus, Júpiter y, en el modelo final, también Mercurio, Marte y Saturno).

### 1) Ecuaciones del problema de N‑cuerpos
El movimiento de cada cuerpo \(i\) está gobernado por la suma de las fuerzas gravitacionales de los demás cuerpos:

$$
\ddot{\mathbf r}_i \;=\; -G\sum_{j\neq i} m_j \frac{\mathbf r_i-\mathbf r_j}{\lVert \mathbf r_i-\mathbf r_j\rVert^3}.
$$

En este proyecto se usan **unidades canónicas** astronómicas (AU, $$M_\odot\$$, UT) y se trabaja con \(G = 1\), como se define explícitamente en los notebooks (por ejemplo en el Modelo 2 y el Modelo 3).  
Modelo 2: [`modelos/modelo2_SEMAJV.ipynb`](modelos/modelo2_SEMAJV.ipynb)  
Modelo 3: [`modelos/modelo3_completo.ipynb`](modelos/modelo3_completo.ipynb)

### 2) Magnitud clave observacional: distancia Tierra–Apophis
La cantidad que se minimiza para encontrar el máximo acercamiento es:

$$
d_{EA}(t)=\left\|\mathbf r_{\text{Apophis}}(t)-\mathbf r_{\text{Tierra}}(t)\right\|.
$$

Todos los modelos reportan el mínimo global de \(d_{EA}(t)\) en la ventana 2028–2029, y luego ejecutan un refinamiento local para aumentar la resolución temporal alrededor del mínimo.  
Modelo 1: [`modelos/modelo1_SEMAJ.ipynb`](modelos/modelo1_SEMAJ.ipynb)  
Modelo 2: [`modelos/modelo2_SEMAJV.ipynb`](modelos/modelo2_SEMAJV.ipynb)  
Modelo 3: [`modelos/modelo3_completo.ipynb`](modelos/modelo3_completo.ipynb)

---

## Metodología 
La metodología se construyó como un **pipeline único** que se repite en los tres notebooks, aumentando el número de cuerpos (y por tanto la física capturada) a medida que el proyecto avanza:

1. **Definir el “sistema” a integrar (qué cuerpos entran al modelo).**  
   El punto de partida fue un modelo reducido que captura lo esencial del encuentro (Sol + Tierra–Luna + Apophis), y luego se añadieron perturbadores relevantes (primero Júpiter; luego Venus; y finalmente Mercurio, Marte y Saturno).  
   La lógica detrás de este crecimiento es que, en un encuentro cercano, la distancia mínima Tierra–Apophis puede ser sensible a pequeñas perturbaciones acumuladas, y por eso interesa ver **cómo converge** el resultado al aumentar la complejidad.

2. **Normalización del problema con unidades canónicas.**  
   Para trabajar con un integrador estable y comparable entre modelos, se adoptó el sistema canónico astronómico donde \(G=1\). En estas unidades:
   - las posiciones quedan en AU,
   - las masas se expresan como razones respecto al Sol (masas “canónicas”),
   - y el tiempo se mide en UT (\(\approx 58.12\) días).  
   Esta normalización reduce escalas extremas y hace que los parámetros (como masas y velocidades) estén en un rango numéricamente más cómodo.

3. **Construcción de condiciones iniciales a partir de efemérides reales.**  
   Para que los resultados tengan interpretación física, los estados iniciales (posición y velocidad) se toman de **NASA Horizons** en la época
   \(\;t_0 =\) **2028-01-01**, en marco **baricéntrico** del Sistema Solar (`location='@0'`).  
   Las velocidades, que Horizons entrega en AU/día, se convierten a AU/UT para ser coherentes con la integración en unidades canónicas.  
   Con esto, cada modelo arranca desde el mismo “punto de partida” temporal y dinámico, cambiando únicamente el conjunto de perturbadores incluidos.

4. **Integración numérica N‑cuerpos en una ventana amplia (2028 → 2029).**  
   Con las masas canónicas y las condiciones iniciales, se resuelve el sistema de ecuaciones diferenciales del problema de N‑cuerpos durante una ventana de meses que cubre el encuentro (aprox. 2028-01-01 → 2029-07-01).  
   El objetivo de esta corrida “base” es localizar de manera robusta **en qué región temporal ocurre el mínimo** de la distancia Tierra–Apophis.

5. **Estimación del máximo acercamiento (mínimo de distancia).**  
   A partir de las trayectorias integradas, se calcula \(d_{EA}(t)\) en cada instante y se toma su mínimo global. Esto produce una primera estimación de:
   - fecha aproximada del encuentro,
   - distancia mínima,
   - y distancia en radios terrestres (cuando se reporta).

6. **Refinamiento local (reintegración con paso fino).**  
   Como el mínimo puede ocurrir en una ventana temporal corta, se ejecuta una segunda integración centrada en el mínimo hallado, típicamente en una ventana de **±30 días**, con un paso temporal mucho más pequeño (del orden de minutos).  
   Esta etapa mejora la precisión de la fecha y de la distancia mínima, y es la que se toma como resultado “final” para cada modelo.

7. **Validación de calidad numérica.**  
   Para asegurar que el resultado no es un artefacto del paso temporal o del integrador, se verifica la **conservación de energía** (y se discute momento angular como invariante). Una buena conservación es un indicador práctico de estabilidad del experimento numérico en el intervalo simulado.

8. **(Modelo 3) Convergencia e incertidumbre.**  
   En el Modelo 3 se consolidan dos extensiones del pipeline:
   - **Convergencia:** se comparan modelos crecientes (SEA → SEMA → SEMAJ → SEMAJV → 9 cuerpos) para cuantificar cuánto cambia \(d_{\min}\) al añadir perturbadores.
   - **Incertidumbre (Monte Carlo):** se evalúa la sensibilidad del resultado ante variaciones en condiciones iniciales de Apophis.
---

## Modelos del proyecto

### Modelo 1: SEMAJ (5 cuerpos)
- **Sistema:** Sol, Tierra, Luna, Júpiter, Apophis  
- **Objetivo:** estimación inicial del máximo acercamiento con el conjunto mínimo de perturbadores dominantes.
Notebook: [`modelos/modelo1_SEMAJ.ipynb`](modelos/modelo1_SEMAJ.ipynb)

### Modelo 2: SEMAJV (6 cuerpos)
- **Sistema:** Sol, Tierra, Luna, Júpiter, Venus, Apophis  
- **Qué mejora:** añade Venus como perturbador adicional y refuerza el pipeline (unidades, condiciones iniciales, integración, validación y refinamiento).
Notebook: [`modelos/modelo2_SEMAJV.ipynb`](modelos/modelo2_SEMAJV.ipynb)

### Modelo 3 (final): 9 cuerpos + convergencia + incertidumbre
- **Sistema de 9 cuerpos:** Sol, Mercurio, Venus, Tierra, Luna, Marte, Júpiter, Saturno, Apophis  
- **Bloques principales del cuaderno:**
  1. integración del modelo final de 9 cuerpos y estimación del máximo acercamiento,
  2. convergencia de \(d_{\min}\) al añadir cuerpos en cadena (SEA → SEMA → SEMAJ → SEMAJV → 9C),
  3. análisis de incertidumbre (Monte Carlo) sobre condiciones iniciales de Apophis.
Notebook: [`modelos/modelo3_completo.ipynb`](modelos/modelo3_completo.ipynb)

---

## Visualizaciones
El repositorio incluye visualizaciones para interpretar el evento:
- curvas de distancia Tierra–Apophis,
- plano XY del sistema completo,
- y una animación 2D centrada en la Tierra (HTML) para el acercamiento local.

Estas salidas están desarrolladas principalmente en el Modelo 3:  
[`modelos/modelo3_completo.ipynb`](modelos/modelo3_completo.ipynb)

---

## Dependencias principales
- `pymcel` (integración N‑cuerpos y utilidades de efemérides)  
  https://github.com/seap-udea/pymcel

---

## Orbit Viewer (nuevo)

Se agregó una implementación reusable tipo **VOP/JPL** en:

- `orbit_viewer.py` — módulo reutilizable (librería)
- `demo_orbit_viewer.ipynb` — notebook ejecutable con el demo completo

### Qué incluye

- API con `OrbitElements`, `propagate_two_body`, `orbital_plane_normal`, `make_orbit_viewer`.
- Propagación Kepleriana **two-body** usando `pymcel.kepler_newton`.
- Visualización 3D interactiva con:
  - órbitas de Mercurio, Venus, Tierra, Marte y **Apophis** con sus colores,
  - **tick marks** (picket-fence) que muestran la inclinación de cada plano orbital,
  - **posición actual** de cada cuerpo al 2026-04-22 00:00 UTC,
  - plano de la eclíptica con grilla (`z=0`),
  - línea de nodos (amarilla) y vector normal al plano orbital (verde) de Apophis,
  - banda sombreada entre eclíptica y plano orbital de Apophis,
  - anotación de distancias (Sol, Tierra) y fecha.

### Ejecutar demo

Abrir y ejecutar `demo_orbit_viewer.ipynb` en JupyterLab, Jupyter Notebook o Google Colab.

La celda final genera una figura interactiva en 3D embebida en el notebook.
