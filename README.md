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

El repositorio incluye además el módulo de investigación **"DART Inverso"**: diseño analítico de una maniobra de impacto cinético sobre Apophis para desviar su órbita de manera que impacte la Tierra, cubriendo desde el escape geocéntrico del impactor hasta el encuentro hiperbólico final.

---

## Estructura del repositorio

```
Apophis-Asteroid-Project/
├── modelos/                              ← Integración N-cuerpos (pipeline principal)
│   ├── modelo1_SEMAJ.ipynb
│   ├── modelo2_SEMAJV.ipynb
│   └── modelo3_completo.ipynb
├── analisis/                             ← Análisis avanzados (mecánica celeste profunda)
│   ├── nb4_kepler_benchmark.ipynb
│   ├── nb5_hodografo.ipynb
│   ├── nb6_crtbp_jacobi_tisserand.ipynb
│   └── nb7_perturbaciones_gauss.ipynb
├── dart_inverso/                         ← NUEVO: Maniobra cinética "DART Inverso"
│   ├── nb00_parches_conicos.ipynb
│   ├── nb01_kepler_solver_clase.ipynb
│   ├── nb02_cuadraturas_invariantes.ipynb
│   ├── nb03_hodografo_velocidad.ipynb
│   ├── nb04_sensibilidad_montecarlo.ipynb
│   └── nb05_gauss_lagrange_perturbaciones.ipynb
├── Orbit_Viewer.ipynb
├── cuadraturas_PRelativo2cuerpos.ipynb
└── README.md
```

### Flujo narrativo recomendado

```
cuadraturas (existente) → nb4_kepler → nb5_hodografo → nb6_jacobi_tisserand → nb7_gauss
       ↑                                                          ↑
  modelo3_completo (existente) ─────────────────────────────────┘

         │
         ▼  (continúa en dart_inverso/)
  nb00_parches_conicos → nb01_kepler_solver → nb02_cuadraturas_invariantes
                              ↓
              nb03_hodografo → nb04_sensibilidad → nb05_gauss_lagrange
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

## Módulo "DART Inverso" (`dart_inverso/`)

Esta carpeta contiene la investigación para diseñar y visualizar una maniobra de **impacto cinético sobre Apophis** con el objetivo de desviar su órbita heliocéntrica de forma que impacte la Tierra. El problema se aborda mediante la **jerarquía de parches cónicos**: una sucesión de sistemas de 2 cuerpos relativos que modela cada fase de la maniobra.

### Arquitectura del Modelo: Jerarquía de Parches Cónicos

La misión se descompone en tres fases que corresponden a sistemas dinámicos distintos:

| Fase | Sistema dominante | Centro de masa | Condición de éxito |
|------|------------------|----------------|-------------------|
| **Fase 0** — Escape del impactor | Tierra–Impactor | Tierra | $\varepsilon \geq 0$ (escape geocéntrico) |
| **Fase 1** — Intercepción heliocéntrica | Sol–Apophis | Sol (baricentro) | $\Delta V$ que lleva la nueva elipse de Apophis hasta la posición futura de la Tierra |
| **Fase 2** — Encuentro cercano | Tierra–Apophis | Tierra | Periapsis $q < R_\oplus$ de la hipérbola de encuentro |

La **transición entre sistemas** ocurre en la Esfera de Hill terrestre:

$$R_H = a_\oplus \left(\frac{M_\oplus}{3 M_\odot}\right)^{1/3} \approx 0.01 \text{ UA} \approx 1.5 \times 10^6 \text{ km}$$

Cuando Apophis cruza $R_H$, la gravedad del Sol pasa a ser perturbación y el problema dominante cambia de heliocéntrico a geocéntrico.

---

### NB00 — Arquitectura y Parches Cónicos

**Archivo:** [`dart_inverso/nb00_parches_conicos.ipynb`](dart_inverso/nb00_parches_conicos.ipynb)

**Descripción física:**  
Notebook de planificación y diagnóstico. Calcula los parámetros orbitales de referencia para el diseño de la maniobra:

- **Fase 0:** A partir del vector de estado geocéntrico del impactor, calcula la velocidad de escape geocéntrica $v_{esc} = \sqrt{2\mu_\oplus / r}$ y verifica que la energía específica $\varepsilon = v^2/2 - \mu/r \geq 0$. Determina el radio de la Esfera de Hill y el tiempo de tránsito para cruzarla.
- **Fase 1 (diagnóstico):** Obtiene el estado heliocéntrico de Apophis en el instante del impacto (usando Horizons), aplica la Ecuación de Vis-Viva $v^2 = \mu_\odot(2/r - 1/a)$ para comparar la velocidad actual con la velocidad requerida en la nueva órbita, y estima el $\Delta V$ de primer orden.
- **Fase 2 (diagnóstico):** Calcula la velocidad hiperbólica de Apophis al cruzar $R_H$ y estima el periapsis $q$ de la hipérbola resultante mediante $q = a_h(e_h - 1)$, verificando si $q < R_\oplus$.

**Funciones `pymcel` a usar/extender:**

| Función | Uso |
|---------|-----|
| `pc.consulta_horizons(id='99942', epochs=..., datos='vectors')` | Estado heliocéntrico de Apophis en $t_{impacto}$ |
| `pc.consulta_horizons(id='399', ...)` | Estado heliocéntrico de la Tierra en $t_{impacto}$ y $t_{llegada}$ |
| `pc.estado_a_elementos(mu, estado)` | Obtener $(p, e, i, \Omega, \omega, f)$ de Apophis antes del impacto |
| `pc.unidades_canonicas(...)` | Definir el sistema de unidades para cada fase |
| `from pymcel.constantes import GM_sun, GM_earth, R_earth, au` | Constantes físicas |

---

### NB01 — Clase `KeplerSolver`: Benchmarking de Métodos

**Archivo:** [`dart_inverso/nb01_kepler_solver_clase.ipynb`](dart_inverso/nb01_kepler_solver_clase.ipynb)

**Descripción física:**  
Implementa una clase `KeplerSolver` orientada a objetos que encapsula tres métodos para resolver la ecuación trascendental de Kepler $M = E - e \sin E$:

- **Newton-Raphson:** Convergencia cuadrática. Iteración $E_{n+1} = E_n - (E_n - e\sin E_n - M)/(1 - e\cos E_n)$. Sensible a la elección del punto inicial cerca del perigeo ($E \approx 0$).
- **Punto Fijo:** Convergencia lineal. Iteración $E_{n+1} = M + e\sin E_n$. Diverge para $e > 1$ y converge lentamente para $e$ grande.
- **Laguerre-Conway ($n = 5$):** Convergencia de orden superior. Robusto globalmente, recomendado para órbitas muy excéntricas.

El benchmark evalúa:
1. Número de iteraciones hasta convergencia ($\delta < 10^{-12}$) en función de la anomalía media $M \in [0, 2\pi]$.
2. Tiempo de CPU por llamada usando `timeit`.
3. Comportamiento cerca del perigeo ($M \to 0$, $e = 0.1914$ para Apophis) donde Newton-Raphson puede necesitar más iteraciones.

La clase también implementa la solución para órbitas hiperbólicas ($e > 1$) mediante la ecuación de Kepler hiperbólica $M_h = e \sinh F - F$.

**Funciones `pymcel` a usar/extender:**

| Función | Uso |
|---------|-----|
| `pc.kepler_newton(M, e, G0, delta)` | Referencia: método Newton ya implementado en pymcel |
| `pc.kepler_semianalitico(M, e)` | Comparación adicional |
| `pc.kepler_eserie(M, e, orden)` | Comparación: series de Fourier |
| `pc.kepler_bessel(M, e, delta)` | Comparación: método de Bessel |
| `pc.metodo_newton(f, x0, delta, args)` | Método genérico de raíz para extender a hiperbólica |
| `pc.metodo_laguerre(f, x0, delta, args, eta=5)` | Base del método Laguerre-Conway |

> **Extensión propuesta:** La clase `KeplerSolver` encapsula `pc.kepler_newton`, `pc.kepler_semianalitico` y `pc.metodo_laguerre` bajo una interfaz unificada `solver.solve(M, e, method='newton'|'fixed_point'|'laguerre')`.

---

### NB02 — Cuadraturas e Invariantes del Movimiento (M.A.R.E.)

**Archivo:** [`dart_inverso/nb02_cuadraturas_invariantes.ipynb`](dart_inverso/nb02_cuadraturas_invariantes.ipynb)

**Descripción física:**  
Valida numéricamente la conservación de las tres constantes primeras del movimiento kepleriano durante la propagación completa de la maniobra, tanto en la fase heliocéntrica (Sol-Apophis) como en la fase geocéntrica (Tierra-Apophis):

1. **M.A.R.E. (Momento Angular Relativo Específico):**  
   $\vec{h} = \vec{r} \times \dot{\vec{r}}$, con $|\vec{h}| = \sqrt{\mu p}$. Se verifica $\|\vec{h}(t) - \vec{h}_0\| / \|\vec{h}_0\| < 10^{-8}$ a lo largo de toda la trayectoria.

2. **E.R.E. (Energía Específica Relativa):**  
   $\varepsilon = \frac{v^2}{2} - \frac{\mu}{r}$. Se muestra que $\varepsilon$ es constante en cada arco kepleriano y cambia discretamente solo en el instante del impulso ($\Delta V$).

3. **Vector de Laplace-Runge-Lenz:**  
   $\vec{e} = \frac{\dot{\vec{r}} \times \vec{h}}{\mu} - \hat{r}$, con $|\vec{e}| = e$ (excentricidad escalar). Este vector apunta al periapsis y su conservación verifica que no hay precesión orbital en el tramo kepleriano.

El notebook también implementa las funciones $f$ y $g$ de Lagrange para propagar el estado sin integración numérica:

$$\vec{r}(t) = f(t, t_0)\,\vec{r}_0 + g(t, t_0)\,\dot{\vec{r}}_0$$

y las usa para obtener el estado de Apophis en un tiempo futuro dado únicamente el estado inicial y $\Delta t$.

**Funciones `pymcel` a usar/extender:**

| Función | Uso |
|---------|-----|
| `pc.estado_a_elementos(mu, estado)` | Extraer elementos en cada instante para verificar conservación |
| `pc.elementos_a_estado(mu, elementos)` | Reconstrucción para validar round-trip |
| `pc.doscuerpos_solucion(mu, r, v, ts)` | Propagación kepleriana de referencia |
| `pc.kepler_newton(M, e, delta=1e-12)` | Resolver $E$ en cada paso de la propagación $f$, $g$ |
| `pc.C(z)`, `pc.S(z)` | Series de Stumpff para las funciones $f$ y $g$ universales |
| `pc.funcion_universal_kepler(x, M, e, q)` | Anomalía universal para trayectorias elípticas e hiperbólicas |

> **Extensión propuesta:** Implementar `propagar_fg(mu, r0, v0, dt)` usando las series de Stumpff (`pc.C`, `pc.S`) para cubrir el caso hiperbólico de la Fase 2 sin cambiar de formulación.

---

### NB03 — Hodógrafo de Velocidad: Geometría del Impacto

**Archivo:** [`dart_inverso/nb03_hodografo_velocidad.ipynb`](dart_inverso/nb03_hodografo_velocidad.ipynb)

**Descripción física:**  
Visualiza el hodógrafo de Apophis (la curva que describe el extremo del vector velocidad en el espacio $v_x$-$v_y$) antes y después del impulso cinético. El Teorema de Hamilton establece que este lugar geométrico es siempre un **círculo** para órbitas keplerianas, con:

- **Radio del círculo:** $\mu_\odot / h$
- **Desplazamiento del centro respecto al origen:** $\mu_\odot \, e / h$ (en la dirección del perihelio)

El notebook muestra en un solo gráfico:
1. El hodógrafo de Apophis **pre-impacto** (círculo de radio $\mu/h_0$, desplazado $\mu e_0/h_0$).
2. El hodógrafo de Apophis **post-impacto** (círculo de radio $\mu/h_f$, desplazado $\mu e_f/h_f$), donde $h_f$ y $e_f$ corresponden a la nueva órbita que pasa por la posición futura de la Tierra.
3. El vector $\Delta \vec{V}$ que conecta el punto de velocidad en el instante del impacto entre los dos círculos.
4. Las velocidades N-cuerpos de `modelo3_completo.ipynb` superpuestas, para verificar que los tramos keplerianos caen sobre los círculos teóricos.

El hodógrafo es la representación dual a las cónicas en el espacio de posición; permite ver de forma inmediata el costo energético de la maniobra y la dirección óptima del impulso.

**Funciones `pymcel` a usar/extender:**

| Función | Uso |
|---------|-----|
| `pc.consulta_horizons(id='99942', epochs=..., datos='vectors')` | Velocidades heliocentricas de Apophis en la ventana de interés |
| `pc.estado_a_elementos(mu, estado)` | Obtener $h$, $e$ para construir los dos círculos |
| `pc.doscuerpos_solucion(mu, r, v, ts)` | Generar la trayectoria kepleriana para superponer al hodógrafo |
| `pc.fija_ejes_proporcionales(ax, datos)` | Ejes proporcionales en la figura $v_x$-$v_y$ |

---

### NB04 — Sensibilidad y Antelación: Simulación Monte Carlo

**Archivo:** [`dart_inverso/nb04_sensibilidad_montecarlo.ipynb`](dart_inverso/nb04_sensibilidad_montecarlo.ipynb)

**Descripción física:**  
Evalúa cómo el **tiempo de antelación del impacto** (Lead Time, $\tau$) afecta la magnitud del $\Delta V$ necesario para que Apophis impacte la Tierra. Cuanto antes se golpea a Apophis, más pequeña es la perturbación necesaria porque la desviación tiene más tiempo para acumularse.

**Metodología:**
1. Para cada Lead Time $\tau \in [10, 30, 100, 365, 3 \times 365]$ días antes del encuentro previsto:
   - Calcular el estado de Apophis en $t_{impacto} = t_{encuentro} - \tau$ usando las funciones $f$ y $g$ (propagación analítica hacia atrás).
   - Resolver el problema de Lambert inverso: encontrar la velocidad $\vec{V}_{nueva}$ en $t_{impacto}$ tal que la nueva trayectoria pase por $\vec{r}_{Tierra}(t_{encuentro})$.
   - Calcular $\Delta V = |\vec{V}_{nueva} - \vec{V}_{Apophis}(t_{impacto})|$.

2. **Monte Carlo:** Para cada $\tau$, agregar perturbaciones gaussianas $\sigma_r = 1$ km, $\sigma_v = 1$ mm/s sobre las condiciones iniciales de Apophis (incertidumbre observacional) y evaluar la distribución de $\Delta V$ resultante ($N_{MC} = 500$ realizaciones).

3. Graficar $\Delta V$ vs $\tau$ (escala logarítmica) con bandas de incertidumbre $\pm 1\sigma$, mostrando la ley de escala $\Delta V \propto \tau^{-2}$ predicha por la teoría de deflexión orbital.

**Funciones `pymcel` a usar/extender:**

| Función | Uso |
|---------|-----|
| `pc.solucion_lambert(P1, P2, tf, mu, direccion)` | Resolver la transferencia Apophis → posición futura de la Tierra |
| `pc.consulta_horizons(id='399', epochs=...)` | Efeméride de la Tierra en $t_{encuentro}$ |
| `pc.consulta_horizons(id='99942', epochs=...)` | Estado de Apophis en el instante del impacto para cada $\tau$ |
| `pc.doscuerpos_solucion(mu, r, v, ts)` | Propagación analítica alternativa (validación) |
| `pc.estado_a_elementos(mu, estado)` | Verificar que la nueva órbita es la correcta |
| `pc.funcion_universal_kepler(x, M, e, q)` | Propagación universal para cada realización Monte Carlo |

> **Nota:** Las funciones $f$ y $g$ implementadas en NB02 son la base de la propagación hacia atrás en este notebook.

---

### NB05 — Ecuaciones de Gauss-Lagrange: Perturbaciones Orbitales

**Archivo:** [`dart_inverso/nb05_gauss_lagrange_perturbaciones.ipynb`](dart_inverso/nb05_gauss_lagrange_perturbaciones.ipynb)

**Descripción física:**  
Implementa las **Ecuaciones de Variación de Parámetros (VOP) de Gauss/Lagrange** para modelar cómo el impulso cinético ($\Delta V$ aplicado como aceleración impulsiva) altera los elementos orbitales de Apophis, con énfasis en el semieje mayor $a$ y el movimiento medio $n = \sqrt{\mu/a^3}$.

Las ecuaciones de Gauss en la base radial-transversal-normal $(R, T, W)$ son:

$$\frac{da}{dt} = \frac{2a^2}{\sqrt{\mu p}} \left[ e \sin f \cdot F_R + \frac{p}{r} F_T \right]$$

$$\frac{de}{dt} = \frac{\sqrt{p/\mu}}{1} \left[ \sin f \cdot F_R + \left(\cos f + \frac{e + \cos f}{1 + e\cos f}\right) F_T \right]$$

$$\frac{di}{dt} = \frac{r \cos(\omega + f)}{h} F_W \qquad \frac{d\Omega}{dt} = \frac{r \sin(\omega + f)}{h \sin i} F_W$$

donde $(F_R, F_T, F_W)$ son las componentes de la fuerza perturbadora por unidad de masa en la base orbital.

El notebook:
1. Modela el impulso como una fuerza impulsiva $\vec{F} = \Delta\vec{V} \cdot \delta(t - t_{impacto})$ y calcula la variación instantánea $\Delta a$, $\Delta e$, $\Delta i$ mediante integración analítica de las ecuaciones de Gauss.
2. Integra numéricamente las ecuaciones de Gauss en el período siguiente al impacto para ver la evolución secular de $a(t)$ y $n(t) = \sqrt{\mu/a(t)^3}$.
3. Compara con la integración directa N-cuerpos del Modelo 3.
4. Menciona las **Variables de Delaunay** como el formalismo canónico equivalente:
   - $L = m\sqrt{\mu a}$ (conjugada a $M$), $G = m\sqrt{\mu a(1-e^2)}$ (conjugada a $\omega$), $H = G\cos i$ (conjugada a $\Omega$).

**Funciones `pymcel` a usar/extender:**

| Función | Uso |
|---------|-----|
| `pc.estado_a_elementos(mu, estado)` | Obtener elementos orbitales en cada instante |
| `pc.elementos_a_estado(mu, elementos)` | Reconstrucción tras integrar las VOP |
| `pc.doscuerpos_solucion(mu, r, v, ts)` | Propagación de referencia (sin perturbación) |
| `pc.ncuerpos_solucion(sistema, ts)` | Solución N-cuerpos para comparación directa |
| `pc.consulta_horizons(id='99942', epochs=...)` | Estado real de Apophis antes y después del encuentro |

---

### Resumen de funciones `pymcel` por módulo

| Módulo `pymcel` | Funciones clave | Notebooks que la usan |
|----------------|-----------------|----------------------|
| **Efemérides** | `consulta_horizons`, `consulta_spice`, `carga_kernels` | NB00, NB01, NB03, NB04 |
| **Kepler** | `kepler_newton`, `kepler_semianalitico`, `kepler_bessel`, `kepler_eserie`, `metodo_laguerre` | NB01, NB02 |
| **Dos cuerpos** | `doscuerpos_solucion`, `estado_a_elementos`, `elementos_a_estado` | NB00, NB02, NB03, NB04, NB05 |
| **N-cuerpos** | `ncuerpos_solucion`, `ncuerpos_a_pandas`, `plot_ncuerpos_3d` | NB05 |
| **Lambert** | `solucion_lambert` | NB04 |
| **CRTBP** | `crtbp_solucion`, `constante_jacobi`, `funcion_puntos_colineales` | NB00 (diagnóstico esfera de Hill) |
| **Stumpff/Universal** | `C(z)`, `S(z)`, `funcion_universal_kepler` | NB02 |
| **Visualización** | `fija_ejes_proporcionales`, `fija_ejes3d_proporcionales`, `plot_ncuerpos_3d`, `plotly_esfera` | NB03 |
| **Constantes** | `GM_sun`, `GM_earth`, `R_earth`, `au`, `year`, `M_sun` | Todos |

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
