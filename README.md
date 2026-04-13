# Apophis-Asteroid-Project
Proyecto de Mecánica Celeste que estudia la aproximación a la Tierra del asteroide Apophis en abril del 2029.

## Modelo SEMAJ (.ipynb):
Modelo Sol-Tierra-Luna-Apopphis_Júpiter
### Metodología

1. **Condiciones iniciales**: obtenidas de NASA Horizons el 2028-01-01 referidas al baricentro del sistema solar.
2. **Integración N-cuerpos**: `pymcel.ncuerpos_solucion()` resuelve las ecuaciones de Newton con G = 1 (unidades canónicas).
3. **Unidades canónicas**: AU · M☉ · UT, donde UT = √(AU³/GM☉) ≈ 58.12 días ≈ año/2π.
4. **Refinamiento**: segunda integración con paso fino (≈ minutos) en la ventana del encuentro.
5. **Validación**: comprobación de conservación de energía y momento angular.

### Próximos pasos

- Hallar el valor de h para Apophis respecto al Sol.
- Hallas las cuadraturas de Apophis y verificar que permanezcan constantes.
- Cuantificar la convergencia de la distancia mínima al ir añadiendo cuerpos.
- Análisis de incertidumbre: Monte Carlo sobre las condiciones iniciales de Apophis.
