# Guion — Defensa TFG (25 minutos)

Distribución objetivo: ~23 min de presentación + 2 min de colchón/transiciones.

Cada observable de la sección de Resultados y de la Comparativa dinámica se muestra ahora en **tres diapositivas grandes** (`fullsizeslides`), una por algoritmo (Metropolis, Glauber, Wolff), sin diapositiva de resumen previa con las tres gráficas pequeñas. Toda la explicación se da sobre la **primera diapositiva (Metropolis)**, que incluye la gráfica algo más pequeña para dejar sitio a los puntos clave debajo; las diapositivas de Glauber y Wolff que siguen se pasan rápido (10–15 s), solo señalando que el comportamiento es análogo, salvo que quieras remarcar algo concreto — se anota en cada bloque con *(Glauber/Wolff — pasar rápido)*.

---

## [Portada] — 0:00–0:30 (30 s)

> "Buenos días. Mi trabajo de fin de grado versa sobre el modelo de Ising 2D y la comparativa de tres algoritmos Montecarlo: Metropolis, Glauber y Wolff."

---

## [Contenidos] — 0:30–1:00 (30 s)

> "La presentación tiene seis bloques: introducción al modelo, fundamentos teóricos, los tres algoritmos, la configuración de la simulación, los resultados de equilibrio y la comparativa dinámica, para acabar con las conclusiones."

---

## SECCIÓN 1 — Introducción y motivación

### [El modelo de Ising] — 1:00–2:30 (1:30)

> "El modelo de Ising es el banco de pruebas clásico de la física estadística computacional. Fue propuesto por Lenz en 1920 y resuelto en 1D por Ising en 1925. La solución exacta en 2D la dio Onsager en 1944, y predice una transición de fase ferromagnética–paramagnética de segundo orden."
>
> "El hamiltoniano es sencillo: solo interacción entre vecinos más próximos con constante de acoplamiento J. Con J=k_B=1, la temperatura crítica exacta es T_c ≈ 2,2692, y ese valor exacto es lo que nos permite validar todos los resultados numéricos."

### [Objetivos] — 2:30–3:30 (1:00)

> "Los objetivos son cuatro: implementar los tres algoritmos, calcular los observables de equilibrio para cuatro tamaños de red —L de 16 a 128—, aplicar escalado de tamaño finito para extraer exponentes críticos y T_c, y comparar la dinámica mediante el tiempo de autocorrelación integrado τ_int."

---

## SECCIÓN 2 — Fundamentos teóricos

### [Transición de fase de segundo orden] — 3:30–4:15 (45 s)

> "Antes de ver los observables concretos, conviene fijar el marco teórico. En una transición de segundo orden como la del Ising, el parámetro de orden se anula de forma continua en T_c y la longitud de correlación diverge como |t|^{−ν}. Esa divergencia arrastra la de la susceptibilidad y la capacidad calorífica. Los exponentes ν, γ, α… caracterizan completamente la transición, y medirlos numéricamente y contrastarlos con los valores exactos de Onsager es precisamente el objetivo del trabajo."

### [Observables del sistema] — 4:15–5:30 (1:15)

> "Los cinco observables que medimos son: la magnetización, que es el parámetro de orden; la capacidad calorífica, que recoge las fluctuaciones de energía; la susceptibilidad, que mide la respuesta ante un campo externo y diverge en T_c; el cumulante de Binder, que es adimensional y tiene la propiedad de que sus curvas para distintos L se cruzan en T_c; y la longitud de correlación ξ_L, que extraemos del factor de estructura mediante el estimador de Caracciolo-Sokal."

### [Clase de universalidad y FSS] — 5:30–6:45 (1:15)

> "Los exponentes críticos no dependen de los detalles microscópicos del modelo, sino solo de su dimensión y la simetría del parámetro de orden. A eso se llama clase de universalidad. El Ising 2D tiene valores exactos: ν=1, γ/ν=7/4, β/ν=1/8 y α=0."
>
> "En una red finita ξ no puede divergir: L hace su papel. Con la variable de escala x = (T−T_c)·L, las curvas de χ·L^{−7/4} y ξ_L/L colapsan sobre funciones universales. Esto nos permite extraer los exponentes incluso en redes pequeñas."

---

## SECCIÓN 3 — Algoritmos Montecarlo

### [Montecarlo: muestreo por importancia] — 6:45–7:45 (1:00)

> "El problema fundamental es que la función de partición tiene 2^N términos, inabarcable para N grande. La solución es el muestreo por importancia: construimos una cadena de Markov que converge a la distribución de Boltzmann. La condición suficiente para garantizarlo es el balance detallado, que equilibra los flujos de probabilidad entre estados."

### [Metropolis y Glauber] — 7:45–9:15 (1:30)

> "Metropolis y Glauber son algoritmos locales: en cada paso se elige un espín al azar y se propone invertirlo. Metropolis acepta con min(1, e^{−βΔE}). Glauber asigna directamente la orientación con una probabilidad sigmoidal. Ambos satisfacen el balance detallado, pero Glauber tiene una probabilidad de aceptación sistemáticamente menor, lo que lo hace entre 1,5 y 2 veces más lento. El exponente dinámico z ≈ 2,17 es el mismo para los dos."

### [Wolff: algoritmo de cúmulo] — 9:15–10:45 (1:30)

> "Wolff opera de forma muy diferente: en lugar de un espín, invierte un cúmulo entero. Se parte de una semilla aleatoria, se añaden vecinos con la misma orientación con probabilidad p = 1 − e^{−2βJ}, y se invierte todo el cúmulo de una vez. La ventaja clave es que cerca de T_c el cúmulo es grande y el algoritmo realiza movimientos globales, eliminando prácticamente el critical slowing down. El exponente dinámico cae a z ≈ 0,25."

---

## SECCIÓN 4 — Configuración de la simulación

### [Parámetros de la simulación] — 10:45–11:45 (1:00)

> "Simulamos L ∈ {16, 32, 64, 128}, con un factor de variación ×8 entre extremos, lo suficiente para ver el escalado. Arrancamos en frío, todos los espines alineados, para evitar metaestabilidad. Usamos casi dos millones de pasos de medida y una malla de temperaturas adaptativa: muy fina alrededor de T_c —con un espaciado de 0,1·T_c/L— y gruesa fuera de esa región."

---

## SECCIÓN 5 — Resultados

### [Capacidad calorífica c — Metropolis] — 11:45–12:30 (45 s)

> "La capacidad calorífica muestra el comportamiento esperado: el máximo se agudiza y se desplaza hacia T_c al crecer L. Pero dado que α=0 —divergencia logarítmica—, extraer T_c de aquí tiene poca precisión. Los tres algoritmos dan resultados prácticamente idénticos, como debe ser para un observable de equilibrio."

*(Glauber / Wolff — pasar rápido, mismo comportamiento)*

### [Magnetización — Metropolis] — 12:30–13:00 (30 s)

> "La magnetización muestra la transición de orden-desorden. La caída se vuelve más abrupta al aumentar L, aproximándose al salto discontinuo del límite termodinámico."

*(Glauber / Wolff — pasar rápido)*

### [Cumulante de Binder — Metropolis] — 13:00–13:45 (45 s)

> "El cumulante de Binder es mucho más limpio para estimar T_c: las curvas de los distintos tamaños se cortan en un punto único, T_c, con valor universal U_4^* ≈ 0,61, consistente con Onsager."

*(Glauber / Wolff — pasar rápido)*

### [Extrapolación de T_c] — 13:45–14:45 (1:00)

> "Para cuantificar T_c tomamos los cruces entre pares de tamaños y los extrapolamos a L→∞ graficando frente a 1/L^2. Solo Wolff muestra la convergencia sistemática esperada; Metropolis y Glauber están dominados por ruido estadístico para L=128, consecuencia directa del critical slowing down que veremos en la sección dinámica."

### [Susceptibilidad magnética — Metropolis] — 14:45–15:25 (40 s)

> "La susceptibilidad presenta el pico más claro de todos los observables. El máximo diverge con L siguiendo χ_max ~ L^{γ/ν}."

*(Glauber / Wolff — pasar rápido)*

### [Escalado de χ_max] — 15:25–16:10 (45 s)

> "En escala log-log la pendiente es ≈ 1,75, perfectamente compatible con el valor exacto γ/ν = 7/4. Esto es una verificación directa del exponente crítico."

### [Colapso FSS de la susceptibilidad — Metropolis] — 16:10–16:55 (45 s)

> "El colapso FSS es donde aparece la diferencia dinámica entre algoritmos. Metropolis y Glauber muestran el pico en x ≈ 1, mientras que Wolff lo tiene en x ≈ 1,5. La causa es que para L=128 las medidas de los algoritmos locales aún están correlacionadas —el sistema no ha decorrelacionado entre muestras—, lo que sesga la estimación del observable."

*(Glauber — pasar rápido, mismo patrón; Wolff — señalar el desplazamiento del pico a x≈1,5 al mostrarla)*

### [Longitud de correlación ξ_L — Metropolis] — 16:55–17:15 (20 s)

> "La longitud de correlación confirma el mismo patrón: los máximos crecen con L y se desplazan a T_c."

*(Glauber / Wolff — pasar rápido)*

### [Ratio ξ_L/L — Metropolis] — 17:15–17:35 (20 s)

> "El ratio ξ_L/L es especialmente útil porque sus curvas también se cruzan en T_c, igual que el cumulante de Binder pero a partir de un observable completamente independiente —es una comprobación de coherencia interna."

*(Glauber / Wolff — pasar rápido)*

### [Escalado de ξ_L] — 17:35–18:15 (40 s)

> "El escalado de ξ_max con L da pendiente ≈ 1 en log-log, directamente ν ≈ 1."

### [Colapso FSS de ξ_L/L — Metropolis] — 18:15–18:35 (20 s)

> "Y el colapso FSS de ξ_L/L con x = (T−T_c)·L^{1/ν} confirma ese valor sin ningún ajuste adicional."

*(Glauber / Wolff — pasar rápido)*

---

## SECCIÓN 6 — Comparativa dinámica

### [Tiempo de autocorrelación τ_int — Metropolis] — 18:35–19:50 (1:15)

> "Pasamos a la parte más distintiva del trabajo: la comparativa dinámica. El tiempo de autocorrelación integrado mide cuántos pasos Montecarlo hay que dar para obtener una muestra estadísticamente independiente. Para Metropolis y Glauber, τ_int tiene un pico enorme en T_c que crece con L —eso es el critical slowing down—. Glauber es entre 1,5 y 2 veces más lento que Metropolis en todo el rango. Wolff, en cambio, mantiene τ_int del orden de 1 en todo el rango de temperaturas, sin estructura apreciable."

*(Glauber — pasar rápido, mismo pico algo mayor; Wolff — señalar la ausencia de estructura al mostrarla)*

### [Exponente dinámico z] — 19:50–20:50 (1:00)

> "Cuantificamos el critical slowing down midiendo τ_max en función de L. La ley de potencias τ_max ~ L^z da z ≈ 2 para Metropolis y Glauber —consistente con la literatura, z ≈ 2,17— y z ≈ 0 para Wolff. La consecuencia práctica es brutal: para L=128, obtener una muestra independiente en T_c cuesta aproximadamente L^2 ≈ 16 000 veces más con algoritmos locales."

---

## SECCIÓN 7 — Conclusiones

### [Conclusiones] — 20:50–22:35 (1:45)

> "Las conclusiones se dividen en dos bloques."
>
> "En equilibrio: los tres algoritmos reproducen correctamente la física. La temperatura crítica estimada por cruces de U_4 es compatible con Onsager. Los exponentes γ/ν ≈ 7/4 y ν ≈ 1 sitúan los datos en la clase de universalidad del Ising 2D."
>
> "En dinámica: Metropolis y Glauber tienen z ≈ 2 y sufren critical slowing down severo cerca de T_c. Glauber es sistemáticamente más lento que Metropolis. Wolff prácticamente elimina el problema, con z ≈ 0. La correlación residual de los algoritmos locales tiene además un efecto visible en el colapso FSS para L grandes."
>
> "Como perspectiva natural está la extensión al Ising 3D, donde no existe solución analítica y los exponentes críticos solo se conocen numéricamente, por lo que la comparativa dinámica entre algoritmos tiene aún más relevancia práctica."

### [Gracias] — 22:35–23:05 (30 s)

> "Muchas gracias. Quedo a vuestra disposición para las preguntas."

---

**Tiempo total estimado: ~23:05**, dejando ~2 minutos de margen para pausas naturales y transiciones.

---

---

# Preguntas probables del tribunal

## Sobre el modelo y la teoría

1. **¿Por qué el cumulante de Binder se cruza exactamente en T_c y no depende de L?**
   — La clave está en que U_4 es una cantidad adimensional construida como cociente de momentos de m. En el punto crítico la dependencia en L se cancela exactamente, dejando solo el valor universal U_4^* ≈ 0,61.

2. **¿Qué significa exactamente que α = 0 para la capacidad calorífica?**
   — Que la divergencia es logarítmica, no algebraica. Es el exponente límite: c ~ log|T − T_c|. Por eso la extracción de T_c desde c es mucho menos precisa que desde χ o U_4.

3. **¿Cómo se obtiene la longitud de correlación ξ_L desde los datos de simulación?**
   — A partir del factor de estructura S(k): se mide S(0) y S(k_min) con k_min = 2π/L, y se aplica la fórmula ξ_L = [1/(2sin(π/L))]·sqrt(S(0)/S(k_min) − 1). Es el estimador de Caracciolo–Sokal.

4. **¿Por qué usáis condiciones de contorno periódicas?**
   — Para minimizar los efectos de frontera. Con condiciones periódicas la red no tiene borde y todos los espines ven el mismo número de vecinos, lo que reduce los artefactos de tamaño finito.

## Sobre los algoritmos

5. **¿Qué garantiza que la cadena de Markov converge a la distribución de Boltzmann?**
   — El balance detallado (condición suficiente) junto con ergodicidad (que cualquier configuración sea alcanzable). Los tres algoritmos satisfacen ambas.

6. **¿Por qué Glauber es más lento que Metropolis si tienen el mismo exponente dinámico z?**
   — Mismo z significa el mismo escalado con L, pero el prefactor es distinto. La probabilidad de aceptación de Glauber es sistemáticamente menor que la de Metropolis para ΔE < 0, por lo que hace más rechazos.

7. **¿Por qué Wolff define "1 paso MC" como girar N espines en total y no como un único cúmulo?**
   — Para que la unidad de tiempo sea comparable entre algoritmos. Un barrido de red con Metropolis toca N espines. Con Wolff, un cúmulo cerca de T_c puede ser muy grande o muy pequeño dependiendo de la temperatura; normalizar a N espines girados hace la comparación justa.

8. **¿Podría Wolff tener z < 0? ¿Qué significaría eso?**
   — Teóricamente no tiene sentido físico: z < 0 implicaría que τ_int *decrece* al aumentar L, lo que violaría la ergodización. En la práctica, z ≈ 0 significa que τ_int es esencialmente constante con L.

9. **¿Por qué solo hay 5 valores posibles de ΔE en Metropolis para el Ising 2D?**
   — Porque cada espín tiene 4 vecinos, y ΔE = 2J·σ_i·(suma de vecinos). La suma de los 4 vecinos solo puede tomar los valores −4, −2, 0, +2, +4, dando ΔE ∈ {−8, −4, 0, +4, +8} (en unidades de J). Esto permite precalcular las probabilidades de aceptación.

## Sobre la simulación y los resultados

10. **¿Cómo habéis estimado los errores estadísticos de los observables?**
    — Con el método de bootstrap o jackknife sobre las series temporales, teniendo en cuenta las correlaciones temporales mediante τ_int. Sin corregir por τ_int los errores quedarían subestimados, especialmente cerca de T_c con algoritmos locales.

11. **¿Por qué arrancáis en frío (todos los espines alineados) y no en caliente (aleatorio)?**
    — Para evitar quedar atrapados en estados metaestables en la fase ordenada. El arranque frío termaliza bien en todo el rango de T; el arranque caliente puede tener problemas de metaestabilidad por debajo de T_c para L grandes.

12. **¿Cuántos pasos de termalización usáis y cómo los determinasteis?**
    — Descartando la parte inicial de la serie temporal hasta que el observable se estabiliza. El criterio práctico es esperar al menos varios τ_int. Se puede verificar comparando estimadores calculados en la primera y la segunda mitad de la simulación.

13. **¿Por qué el colapso FSS de la susceptibilidad está desplazado para Metropolis/Glauber respecto a Wolff?**
    — Porque para L=128 los algoritmos locales no han decorrelacionado suficientemente entre muestras, aunque se tome una cada 50 pasos. Las muestras siguen correlacionadas, lo que introduce un sesgo sistemático en la estimación de χ y desplaza el pico del colapso.

14. **¿Por qué no habéis simulado L = 256?**
    — El coste computacional escala como L^{2+z}: para Metropolis/Glauber eso es L^4, multiplicar L por 2 supone un coste ×16. Con Wolff el coste es mucho más razonable, pero los datos para L=128 ya muestran convergencia clara de los exponentes.

15. **La extrapolación de T_c desde los cruces del Binder da resultados ruidosos para Metropolis/Glauber. ¿Por qué no usáis directamente el máximo de χ?**
    — El máximo de χ depende de L y hay que extrapolar de todas formas. El Binder es teóricamente más limpio porque el cruce converge a T_c desde ambos lados monótonamente; el problema aquí no es el método sino la calidad de los datos de entrada, que el critical slowing down degrada para L grandes.

16. **¿Cuál es la diferencia física entre la longitud de correlación ξ_L y el ratio ξ_L/L?**
    — ξ_L diverge en el límite termodinámico; ξ_L/L es un ratio adimensional que, como U_4, tiene la propiedad de cruzarse en T_c independientemente de L, proporcionando una estimación alternativa e independiente de T_c.

17. **¿Habéis comprobado la invarianza bajo el grupo de renormalización explícitamente?**
    — De forma implícita sí: el colapso FSS es precisamente la expresión numérica de que la física cerca de T_c es invariante de escala. Los exponentes que medimos son los que parametrizan esa invarianza.

18. **¿Qué haríais diferente si repetiráis el trabajo?**
    — Paralelizar las simulaciones por temperatura (trivialmente paralelizable), añadir L=256 al menos con Wolff, y aplicar análisis de errores más sofisticado como MCMC con correcciones de autocorrelación (método de Madras-Sokal o método de Γ).
