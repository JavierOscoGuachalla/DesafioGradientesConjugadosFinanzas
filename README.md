# Simulador de Ajustes Presupuestarios mediante Gradientes Conjugados Precondicionados

**Autor: Javier Osco Guachalla**
**Mencion: Desarrollo de Software**

## 1. Descripción del proyecto

Aplicación web desarrollada como proyecto académico de la asignatura **Métodos Numéricos**, cuyo objetivo es simular la distribución y el ajuste de un presupuesto entre diferentes áreas de una empresa de desarrollo de software.

El simulador utiliza el **Método de Gradientes Conjugados Precondicionado (GCP)** para resolver sistemas de ecuaciones lineales y encontrar una asignación presupuestaria que se aproxime a las necesidades y prioridades de cada área.

El proyecto permite modificar los parámetros financieros, ejecutar diferentes escenarios y comparar los resultados del método precondicionado con el método de Gradientes Conjugados convencional.

## 2. Objetivos

- Aplicar métodos numéricos a un problema de distribución presupuestaria.
- Formular un sistema de ecuaciones lineales \(Ax=b\).
- Implementar el método de Gradientes Conjugados Precondicionado.
- Comparar su rendimiento con el método de Gradientes Conjugados convencional.
- Evaluar la precisión de los resultados y el condicionamiento de las matrices.
- Permitir la modificación de parámetros y la simulación de distintos escenarios.

## 3. Modelo matemático

El problema se representa mediante un sistema de ecuaciones lineales:

\[
Ax=b
\]

Donde:

- **A:** matriz de coeficientes que representa las relaciones y restricciones del modelo financiero.
- **x:** vector de variables que representa los presupuestos asignados a cada área.
- **b:** vector de términos independientes que contiene los valores de referencia del modelo.

La aplicación utiliza seis áreas presupuestarias como variables del sistema.

La matriz debe ser simétrica y definida positiva (SPD), condiciones necesarias para aplicar correctamente el método de Gradientes Conjugados.

## 4. Métodos numéricos implementados

La aplicación incorpora los siguientes procedimientos:

- **Gradientes Conjugados (GC):** método iterativo utilizado como referencia de comparación.
- **Gradientes Conjugados Precondicionados (GCP):** método principal para mejorar la convergencia del sistema.
- **Precondicionadores:** Jacobi, IC(0) y SSOR.
- **Métodos directos:** Cholesky, eliminación de Gauss y factorización LU.

El criterio de convergencia utilizado es:

\[
\frac{\|r_k\|_2}{\|b\|_2}<10^{-6}
\]

Donde \(r_k\) representa el residuo en la iteración actual.

## 5. Escenarios de simulación

El simulador permite evaluar tres situaciones financieras:

1. **Escenario ideal:** representa condiciones presupuestarias estables.
2. **Escenario de estrés financiero:** representa modificaciones en las prioridades y disponibilidad de recursos.
3. **Escenario mal condicionado:** permite analizar cómo afecta el condicionamiento de la matriz a la convergencia de los métodos numéricos.

## 6. Funcionalidades principales

- Visualización y edición de la matriz A y el vector b.
- Configuración del vector inicial \(x_0\).
- Selección del método de resolución directa.
- Selección del precondicionador.
- Configuración de tolerancia y número máximo de iteraciones.
- Validación de las condiciones de simetría y definición positiva.
- Ejecución de los tres escenarios financieros.
- Comparación entre GC y GCP.
- Cálculo de errores absolutos, relativos y porcentuales.
- Evaluación del residuo y número de condición.
- Pruebas de validación del sistema.

## 7. Tecnologías utilizadas

- **HTML5:** estructura de la aplicación.
- **CSS3:** diseño y presentación de la interfaz.
- **JavaScript:** implementación de los algoritmos numéricos, validaciones y simulaciones.

No requiere instalación de dependencias ni configuración de un servidor para su ejecución local.

## 8. Estructura del proyecto

```text
simulador-presupuesto-gradientes-conjugados/
│
├── index.html
├── styles.css
├── script.js
└── README.md
```

- `index.html`: estructura y contenido de la aplicación.
- `styles.css`: estilos visuales y diseño de la interfaz.
- `script.js`: modelo matemático, algoritmos numéricos, validaciones y cálculos.
- `README.md`: documentación del proyecto.

## 9. Ejecución

Para ejecutar la aplicación localmente:

1. Descargar o clonar este repositorio.
2. Abrir la carpeta del proyecto.
3. Ejecutar el archivo `index.html` en un navegador web moderno.

También es posible acceder a la versión pública mediante el enlace de despliegue del proyecto.

## 10. Alcance académico

Este simulador fue desarrollado con fines educativos para demostrar la aplicación de métodos numéricos en un problema financiero.

Los valores y escenarios utilizados son datos ilustrativos y no representan registros contables reales de una empresa.

El proyecto corresponde al desafío académico de la asignatura **Métodos Numéricos**, en la mención de Desarrollo de Software.

## 11. Información académica

**Autor:** Javier Osco Guachalla  
**Universidad:** Universidad Mayor de San Andrés (UMSA)  
**Carrera:** Informática  
**Mención:** Desarrollo de Software  
**Asignatura:** Métodos Numéricos
