## Descripción de los archivos

De la carpeta 'software' cuelgan los archivos que implementan la práctica, aquellos que contienen los datos,... A continuación se explica cada elemento:

- Carpeta 'BIN':
    - Carpeta 'DATA': Contiene los ficheros con los datos que se van a tratar.
    - 'practica3': Ejecutable de la práctica. Se debe indicar la semilla con la que se ejecuta el programa.
    - 'LEEME.md': Fichero estilo MarkDown que describe los elementos que cuelgan de la carpeta 'software'.

- Carpeta 'FUENTES':
    - Carpeta 'INCLUDE':
        - 'INCLUDE/aux.h': Declaración de cabeceras de funciones, estructuras de datos y constantes que se usarán en el programa. Funciones como la lectura de ficheros, normalización, distancias,..., se encuentran aquí.
        - 'INCLUDE/practica1.hpp': Declaración de las funciones asociadas a la práctica 1. Las funciones del algoritmo 1-NN, búsqueda local y greedy RELIEF se encuentran aquí, junto a funciones para mostrar los resultados.
        - 'INCLUDE/practica2.hpp': Declaración de las funciones asociadas a la práctica 2, tales como los algoritmos genéticos generacionales, algoritmos genéticos estacionarios y algoritmos meméticos.
        - 'INCLUDE/practica3.hpp': Declaración de las funciones asociadas a la práctica 3, tales como BMB, ILS, ES e ILS-ES, junto a las correspondientes funciones auxiliares.
        - 'INCLUDE/practicaAlternativa.hpp': Declaración de las funciones asociadas a GRSA y sus mejoras.
        - 'INCLUDE/random.hpp': Fichero ofrecido por el profesor para tratar con aleatoriedad de forma cómoda.

    - Carpeta 'SRC': 
        - 'SRC/aux.cpp': Implementación de las funciones declaradas en 'INCLUDE/aux.h'.
        - 'SRC/practica1.cpp': Implementación de las funciones declaradas en 'INCLUDE/practica1.hpp'.
        - 'SRC/practica2.cpp': Implementación de las funciones declaradas en 'INCLUDE/practica2.hpp'.
        - 'SRC/practica3.cpp': Implementación de las funciones declaradas en 'INCLUDE/practica3.hpp'.
        - 'SRC/practicaAlternativa.cpp': Implementación de las funciones declaradas en 'INCLUDE/practicaAlternativa.hpp'.
        - 'SRC/main.cpp': Función principal del programa, recibe los parámetros y llama a las funciones necesarias para calcular y mostrar resultados.

    - Carpeta 'CMakeFiles': Carpeta generada a partir del archivo de generación de ficheros para la compilación.
    - Fichero 'CMakeLists.txt': Fichero con las órdenes necesarias para realizar los enlaces de bibliotecas y generar el archivo 'Makefile' para la compilación.
    - Fichero 'Makefile': Usado para la generación del ejecutable del programa.
    - Fichero 'CMakeCache.txt': Generado por 'cmake .' para obtener el archivo 'Makefile'.
    - Fichero 'cmake_install.cmake': Generado por 'cmake .' para obtener el archivo 'Makefile'.



## Obtención del ejecutable y ejecución

Seguir los siguientes pasos:

- **Paso 1**: Abrir la terminal, o ir, a la carpeta 'software/FUENTES'.

- **Paso 2**: Ejecutar el siguiente comando.
```
cmake .
```

- **Paso 3**: Ejecutar el comando siguiente, que creará 'practicaAlternativa' en la carpeta 'software/BIN'.
```
make
```

- **Paso 4**: Ir a la carpeta 'software/BIN'.

- **Paso 5**: Ejecutar el siguiente comando, siendo ``<semilla>`` un número que indica la semilla que se usará.
```
./practicaAlternativa <semilla>
```
