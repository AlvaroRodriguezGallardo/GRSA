#PROJECT(CXX)
#cmake_minimum_required(VERSION 3.16)
#ADD_EXECUTABLE(random SRC/ejemplo_random.cc)
#TARGET_LINK_LIBRARIES(random)
#SET(TIEMPOS_SRC SRC/fact.cc SRC/tiempos.cc)
#ADD_EXECUTABLE(tiempos ${TIEMPOS_SRC})
#TARGET_LINK_LIBRARIES(tiempos m)

cmake_minimum_required(VERSION 3.16)

# Establecer la variable de compilación CMAKE_CXX_FLAGS con la opción de optimización -O2
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Nombre del proyecto
project(practicaAlternativa)

# Archivos de cabecera
include_directories(INCLUDE)

# Lista de archivos fuente
set(SOURCE_FILES
	SRC/aux.cpp
	SRC/cec17.c
	SRC/cec17_test_func.c
	SRC/practica1.cpp
	SRC/practicaAlternativa.cpp
        SRC/main.cpp)

# Ejecutable
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

# Carpeta de destino para el ejecutable
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../BIN
)
