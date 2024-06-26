// Alumno: Álvaro Rodríguez Gallardo
// DNI: 77034155W
// Correo: alvaro155w@correo.ugr.es
// Grupo 2 (Miércoles, 17.30-19.30)

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>
#include <utility>
#include <iomanip> 

using namespace std;


//*******************************************************************************************************************************************************
//************************************************************** CONSTANTES *****************************************************************************
//*******************************************************************************************************************************************************

// Constante que pondera la importancia entre acierto y reducción de la solución encontrada
const float ALPHA = 0.75;

// Número de datasets: 'breast_cancer', 'ecoli', 'parkinsons'
const int NUM_CONJUNTOS_DATOS = 3;

// Para cada dataset, el número de archivos en que se particiona
const int NUM_PARTICIONES = 5;

// Nombres de los datasets y extensión, para facilitar la lectura
const string BCANCER = "breast-cancer_";
const string ECOLI = "ecoli_";
const string PARKINSON = "parkinsons_";
const string EXTENSION = ".arff";

//*******************************************************************************************************************************************************
//************************************************************* REPRESENTACIÓN SOLUCIÓN *****************************************************************
//*******************************************************************************************************************************************************

/**
    @struct Atributo
    @brief Representación del atributo: el nombre y el tipo de dato
*/
struct Atributo {
    string nombre;  // Nombre de un atributo, p.e., 'radius_1'
    string tipo;    // Tipo del atributo, p.e., 'real'
};

/**
    @struct Solucion
    @brief Representación de pesos y clasificación de una muestra
*/
struct Muestra {
    vector<double> caracteristicas; // Valores de las características de una muestra
    string clase;                   // Clasificación de una muestra
};

/**
    @struct Dataset
    @brief Representa un dataset dado por un archivo .arff
*/
struct Dataset {
    vector<Atributo> atributos;     // Atributos que hay en un archivo .arff
    vector<Muestra> muestras;       // Todas las muestras que hay en un archivo .arff
};

/**
    @struct Pesos
    @brief Representa los pesos del problema que se afronta
*/
struct Pesos {
    vector<double> valores;
};


//*******************************************************************************************************************************************************
//******************************************************* PREPROCESAMIENTO DE DATOS *********************************************************************
//*******************************************************************************************************************************************************

/**
    @fn lecturaFichero
    @brief Función que lee un fichero de extensión .arff
    @param fichero. Nombre del archivo que se lee
    @return Dataset. Representación de los datos del fichero para ser usado en el programa
*/
Dataset lecturaFichero(const string fichero);

/**
    @fn normalizarDatos
    @brief Dada una lista de dataset, se normalizan estos
    @param datasets. Vector de datasets que se van a normalizar
*/
void normalizarDatos(vector<Dataset> & vect_ds);

/**
    @fn unirDatasets
    @brief Dado un vector de datasets (posibles particiones), se unifican en un único dataset
    @param vect_ds. Vector de datasets a unir
    @return dataset unificado
*/
Dataset unirDatasets(const vector<Dataset>& vect_ds);

//*******************************************************************************************************************************************************
//************************************************************* DISTANCIAS ******************************************************************************
//*******************************************************************************************************************************************************

/**
    @fn distanciaEuclidea
    @brief Dadas dos muestras de un dataset, se calcula la distancia (euclídea) entre ellos
    @param m1. Una muestra del dataset
    @param m2. Una muestra del dataset
    @return Distancia entre m1 y m2
*/
double distanciaEuclidea(const Muestra & m1, const Muestra & m2);

/**
    @fn distanciaEuclideaPonderada
    @brief Dadas dos muestras de un dataset, se calcula la distancia (euclídea) ponderada entre ellos
    @param m1. Una muestra del dataset
    @param m2. Una muestra del dataset
    @param pesos. Pesos por el que se pondera la distancia
    @return Distancia ponderada entre m1 y m2
*/
double distanciaEuclideaPonderada(const Muestra & m1, const Muestra & m2, const Pesos & pesos);

/**
    @fn buscarMuestraAMenorDistancia
    @brief Dado un dataset, busca la muestra con menor distancia a la dada
    @param ds. Dataset en el que buscar
    @param m. Muestra a la que se busca la otra muestra más cercana
    @return muestra más cercana a m
*/
int buscarMuestraAMenorDistancia(const Dataset & ds, const Muestra & m);

//*******************************************************************************************************************************************************
//******************************************************************* AUXILIARES ************************************************************************
//*******************************************************************************************************************************************************

/**
    @fn mostrarDataset
    @brief Muestra las muestras de un dataset. Posible uso en depuración
    @param ds. Dataset a mostrar
*/
void mostrarDataset(const Dataset & ds);

/**
    @fn escribirCSV
    @brief Escribir datos en un .csv con separador ; para escribir resultados más rápido en Excel
    @param tasa_cla tasa de clasificación
    @param tasa_red tasa de reducción
    @param fitness fitness
    @param tiempo tiempo de ejecución
    @param media se presenta la media de valores o no
*/
void escribirCSV(double tasa_cla,double tasa_red,double fitness,double tiempo,bool media);

/**
    @fn escribirFitnessParaConvergencia
    @brief Escribir el fitness en archivos para estudiar convergencia de algoritmos
    @param num_part número de la partición
    @param fitness valor a escribir
    @param n_generacion número de la generación
    @param algoritmo
    @param nombre_dataset
*/
void escribirFitnessParaConvergencia(int num_part, double fitness, int n_generacion, string algoritmo, string nom_carpeta);