// Alumno: Álvaro Rodríguez Gallardo
// DNI: 77034155W
// Correo: alvaro155w@correo.ugr.es
// Grupo 2 (Miércoles, 17.30-19.30)
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>
#include "aux.h"
#include "random.hpp"

using Random = effolkronium::random_static;

//*******************************************************************************************************************************************************
//************************************************************** CONSTANTES P1 ***************************************************************************
//*******************************************************************************************************************************************************

// Varianza de la distribución normal
const double VARIANZA = 0.3;

// Criterio de parada por iteraciones
const int MAXIMO_ITERACIONES = 15000;

// Criterio de parada por límite
const int LIMITE = 20;


//*******************************************************************************************************************************************************
//************************************************************** ESTADÍSTICOS ***************************************************************************
//*******************************************************************************************************************************************************

/**
    @fn tasaClasificacion
    @brief Dadas las particiones de un conjunto de datos, y los pesos, se calcula la tasa de clasificación
    @param entrenamiento. Dataset que sirve para entrenar al algoritmo
    @param test. Dataset usado para el test
    @param pesos. Vector de pesos para el cálculo de la tasa
    @return valor de la tasa de clasificación
*/
double tasaClasificacion(const Dataset & entrenamiento, const Dataset & test, const Pesos & pesos);

/**
    @fn tasaClasificacionLeaveOneOut
    @brief Calcula la tasa de aciertos en clasificación usando leave one out
    @param entrenamiento. Dataset de entrenamiento
    @param pesos. Vector de pesos
    @return tasa de clasificación usando leave one out
*/

double tasaClasificacionLeaveOneOut(const Dataset & entrenamiento, const Pesos & pesos);

/**
    @fn tasaReduccion
    @brief Dados los pesos en un momento del problema, se calcula la tasa de reducción
    @param pesos. Pesos para el cálculo de la tasa
    @return valor de la tasa de reducción
*/
double tasaReduccion(const Pesos & pesos);

/**
    @fn fitness
    @brief Función objetivo que se usa para buscar el máximo
    @param tasa_cla. Tasa de clasificación
    @param tasa_red Tasa de reducción
    @return valor de la función fitness
*/

double fitness(double tasa_cla, double tasa_red);


//*******************************************************************************************************************************************************
//****************************************************************** 1-NN *******************************************************************************
//*******************************************************************************************************************************************************

/**
    @fn algoritmo1NN
    @brief Dada una muestra, la va a clasificar, teniendo en cuenta los pesos
    @param m. Muestra tomada que se va a clasificar
    @param entrenamiento. Dataset para usar en el algoritmo
    @param pesos. Pesos
    @cond m.caracteristicas.size() == entrenamiento[j].muestras[i].caracteristicas.size()==pesos.valores.size() para todo i, para todo j
    @return clase de m
*/
string algoritmo1NN(const Muestra & m, const Dataset & entrenamiento, const Pesos & pesos);

//*******************************************************************************************************************************************************
//************************************************************** BÚSQUEDA LOCAL (BL) ********************************************************************
//*******************************************************************************************************************************************************

/**
    @fn BusquedaLocal
    @brief Aplica búsqueda local con esquema El Primero Mejor para calcular los pesos
    @param entrenamiento. Dataset para el algoritmo
    @cond entrenamiento[j].muestras[i].caracteristicas.size() == pesos.valores.size() para todo i, para todo j
    @return Vector de pesos y fitness asociado
*/
std::pair<Pesos,double> BusquedaLocal(const Dataset & entrenamiento);

/**
    @fn mutacionBL
    @brief Dado un vector de pesos, se realiza la mutación de una característica
    @param pesos. Vector que va a mutar
    @param index. Posición del elemento a mutar
    @param normal. Distribución normal para mutar
*/

void mutacionBL(Pesos & pesos, int index, std::normal_distribution<double> normal);

/**
    @fn mezclarIndices
    @brief Generar el vector de índices y mezclar el vector aleatoriamente
    @param n_pesos. Cantidad de pesos
    @param indexes. Vector de índices que se quiere mezclar
    @param inicializar. Recorrer o no el bucle para inicializar el vector de índices

*/

void mezclarIndices(int n_pesos,vector<int> & indexes, bool inicializar);

//*******************************************************************************************************************************************************
//******************************************************************* RELIEF ****************************************************************************
//*******************************************************************************************************************************************************

/**
    @fn GreedyRelief
    @brief Algoritmo greedy relief para APC
    @param entrenamiento. Dataset para el algoritmo
    @return Vector de pesos y fitness asociado
*/

std::pair<Pesos,double> GreedyRelief(const Dataset & entrenamiento);

/**
    @fn buscarEnemigoMasCercano
    @brief Dada una muestra, busca su enemigo más cercano (muestra en entrenamiento con menor distancia a m con clase diferente)
    @param entrenamiento. Dataset de entrenamiento del que sacar el enemigo
    @param m. Muestra a la que se le busca un enemigo más cercano
    @return enemigo más cercano
*/

Muestra buscarEnemigoMasCercano(const Dataset & entrenamiento, const Muestra & m);

/**
    @fn buscarAmigoMasCercano
    @brief Dada una muestra, busca su amigo más cercano (muestra con misma clase y mínima distancia, distinto a la muestra)
    @param entrenamiento. Dataset de entrenamiento en el que se busca el amigo
    @param m. Muestra a la que se le busca un amigo más cercano
    @return amigo más cercano
*/

Muestra buscarAmigoMasCercano(const Dataset & entrenamiento, const Muestra & m);

//*******************************************************************************************************************************************************
//******************************************************************* RESULTADOS ************************************************************************
//*******************************************************************************************************************************************************

/**
    @fn resultadosP1
    @brief Muestra los resultados de la práctica 1 en función del parámetro
    @param caso. indica cómo se obtiene el vector de pesos
*/
void resultadosP1(string caso);