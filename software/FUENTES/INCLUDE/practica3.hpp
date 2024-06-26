// Alumno: Álvaro Rodríguez Gallardo
// DNI: 77034155W
// Correo: alvaro155w@correo.ugr.es
// Grupo 2 (Miércoles, 17.30-19.30)

#include "practica2.hpp"
#include <unordered_set>

//************************************************************************ CONSTANTES ***************************************************************************

//************************************************************************ CONSTANTES ES ***************************************************************************


const double PROBAB_ACEPTAR_SOLUC_MU_POR_1_PEOR = 0.3;

const double MU = 0.1;

const int LIMITE_MAX_VECINOS = 10;

const double LIMITE_MAX_EXITOS = 0.1;

const int NUM_ITERACIONES_ES = 15000;

//************************************************************************ CONSTANTES BMB ***************************************************************************

const int N_SOLS_INICIALES_BMB_ILS = 20;

//************************************************************************ CONSTANTES ILS ***************************************************************************

const double LIMITE_MUTACION_ILS = 0.2;

const int N_CARACTS_DEBEN_CAMBIAR = 3;

const double LIM_INF_DISTRIBUCION = -0.25;

const double LIM_SUP_DISTRIBUCION = 0.25;

//************************************************************************ CONSTANTES BMB e ILS ***************************************************************************

const int MAX_EVALS_BL_EN_BMB_ILS = 750;

//************************************************************************ ESTRUCTURAS DE DATOS ***************************************************************************

struct Solucion{
    Pesos pesos;    // Vector de pesos
    double fitness; // Valor de la función objetivo
};

//************************************************************************ FUNCIONES AUXILIARES GENERAL ***************************************************************************

/**
    @fn BL_dada_sol_inicial
    @brief Implementa búsqueda local donde la solución inicial es dada
    @param entrenamiento. Conjunto de entrenamiento
    @param sol_inicial. Solución inicial para BL
    @return Solución encontrada por BL
*/
std::pair<Pesos,double> BL_dada_sol_inicial(const Dataset & entrenamiento, const Pesos &  sol_inicial);

/**
    @fn generarSolucionAleatoria
    @brief Genera una solución aleatoria según una uniforme en [0,1]
    @param solucion. Solución a inicializar
    @param n_pesos
*/
void generarSolucionAleatoria(Pesos & solucion, int n_pesos);

//************************************************************************ FUNCIONES AUXILIARES BMB ***************************************************************************

/**
    @fn compararFitness
    @brief compara el fitness de dos soluciones para ordenar de forma descendente
    @param s1. Solución 1
    @param s2. Solución 2
    @return True si s1 es mejor, false en otro caso
*/
bool compararFitness(const Solucion& s1, const Solucion& s2);

//************************************************************************ FUNCIONES AUXILIARES ILS ***************************************************************************

/**
    @fn mutacionFuerteBL
    @brief Aplica mutación sobre un peso según una uniforme (-0.25,0.25)
    @param pesos
    @param index
    @param uniforme
*/
void mutacionFuerte(Pesos & pesos, int index, std::uniform_real_distribution<double> uniforme);

//************************************************************************ FUNCIONES AUXILIARES ES ***************************************************************************

/**
    @fn obtenerTemperaturaInicial
    @brief Obtención de la temperatura incial según el esquema de Cauchy modificado
    @param coste_sol_inicial. Valor fitness de la solución inicial
    @return valor inicial de la temperatura
*/
double obtenerTemperaturaInicial(double coste_sol_inicial);

/**
    @fn obtenerTemperaturaCauchy
    @brief Dada la temperatura T_k, se obtiene T_{k+1} según el esquema de Cauchy modificado
    @param t_k
    @param beta
    @return valor de T_{k+1}
*/
double obtenerTemperaturaCauchy(double t_k, double beta);

/**
    @fn obtenerBeta
    @brief Obtención del parámetro beta
    @param t_inicial
    @param max_vecinos
    @param t_final
    @return beta
*/
double obtenerBeta(double t_inicial,int max_vecinos, double t_final);

/**
    @fn limiteExponencial
    @brief Aplica la fórmula de la exponencial para el límite de aceptación
    @param delta
    @param T
    @return valor de la función
*/
double limiteExponencial(double delta, double T);

//************************************************************************ ALGORITMOS ***************************************************************************

/**
    @fn BMB
    @brief Implementa el algoritmo de multiarranque básico
    @param entrenamiento. Dataset de entrenamiento
    @return solución y fitness
*/
std::pair<Pesos,double> BMB(const Dataset & entrenamiento);

/**
    @fn ES
    @brief Implementa el algoritmo de enfriamiento simulado
    @param entrenamiento. Dataset de entrenamiento
    @return solución y fitness
*/
std::pair<Pesos,double> ES(const Dataset & entrenamiento);

/**
    @fn ES_dada_sol_inicial
    @brief Aplica enfriamiento simulado a una solución dada
    @param entrenamiento. Dataset de entrenamiento
    @param sol_inicial Solución inicial
    @return solución de enfriamiento simulado y su fitness
*/
std::pair<Pesos,double> ES_dada_sol_inicial(const Dataset& entrenamiento, const Pesos& sol_inicial);

/**
    @fn ILS
    @brief Implementa el algoritmo iterativo de búsqueda local de la P1
    @param entrenamiento. Dataset de entrenamiento
    @return solución y fitness
*/
std::pair<Pesos,double> ILS(const Dataset & entrenamiento);

/**
    @fn BMB
    @brief Implementa el ILS usando ES
    @param entrenamiento. Dataset de entrenamiento
    @return solución y fitness
*/
std::pair<Pesos,double> ILS_ES(const Dataset & entrenamiento);

//************************************************************************ RESULTADOS ***************************************************************************

/**
    @fn resultadosP3
    @brief expone los resultados de los algoritmos de la práctica 3
    @param caso
*/
void resultadosP3(const string & caso);