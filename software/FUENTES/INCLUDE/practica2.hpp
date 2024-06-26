// Alumno: Álvaro Rodríguez Gallardo
// DNI: 77034155W
// Correo: alvaro155w@correo.ugr.es
// Grupo 2 (Miércoles, 17.30-19.30)
//#include <iomanip>
//#include <string>
//#include <cmath>
//#include <algorithm>
//#include "aux.h"
#include "practica1.hpp"
//#include "random.hpp"

//****************************************************** Constantes Problemas *************************************************************************************

const int N_CROMOSOMAS=50;

const int MAX_EVALUACIONES = 15000;

//const double VARIANZA = 0.3;

//****************************************************** Constantes Genéticos *************************************************************************************

const double ALPHA_BLX = 0.3;

const double PROB_CRUCE_AGG = 0.7;

const double PROB_CRUCE_AGE = 1;

const double PROB_MUTAR_INDIVIDUO_GENETICOS = 0.08;

//****************************************************** Constantes Meméticos *************************************************************************************

const int N_GENERACIONES = 10;

const double PROB_CRUCE_MEMETICOS = 0.7;

const double PROB_MUTAR_INDIVIDUO_MEMETICOS = 0.08;

const double PROB_MEMETICOS_BL = 0.1;

const int LIMITE_BL = 2;

//****************************************************** EDs AUXILIARES *************************************************************************************

class Cromosoma{
    private:
        Pesos pesos;
        double valor;
    public:
        // Constructor con parámetros
        Cromosoma(const Dataset& entrenamiento, int N_pesos);

        // Constructor con parámetros
        Cromosoma(const Pesos& pesos, double valor);

        /**
            @fn esMejorQue
            @brief Compara si el cromosoma tiene mejor valor que otro
            @param crom. Cromosoma con el que se compara
            @return True si es mejor, False en otro caso
        */
        bool esMejorQue(const Cromosoma& crom) const;

        /**
            @fn getPesos
            @brief obtener el vector de pesos
            @return pesos, el vector de pesos
        */
        Pesos getPesos() const;

        /**
            @fn getValor
            @brief obtener valor de fitness
            @return valor, el valor de fitness
        */
        double getValor() const;

        /**
            @fn mutacion
            @brief Realiza la mutación del vector de pesos
            @param indice. Posición que se va a mutar
            @param entrenamiento. Dataset de entrenamiento para volver a calcular fitness
        */
        void mutacion(int indice,const Dataset& entrenamiento);
};

struct Poblacion {
    vector<Cromosoma> cromosomas;
};

//****************************************************** FUNCIONES AUXILIARES *************************************************************************************

/**
    @fn inicializacion
    @brief Inicializar la población con cromosomas según datos de entrenamiento
    @param entrenamiento. Datos para inicializar la población
*/
void inicializacion(const Dataset& entrenamiento);

//****************************************************** ALGORITMOS GENÉTICOS *************************************************************************************

/**
    @fn Cruce_BLX
    @brief Hace la mutación BLX-alpha de dos cromosomas
    @param c1. Cromosoma 1
    @param c2. Cromosoma 2
    @return Vector de dos cromosomas mutados
*/
vector<Cromosoma> Cruce_BLX(const Cromosoma& c1, const Cromosoma& c2,const Dataset& entrenamiento);

/**
    @fn Cruce_Arit
    @brief Hace la mutación aritmética aleatoria de dos cromosomas
    @param c1. Cromosoma 1
    @param c2. Cromosoma 2
    @param entrenamiento
    @return Dos cromosomas mutados aleatoriamente
*/
vector<Cromosoma> Cruce_Arit(const Cromosoma& c1, const Cromosoma& c2, const Dataset& entrenamiento);

/**
    @fn torneo
    @brief Operador de selección de una población. Implementa un torneo en que se cogen 3 cromosomas aleatorios y se coge el mejor
    @param poblacion. Población de cromosomas
    @return Cromosoma que gana al resto
*/
Cromosoma torneo(const Poblacion& poblacion);

/**
    @fn reemplazoGeneracional
    @brief Esquema de reemplazo para AGG. Para conservar elitismo, SI la mejor solución de población no sobrevive, esta sustituye a la peor de la nueva población, A MENOS QUE se encuentre en la nueva población  una solución mejor
    @param poblacion. Población en el instante n
    @param hijos. Población en el instante n+1
    @return Población tras el reemplazo elitista
*/
Poblacion reemplazoGeneracional(const Poblacion& poblacion, const Poblacion& hijos);

/**
    @fn reemplazoEstacionario
    @brief Esquema de reemplazo para AGE. Los hijos sustituyen (SI SON MEJORES) a las dos peores
    @param poblacion. Población de cromosomas
    @param hijos. Hijos de dos cromosomas de la población
*/
void reemplazoEstacionario(Poblacion& poblacion, const vector<Cromosoma>& hijos);

/**
    @fn estaDentro
    @brief Comprueba si un cromosoma está dentro de una población
    @param poblacion
    @param crom
    @return True si está dentro, False si no
*/
bool estaDentro(const Poblacion & poblacion, const Cromosoma& crom);

/**
    @fn mejorCromosoma
    @brief De una población, escoge el mejor cromosoma
    @param poblacion
    @return el mejora cromosoma de la población
*/
std::pair<Cromosoma,int> mejorCromosoma(const Poblacion& poblacion);

/**
    @fn peorCromosoma
    @brief En una población, devuelve el peor cromosoma
    @param poblacion
    @return peor cromosoma de la población
*/
std::pair<Cromosoma,int> peorCromosoma(const Poblacion& poblacion);
/**
    @fn AGG_BLX
    @brief Aplicar AGG-BLX para obtener vector de pesos tras entrenar
    @param entrenamiento. Dataset con las muestras para entrenar el modelo
    @param nombre del dataset
    @param n_part número de partición
    @return vector de pesos entrenados
*/
std::pair<Pesos,double> AGG_BLX(const Dataset& entrenamiento, string nombre, int n_part);

/**
    @fn AGG_Arit
    @brief Aplicar AGG-CA para obtener vector de pesos tras entrenar
    @param entrenamiento. Dataset con las muestras para entrenar el modelo
    @param nombre del dataset
    @param n_part número de partición
    @return vector de pesos entrenados
*/
std::pair<Pesos,double> AGG_Arit(const Dataset& entrenamiento, string nombre, int n_part);

/**
    @fn AGE_BLX
    @brief Aplicar AGE-BLX para obtener vector de pesos tras entrenar
    @param entrenamiento. Dataset con las muestras para entrenar el modelo
    @param nombre del dataset
    @param n_part número de partición
    @return vector de pesos entrenados
*/
std::pair<Pesos,double> AGE_BLX(const Dataset& entrenamiento, string nombre, int n_part);

/**
    @fn AGE_Arit
    @brief Aplicar AGE-CA para obtener vector de pesos tras entrenar
    @param entrenamiento. Dataset con las muestras para entrenar el modelo
    @param nombre del dataset
    @param n_part número de partición
    @return vector de pesos entrenados
*/
std::pair<Pesos,double> AGE_Arit(const Dataset& entrenamiento, string nombre, int n_part);

//****************************************************** ALGORITMOS MEMÉTICOS *************************************************************************************

// NOTA: Se junta el mejor resultado de BL (práctica 1) con el mejor AGG (de los dos que hay)

/**
    @fn devolverSubconjuntoCromosomas
    @brief Dada una población, se devuelve el subconjunto de tamaño tam_subconjunto. Si elMejor=false, devuelve uno según aleatoriedad. En otro caso, devuelve los mejores según fitness
    @param poblacion
    @param tam_subconjunto. Tamaño del subconjunto a devolver
    @param esMejor. Usar aleatoriedad o según fitness
    @return subconjunto de cromosomas pos sus índices
*/
vector<int> devolverSubconjuntoCromosomas(const Poblacion& poblacion, int tam_subconjunto=-1, bool elMejor=false);

/**
    @fn busquedaLocalPractica2
    @brief Realiza, para el vector de pesos de un cromosoma, búsqueda local sobre el mismo bajo ciertas condiciones
    @param entrenamiento. Para las evaluaciones de fitness
    @param crom. Cromosoma al que cambiar vector de pesos
    @return número de evaluaciones de fitness
*/
int busquedaLocalPractica2(const Dataset& entrenamiento, Cromosoma& crom);

// TRAS ALGUNAS EJECUCIONES, SE IMPLEMENTA AGG-CA AL DAR LIGEROS RESULTADOS MEJORES QUE AGG-BX

/**
    @fn AM_All
    @brief Aplicar AM-(10,1) para obtener vector de pesos tras entrenar
    @param entrenamiento. Dataset con las muestras para entrenar el modelo
    @param nombre del dataset
    @param n_part número de partición
    @return vector de pesos entrenados
*/
std::pair<Pesos,double> AM_All(const Dataset& entrenamiento,string nombre, int n_part);

/**
    @fn AM_Rand
    @brief Aplicar AM-(10,0.1) para obtener vector de pesos tras entrenar
    @param entrenamiento. Dataset con las muestras para entrenar el modelo
    @param nombre del dataset
    @param n_part número de partición
    @return vector de pesos entrenados
*/
std::pair<Pesos,double> AM_Rand(const Dataset& entrenamiento, string nombre, int n_part);

/**
    @fn AM_Best
    @brief Aplicar AM-(10,0.1mej) para obtener vector de pesos tras entrenar
    @param entrenamiento. Dataset con las muestras para entrenar el modelo
    @param nombre del dataset
    @param n_part número de partición
    @return vector de pesos entrenados
*/
std::pair<Pesos,double> AM_Best(const Dataset& entrenamiento, string nombre, int n_part);


//************************************************************* RESULTADOS *************************************************************************************

/**
    @fn resultadosP2
    @brief Ejecutar y mostrar resultados de un algoritmo
    @param caso. Algoritmo que se quiere ejecutar
*/
void resultadosP2(const string& caso);