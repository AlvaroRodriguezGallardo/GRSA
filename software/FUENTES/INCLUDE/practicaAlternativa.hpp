// Alumno: Álvaro Rodríguez Gallardo
// DNI: 77034155W
// Correo: alvaro155w@correo.ugr.es
// Grupo 2 (Miércoles, 17.30-19.30)

#include "practica1.hpp"
extern "C" {
#include "cec17.h"
}

//******************************************************************************************************************************************************************
//******************************************************************** CONSTANTS ***********************************************************************************
//******************************************************************************************************************************************************************

// Some constatns depend on input data
// This dependencies are because of a good performance of GRSA

// Low limit of geometric coefficients
const double GM_MIN = 0;

// High limit of geometric coefficients
const double GM_MAX = 1;

// Geometric coeffient 1
const double GM_1 = 0.01;

// Geometric coefficient 2
const double GM_2 = 0.99;

// Max evaluations of Local Search for Memetic Algorithm
const int MAX_EVALS_BL = 750;

//******************************************************************************************************************************************************************
//********************************************************************** DATA STRUCTS ******************************************************************************
//******************************************************************************************************************************************************************

struct Particula{
    Pesos posicion;         // Position in space-time
    Pesos mejor_posicion;   // Best position particle has got, for a posible mutation
    double fitness;         // Fitness of the particle

    bool operator==(const Particula& other) const {
        return (posicion.valores == other.posicion.valores && fitness == other.fitness);
    }
};

struct EspacioTiempo{
    vector<Particula> T;    // Set of particles, named 'Tensor'
};

//******************************************************************************************************************************************************************
//************************************************************** FUNCTIONS OF GRSA**********************************************************************************
//******************************************************************************************************************************************************************

/**
    @fn inicializarParametros
    @brief initialize dependent parameters like T_max
    @param dim
 */
void inicializarParametros(int dim);

/**
    @fn inicializacionDistinta
    @brief initialize another population for geodesics
    @param tensor
    @param n_evals
 */
void inicializacionDistinta(EspacioTiempo & tensor, int & n_evals);

/**
    @fn inicializarTensor
    @brief initialize a tensor as it was described in paper
    @param tensor. Tensor that must be initialized
    @param n_evals
*/
void inicializarTensor(EspacioTiempo& tensor, int & n_evals);

/**
    @fn particionarSubespacios
    @brief Randomly it separes a tensor in S subspaces
    @return vector<vector<int>> vector of subspaces (indexes of positions within global tensor)
*/
vector<vector<int>> particionarSubespacios();

// The following functions obtain characteristics of the LOCAL interaction between particles

/**
    @fn obtenerRadiosEnergiaCinetica
    @brief It gets kinetic energy radios for particles within a subtensor. Because of intuitive equation in paper, it does not need about global tensor
    @param subtensor. Positions in global tensor
    @param index_de_subtensor. Number of partition
    @return kinetic energy radios per particle
*/
vector<double> obtenerRadiosEnergiaCinetica(const vector<int>& subtensor, int index_de_subtensor);

/**
    @fn obtenerVectoresTransformadaLorentz
    @brief First it computes an scalar and with it, Lorentz transform vector is simulated per particle
    @param radios_e_cinetica. Values of kinetic energy radios per perticle
    @return set of Lorentz transform vector per particle
*/
vector<vector<double>> obtenerVectoresTransformadaLorentz(const vector<double>& radios_e_cinetica);

/**
    @fn obtenerK_V_dadaUnaParticula
    @brief It computes the spacetime coefficient
    @param mejorParticula best particle within subspace s
    @param particulaAleatoria random particle within subspace s
    @return spacetime coefficients vector. If it returns 0 vector, then particle i does not move
*/
vector<double> obtenerK_V_dadaUnaParticula(const Particula & mejorParticula, const Particula & particulaAleatoria);

/**
    @fn obtenerVelocidades
    @brief it computes, per particle, a vector of velocities, that asociates a vector in R^{DIMENSION} with which a local influenced particle is moved
    @param tensor. Global tensor
    @param subtensor. Particles in this tensor will get its asociated velocity tensor
    @param transformadasLorentz. Lorentz transform vectors for movement
    @return vectors of velocities per particle within subtensor
*/
vector<vector<double>> obtenerVelocidades(const EspacioTiempo& tensor, const vector<int>& subtensor,vector<vector<double>> transformadasLorentz);

// The following functions obtain INDEPENDENT or GLOBAL characteristics to the close interaction of particles

/**
    @fn obtenerDistanciasARecorrer
    @brief Using local influenced velocities per particle, it computes what distances each particle runs
    @param velocidades. Velocities vectors per particle
    @param iteracion. Iteration of algorithm. Used for a linear regression proposed wihtin paper
    @return distances each particle should run
*/
vector<vector<double>> obtenerDistanciasARecorrer(const vector<vector<double>>& velocidades, int iteracion);

/**
    @fn sgn
    @brief Implementation of sign function
    @param valor
    @return 1 if valor>0, -1 if valor<0 and 0 in other case
*/
int sgn(double valor);

/**
    @fn obtenerDeltasGeodesicas
    @brief it computes null geodesic, time geodesic and space geodesic (direction and orientation of movement) per particle
    @param tensor. Geodesics must be compute for them
    @param tensor_anterior. Needed for time geodesic
    @return deltas for orientation and direction per particle
*/
vector<vector<int>> obtenerDeltasGeodesicas(const EspacioTiempo& tensor, const EspacioTiempo& tensor_anterior);

/**
    @fn obtenerSiguientePoblacion
    @brief given a population of particles, it updates next population
    @param tensor. What should be updated
    @param distancias_recorrer. How much a particle should move in space-time
    @param deltas_geodesicas. Direction and orientation per particle
    @return new population
*/
EspacioTiempo obtenerSiguientePoblacion(const EspacioTiempo& tensor, const vector<vector<double>>& distancias_recorrer, const vector<vector<int>>& deltas_geodesicas, int & n_evals);

/**
    @fn mutarPeorParticulaCadaSubespacio
    @brief each worst particle within each subspace is mutated
    @param tensor. Population is goint to be mutate
    @param subespacios. Subspaces of tensor
    @param n_evals
*/
void mutarPeorParticulaCadaSubespacio(EspacioTiempo& tensor, const vector<vector<int>>& subespacios, int & n_evals);

// Auxiliary functions to the algorithm

/**
    @fn buscarParticulasPeorFitnessTodoEspacio
    @brief it looks for S worst particles
    @param tensor. Global tensor
    @return S worst particle positions
*/
vector<int> buscarParticulasPeorFitnessTodoEspacio(const EspacioTiempo& tensor);

/**
    @fn obtenerMejorParticulaPorSubespacio
    @brief it looks for best particle per subspace
    @param tensor. Global tensor
    @param subespacios. Subspaces
    @return best particles per subspace
*/
vector<Particula> obtenerMejorParticulaPorSubespacio(const EspacioTiempo& tensor, const vector<vector<int>> subespacios);

/**
    @fn obtenerPosicionesMejorParticulaPorSubespacion
    @brief it looks for the index of the best particle per subspace
    @param tensor. Global tensor
    @param subespacios. Subspaces
    @return best particles per subspace
*/
vector<int> obtenerPosicionesMejorParticulaPorSubespacion(const EspacioTiempo& tensor, const vector<vector<int>> subespacios);

/**
    @fn obtenerMejorParticula
    @brief given a tensor, it looks for the best particle (maximize fitness)
    @param tensor. Population
    @return best particle and its index
*/
std::pair<Particula,int> obtenerMejorParticula(const EspacioTiempo& tensor);

/**
    @fn devolverSubtensor
    @brief it looks for a subtensor given a set of indexes
    @param tensor. Global tensor
    @param indexes_subtensor. Indexes of particles that live within the new subtensor
    @return a subtensor
*/
EspacioTiempo devolverSubtensor(const EspacioTiempo& tensor, const vector<int>& indexes_subtensor);

/**
    @fn rectaDecreciente
    @brief For step length, it returns a value of a straight that has points (1,0.9) and(T_MAX,0.1)
    @param iteracion. Value X of straight
    @return Value Y of straight. This behaviour will cause smaller step length increasing iterations in GRSA
*/
double rectaDecreciente(int iteracion);

/**
    @fn compararFitnessOrdenarMenorMayorAlternativa
    @brief compares fitness between two particles
    @param p1
    @param p2
    @return True if p1 is worse than p2, False in other case
*/
bool compararFitnessOrdenarMenorMayorAlternativa(const Particula& p1, const Particula& p2);
/**
    @fn mejoraPropuesta
    @brief proposed improvement. With this I pretend to enhance diversity
    @param tensor. Space-Time which is going to be diversified
    @param n_evals
 */
void mejoraPropuesta(EspacioTiempo & tensor, int & n_evals);

/**
    @fn BL_dada_solucion_inicial
    @brief It applies Local Search with an initial solution. Parameters used are in 'aux.h'
    @param particula_inicial
    @param n_evals. Controlling max evaluations for CEC17
    @return particle and its fitness
*/
std::pair<Pesos,double> BL_dada_solucion_inicial(const Particula& particula_inicial, int & n_evals);

// Algorithm to implement
/**
    @fn GRSA
    @brief General Relativity Search Algorithm, proposed in Beiranvand et al.
    @param num_part
    @return best solution given by GRSA
*/
std::pair<Pesos,double> GRSA(int num_part);

/**
    @fn GRSA_BL
    @brief GRSA which applies Local Search each T_MAX/10 iterations to the best particule per subspace
    @param num_part
    @return best solution
*/
std::pair<Pesos,double> GRSA_BL(int num_part);

/**
    @fn GRSA_Mej
    @brief GRSA with my improvement
    @param num_part
    @return best solution
*/
std::pair<Pesos,double> GRSA_Mej(int num_part);

/**
    @fn GRSA_Todo
    @brief GRSA with my improvement, merged with BL
    @param num_part
    @return best solution
*/
std::pair<Pesos,double> GRSA_Todo(int num_part);

/**
    @fn mostrarResultadosPAlternativa
    @brief it shows results of each execution of algorithm
    @param caso
*/
void mostrarResultadosPAlternativa(const string & caso);