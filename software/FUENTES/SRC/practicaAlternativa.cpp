// Alumno: Álvaro Rodríguez Gallardo
// DNI: 77034155W
// Correo: alvaro155w@correo.ugr.es
// Grupo 2 (Miércoles, 17.30-19.30)

#include "practicaAlternativa.hpp"

// No constant global variables
// Number of particles
int TAM_POBLACION; //= S*h

// Dimension (number of weights)
int DIMENSION;  // depends on number of characteristics

// Number of ubspaces
int S;          // E(DIMENSION/3)

// Number of particles per subspace
int h;          // 10*DIMENSION

// Maximum of iterations
int T_MAX;       // 10*TAM_POBLACION

// T_max per particle
vector<double> T_maxs;

// T_min per particle
vector<double> T_mins;

int MAX_EVALUACIONES;

void inicializarParametros(int dim){
    DIMENSION = dim;
    S = DIMENSION / 3;
    h = DIMENSION * 10;
    TAM_POBLACION = S*h;
    T_MAX = S*h;
    MAX_EVALUACIONES = 10000*dim;
}

// First option: null position
void inicializacionDistinta(EspacioTiempo & tensor, int & n_evals){
    for(int i=0;i<TAM_POBLACION;i++){
        Particula part;

        for(int j=0;j<DIMENSION;j++){
            double valor = 0.0;
            part.posicion.valores.push_back(valor);
            part.mejor_posicion.valores.push_back(valor);
        }

        part.fitness = cec17_fitness(&part.posicion.valores[0]);
        n_evals++;
        tensor.T.push_back(part);
    }
}

void inicializarTensor(EspacioTiempo& tensor, int & n_evals){
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    double valor, T_max, T_min;

    // Limits of each value--> data normalised! Initializing T max and min values
    for(int i=0;i<DIMENSION;i++){
        T_max = Random::get(distribution);

        if(T_max == 0.0){
            T_max = 1.0;
        }
        T_min = T_max-1;

        if(T_min < 0.0){
            T_min = 0.0;
        }
        T_maxs.push_back(T_max);
        T_mins.push_back(T_min);
    }

    // Initializing population
    for(int i=0;i<TAM_POBLACION;i++){
        Particula part;

        for(int j=0;j<DIMENSION;j++){
            valor = T_mins[j] + Random::get(distribution)*(T_maxs[j]-T_mins[j]);
            
            if(valor < 0.0){
                valor = 0.0;
            }
            if(valor > 1.0){
                valor = 1.0;
            }
            part.posicion.valores.push_back(valor);
            part.mejor_posicion.valores.push_back(valor);
        }

        part.fitness = cec17_fitness(&part.posicion.valores[0]);//fitness(tasaClasificacionLeaveOneOut(entrenamiento,part.posicion),tasaReduccion(part.posicion));
        n_evals++;
        tensor.T.push_back(part);
    }

}

vector<vector<int>> particionarSubespacios(){
    vector<int> indexes;
    vector<vector<int>> subtensores;
    int num_index = 0;

    for(int i=0;i<TAM_POBLACION;i++){
        indexes.push_back(i);
    }

    // Without repeat indexes
    Random::shuffle(indexes);

    for(int i=0;i<S;i++){
        vector<int> subtensor;

        for(int j=0;j<h;j++){
            subtensor.push_back(indexes[num_index]);
            num_index++;
        }
        subtensores.push_back(subtensor);
    }
    assert(num_index==TAM_POBLACION);

    return subtensores;
}

vector<double> obtenerRadiosEnergiaCinetica(const vector<int>& subtensor, int index_de_subtensor){
    vector<double> radios_e_cinetica;
    uniform_int_distribution<int> distr(1,h);

    for(int i=0;i<h;i++){
        int I_rand = h*(index_de_subtensor-1)+Random::get(distr);

        double radio = (1.0)*(I_rand-TAM_POBLACION)/(TAM_POBLACION-1)*GM_1 + (1.0)*(1-I_rand)/(TAM_POBLACION-1)*GM_2;
        radios_e_cinetica.push_back(radio);
    }
    return radios_e_cinetica;
}

vector<vector<double>> obtenerVectoresTransformadaLorentz(const vector<double>& radios_e_cinetica){
    vector<vector<double>> vectoresTrLorentz;
    uniform_real_distribution<double> distr(0.0,1.0);
    vector<double> K_g;

    // Initialize K_g random vector
    for(int i=0;i<DIMENSION;i++){
        K_g.push_back(Random::get(distr));
    }

    for(int i=0;i<h;i++){
        double escalarLorentz = 1.0+radios_e_cinetica[i];
        vector<double> vectorParticula;
        
        for(int j=0;j<DIMENSION;j++){
            vectorParticula.push_back(escalarLorentz+(1.0-escalarLorentz)*K_g[j]);
        }
        vectoresTrLorentz.push_back(vectorParticula);
    }

    return vectoresTrLorentz;
}

vector<double> obtenerK_V_dadaUnaParticula(const Particula& mejorParticula, const Particula& particulaAleatoria) {
    assert(mejorParticula.posicion.valores.size() == particulaAleatoria.posicion.valores.size());
    
    vector<double> K_V_i;

    for (int i = 0; i < mejorParticula.posicion.valores.size(); i++) {
        K_V_i.push_back(fabs(mejorParticula.posicion.valores[i] - particulaAleatoria.posicion.valores[i]));
    }

    return K_V_i;
}


vector<vector<double>> obtenerVelocidades(const EspacioTiempo& tensor, const vector<int>& subtensor,vector<vector<double>> transformadasLorentz){
    vector<vector<double>> velocidades;
    std::uniform_real_distribution<double> distribucion_subtensor(0,subtensor.size()-1);

    // Initialize subtensor with particles;
    EspacioTiempo subtensorP;
    for(int i=0;i<subtensor.size();i++){
        subtensorP.T.push_back(tensor.T[subtensor[i]]);
    }

    for(int i=0;i<subtensorP.T.size();i++){
        vector<double> V_i;
        vector<double> transformada_i = transformadasLorentz[i];
        vector<double> K_V_i = obtenerK_V_dadaUnaParticula(obtenerMejorParticula(subtensorP).first,subtensorP.T[Random::get(distribucion_subtensor)]);

        assert(transformada_i.size()==K_V_i.size());

        for(int j=0;j<transformada_i.size();j++){
            V_i.push_back(K_V_i[j] * sqrt((1/(transformada_i[j]*transformada_i[j])) - 1));
        }

        velocidades.push_back(V_i);
    }
    return velocidades;
}

vector<vector<double>> obtenerDistanciasARecorrer(const vector<vector<double>>& velocidades, int iteracion){
    double peso = rectaDecreciente(iteracion);  // It depends only on iteration
    vector<vector<double>> lambdas;

    for(int i=0;i<velocidades.size();i++){
        vector<double> lambdas_particle_i;

        for(int j=0;j<velocidades[0].size();j++){
            lambdas_particle_i.push_back(peso*velocidades[i][j]);
        }

        lambdas.push_back(lambdas_particle_i);
    }

    return lambdas;
}

vector<vector<int>> obtenerDeltasGeodesicas(const EspacioTiempo& tensor, const EspacioTiempo& tensor_anterior){
    assert(tensor.T.size()==tensor_anterior.T.size());
    Particula mejor_global = obtenerMejorParticula(tensor).first;
    int aux_time_geo, aux_space_geo, aux_null_geo;
    vector<vector<int>> deltas_particulas;
    uniform_int_distribution<int> distrib_0_1(0,1);

    for(int i=0;i<tensor.T.size();i++){
        Particula T_i = tensor.T[i];
        Particula T_i_anterior = tensor_anterior.T[i];
        Pesos mejor_posicion_hasta_ahora_i = T_i.mejor_posicion;
        vector<int> delta_i;

        // Random vector K_f per particle, where K_f[j] \in {0,1}
        int K_f_j;

        for(int j=0;j<T_i.posicion.valores.size();j++){
            // Random value \in {0,1}
            K_f_j = Random::get(distrib_0_1);

            // Time geodesic
            aux_time_geo = sgn(T_i.posicion.valores[j]-T_i_anterior.posicion.valores[j]);
            
            // Space geodesic
            aux_space_geo = sgn(T_i.posicion.valores[j]-mejor_global.posicion.valores[j]);

            // Null geodesic
            aux_null_geo = sgn(T_i.posicion.valores[j]-mejor_posicion_hasta_ahora_i.valores[j]);
        
            // Position j of delta of particle i
            // If particle i moves within space geodesic, then it does not move within null geodesic because of Gemeral Relativity
            delta_i.push_back(-sgn(aux_time_geo + K_f_j * aux_space_geo + (1-K_f_j) * aux_null_geo));
        }

        deltas_particulas.push_back(delta_i);

    }

    return deltas_particulas;
}

int sgn(double valor){
    if(valor > 0)
        return 1;
    if(valor < 0)
        return -1;
    return 0;
}

EspacioTiempo obtenerSiguientePoblacion(const EspacioTiempo& tensor, const vector<vector<double>>& distancias_recorrer, const vector<vector<int>>& deltas_geodesicas, int & n_evals){
    EspacioTiempo nuevo_tensor = tensor;

    for(int i=0;i<tensor.T.size() && n_evals<MAX_EVALUACIONES;i++){
            for(int j=0;j<DIMENSION;j++){
                nuevo_tensor.T[i].posicion.valores[j] = tensor.T[i].posicion.valores[j]+distancias_recorrer[i][j]*deltas_geodesicas[i][j];

                if(nuevo_tensor.T[i].posicion.valores[j] < 0.0){
                    nuevo_tensor.T[i].posicion.valores[j] = 0.0;
                }
                if(nuevo_tensor.T[i].posicion.valores[j] > 1.0){
                    nuevo_tensor.T[i].posicion.valores[j] = 1.0;
                }

            }
        
            nuevo_tensor.T[i].fitness = cec17_fitness(&nuevo_tensor.T[i].posicion.valores[0]);//fitness(tasaClasificacionLeaveOneOut(entrenamiento,nuevo_tensor.T[i].posicion),tasaReduccion(nuevo_tensor.T[i].posicion));
            n_evals++;

            if(nuevo_tensor.T[i].fitness > tensor.T[i].fitness){
                nuevo_tensor.T[i].mejor_posicion = nuevo_tensor.T[i].posicion;
            }
        
    }

    return nuevo_tensor;
}

void mutarPeorParticulaCadaSubespacio(EspacioTiempo& tensor, const vector<vector<int>>& subespacios, int & n_evals){
    std::uniform_int_distribution<int> distr(0,1);

    vector<int> alphas_1, alphas_2;

    // Initializing alpha
    for(int i=0;i<DIMENSION;i++){
        alphas_1.push_back(Random::get(distr));
        alphas_2.push_back(1-alphas_1[i]);
    }

    // Auxiliar data structures. R: the worst particles within ALL spacetime. g: best particle per subspace
    std::pair<vector<Particula>,vector<int>> R;
    vector<Particula> g = obtenerMejorParticulaPorSubespacio(tensor,subespacios);

    R.second = buscarParticulasPeorFitnessTodoEspacio(tensor);

    for(int i=0;i<R.second.size();i++){
        R.first.push_back(tensor.T[R.second[i]]);
    }
    
    for(int i=0; i<S && n_evals<MAX_EVALUACIONES;i++){
            for(int j=0;j<DIMENSION;j++){
                tensor.T[R.second[i]].posicion.valores[j] = alphas_1[j] * R.first[i].posicion.valores[j] + alphas_2[j] * g[i].posicion.valores[j];
            }

            double nuevo_fitness = cec17_fitness(&tensor.T[R.second[i]].posicion.valores[0]);//fitness(tasaClasificacionLeaveOneOut(entrenamiento,tensor.T[R.second[i]].posicion),tasaReduccion(tensor.T[R.second[i]].posicion));
            n_evals++;

            if(nuevo_fitness > tensor.T[R.second[i]].fitness){
                tensor.T[R.second[i]].mejor_posicion = tensor.T[R.second[i]].posicion;
            }

            tensor.T[R.second[i]].fitness = nuevo_fitness;
        
    }

}

bool compararFitnessOrdenarMenorMayorAlternativa(const Particula& p1, const Particula& p2){
    return p1.fitness < p2.fitness;
}

vector<int> buscarParticulasPeorFitnessTodoEspacio(const EspacioTiempo& tensor){
    vector<Particula> copia = tensor.T;

    std::stable_sort(copia.begin(),copia.end(),compararFitnessOrdenarMenorMayorAlternativa);
    vector<int> indexes;

    for(int i=0;i<S;i++){
        indexes.push_back(std::distance(tensor.T.begin(), std::find(tensor.T.begin(), tensor.T.end(), copia[i])));
    }

    return indexes;
}

vector<Particula> obtenerMejorParticulaPorSubespacio(const EspacioTiempo& tensor, const vector<vector<int>> subespacios){
    vector<Particula> mejores;

    for(int i=0;i<S;i++){
        EspacioTiempo subtensor = devolverSubtensor(tensor,subespacios[i]);
        mejores.push_back(obtenerMejorParticula(subtensor).first);
    }

    return mejores;
}

vector<int> obtenerPosicionesMejorParticulaPorSubespacion(const EspacioTiempo& tensor, const vector<vector<int>> subespacios){
    vector<int> mejores;

    for(int i=0;i<S;i++){
        EspacioTiempo subtensor = devolverSubtensor(tensor,subespacios[i]);
        mejores.push_back(obtenerMejorParticula(subtensor).second);
    }

    return mejores;
}

std::pair<Particula,int> obtenerMejorParticula(const EspacioTiempo& tensor){
    Particula mejor = tensor.T[0];
    int index_mejor = 0;

    for(int i=1;i<tensor.T.size();i++){
        if(mejor.fitness < tensor.T[i].fitness){
            mejor = tensor.T[i];
            index_mejor = i;
        }
    }

    return make_pair(mejor,index_mejor);
}

EspacioTiempo devolverSubtensor(const EspacioTiempo& tensor, const vector<int>& indexes_subtensor){
    EspacioTiempo subtensor;

    for(int i=0;i<indexes_subtensor.size();i++){
        subtensor.T.push_back(tensor.T[indexes_subtensor[i]]);
    }

    return subtensor;
}

double rectaDecreciente(int iteracion){
    // (1,0.9) to (T_MAX,0.1)
    // 0.9=m+n; 0.1=T_MAX*m+n--> 0.8=(1-T_MAX)*m--> m=0.8/(1-T_MAX) --> n=0.9-m=0.9-
    double m = (0.9-0.1)/(1-T_MAX);  // < 0. Decreases
    double n = 0.9-m;

    return  (m*iteracion)+n;
}

void mejoraPropuesta(EspacioTiempo &tensor, int & n_evals) {
    // Firstly maximum radius is computed
    double max_radio = 0.0;
    for (const auto& particula : tensor.T) {
        if (particula.fitness > max_radio) {
            max_radio = particula.fitness;
        }
    }

    // If a particle is eliminated, we will create another one randomly, but we do not apply this enhance to it
    std::vector<int> indices_a_eliminar;

    for (int i = 0; i < tensor.T.size(); ++i) {
        const Particula& particula_referencia = tensor.T[i];
        const auto& centro_bola = particula_referencia.posicion.valores;
        double radio_bola = std::max(0.05, (particula_referencia.fitness / max_radio) - (T_maxs[i] / 2.0));
        
        for (int j = 0; j < tensor.T.size(); ++j) {
            if (i != j && particula_referencia.fitness > tensor.T[j].fitness) {
                const auto& pos_j = tensor.T[j].posicion.valores;
                bool dentro_de_bola = true;
                
                for (int k = 0; k < DIMENSION && dentro_de_bola; ++k) {
                    double distancia = std::abs(centro_bola[k] - pos_j[k]);     // Used norm = l1
                    if (distancia > radio_bola) {
                        dentro_de_bola = false;
                    }
                }
                
                if (dentro_de_bola) {
                    indices_a_eliminar.push_back(j);
                }
            }
        }
    }

    // Filling missing particles
    for (int idx : indices_a_eliminar) {
        if(n_evals<MAX_EVALUACIONES){
            Particula nueva_particula;

            for(int p=0;p<DIMENSION;p++){
                std::uniform_real_distribution<double> distribution(T_mins[p], T_maxs[p]);
                nueva_particula.posicion.valores.push_back(Random::get(distribution));
                nueva_particula.mejor_posicion.valores.push_back(nueva_particula.posicion.valores[p]);
            }

            nueva_particula.fitness = cec17_fitness(&nueva_particula.posicion.valores[0]);
            n_evals++;
            tensor.T[idx] = nueva_particula;
        }
    }
}

std::pair<Pesos,double> BL_dada_solucion_inicial(const Particula& particula_inicial, int & n_evals){
    Pesos solucion_actual = particula_inicial.posicion;
    // Distributions to use
    std::normal_distribution<double> normal(0.0,std::sqrt(VARIANZA));

    double fitness_actual = particula_inicial.fitness;

    int num_caracteristicas = DIMENSION;

    double tasa_cla_vecino;
    double tasa_red_vecino;
    double fitness_vecino;
    vector<int> indexes(num_caracteristicas);
    bool hubo_mejora = true;
    int num_evaluaciones = 1;
    bool salir = false;     // Extreme situation control

    mezclarIndices(num_caracteristicas,indexes,true);

    while(hubo_mejora && num_evaluaciones<=MAX_EVALS_BL && n_evals<MAX_EVALUACIONES){
        // Every time you enter the loop, the indices are shuffled. Thus, randomness is added to the generation, avoiding repeating the index
        hubo_mejora = false;
        salir = false;
        for(int j=0;j<indexes.size() && !hubo_mejora && !salir && n_evals<MAX_EVALUACIONES;j++){
            Pesos vecino = solucion_actual;
            mutacionBL(vecino,indexes[j],normal);
           

           // tasa_cla_vecino = //tasaClasificacionLeaveOneOut(entrenamiento,vecino);
           // tasa_red_vecino = tasaReduccion(vecino);
            fitness_vecino = cec17_fitness(&vecino.valores[0]);//fitness(tasa_cla_vecino,tasa_red_vecino);
            n_evals++;
            num_evaluaciones++;

            if(fitness_vecino > fitness_actual){
                fitness_actual = fitness_vecino;
                solucion_actual = vecino;
                hubo_mejora = true;
            }
        }
        if (hubo_mejora){
            mezclarIndices(num_caracteristicas,indexes,false);
        }

    }

    return std::make_pair(solucion_actual,fitness_actual);
}

std::pair<Pesos,double> GRSA(int dim,int num_part){
    EspacioTiempo T_anterior, tensor;
    int t=1;
    vector<vector<int>> subespacios;
    int n_evals = 0;

    inicializarParametros(dim);
    inicializacionDistinta(T_anterior,n_evals);
    inicializarTensor(tensor,n_evals);

    subespacios = particionarSubespacios();

    while(t <= T_MAX && n_evals<=MAX_EVALUACIONES){
        vector<vector<double>> velocs_iteracion;

        for(int i=0;i<S;i++){
            vector<double> radios_energia_cinetica_sub_s = obtenerRadiosEnergiaCinetica(subespacios[i],i+1);
            vector<vector<double>> transformadas_Lorentz = obtenerVectoresTransformadaLorentz(radios_energia_cinetica_sub_s);
            vector<vector<double>> velocidades_sub_s = obtenerVelocidades(tensor,subespacios[i],transformadas_Lorentz);

            for(const auto& veloc : velocidades_sub_s){
                velocs_iteracion.push_back(veloc);
            }
        }
        vector<vector<double>> step_lengths = obtenerDistanciasARecorrer(velocs_iteracion,t);
        vector<vector<int>> deltas = obtenerDeltasGeodesicas(tensor,T_anterior);
        T_anterior = tensor;

        if(n_evals<MAX_EVALUACIONES){
            // Mutación propuesta por artículo  
            tensor = obtenerSiguientePoblacion(tensor,step_lengths,deltas,n_evals);
            if(n_evals<MAX_EVALUACIONES)
                mutarPeorParticulaCadaSubespacio(tensor,subespacios,n_evals); 
        }
        escribirFitnessParaConvergencia(num_part,obtenerMejorParticula(tensor).first.fitness,t,"GRSA","convergencia/GRSA/");

        t++;
    }
    Particula mejorParticula = obtenerMejorParticula(tensor).first;

    return std::make_pair(mejorParticula.posicion,mejorParticula.fitness);
}

std::pair<Pesos,double> GRSA_BL(int dim, int num_part){
    EspacioTiempo T_anterior, tensor;
    int t=1;
    vector<vector<int>> subespacios;
    int n_evals = 0;
    inicializarParametros(dim);
    inicializacionDistinta(T_anterior,n_evals);    
    inicializarTensor(tensor,n_evals);
    int distance = T_MAX / 10;
    
    subespacios = particionarSubespacios();

    while(t <= T_MAX && n_evals<=MAX_EVALUACIONES){
        vector<vector<double>> velocs_iteracion;

        for(int i=0;i<S;i++){
            vector<double> radios_energia_cinetica_sub_s = obtenerRadiosEnergiaCinetica(subespacios[i],i+1);
            vector<vector<double>> transformadas_Lorentz = obtenerVectoresTransformadaLorentz(radios_energia_cinetica_sub_s);
            vector<vector<double>> velocidades_sub_s = obtenerVelocidades(tensor,subespacios[i],transformadas_Lorentz);
            //velocs_iteracion.push_back(velocidades_sub_s);
            for(const auto& veloc : velocidades_sub_s){
                velocs_iteracion.push_back(veloc);
            }
        }

        vector<vector<double>> step_lengths = obtenerDistanciasARecorrer(velocs_iteracion,t);
        vector<vector<int>> deltas = obtenerDeltasGeodesicas(tensor,T_anterior);
        T_anterior = tensor;

        if(n_evals<MAX_EVALUACIONES){
            // Mutación propuesta por artículo  
            tensor = obtenerSiguientePoblacion(tensor,step_lengths,deltas,n_evals);
            if(n_evals<MAX_EVALUACIONES)
                mutarPeorParticulaCadaSubespacio(tensor,subespacios,n_evals); 
        } 

        if(t % distance == 0 && n_evals<MAX_EVALUACIONES){
            // Local search to the best particle per subspace
            vector<int> mejores_indices = obtenerPosicionesMejorParticulaPorSubespacion(tensor,subespacios);

            for(int i=0;i<mejores_indices.size() && n_evals<MAX_EVALUACIONES;i++){
                std::pair<Pesos,double> tras_BL = BL_dada_solucion_inicial(tensor.T[mejores_indices[i]],n_evals);
                tensor.T[mejores_indices[i]].posicion = tras_BL.first;
                if(tras_BL.second > tensor.T[mejores_indices[i]].fitness){
                    tensor.T[mejores_indices[i]].mejor_posicion = tras_BL.first;
                }
                tensor.T[mejores_indices[i]].fitness = tras_BL.second;
            }

        }

        escribirFitnessParaConvergencia(num_part,obtenerMejorParticula(tensor).first.fitness,t,"GRSA-BL","convergencia/GRSA_BL/");
        
        t++;
        
    }

    Particula mejorParticula = obtenerMejorParticula(tensor).first;

    return std::make_pair(mejorParticula.posicion,mejorParticula.fitness);
}

std::pair<Pesos,double> GRSA_Mej(int dim, int num_part){
    EspacioTiempo T_anterior, tensor;
    int t=1;
    vector<vector<int>> subespacios;
    int n_evals = 0;

    inicializarParametros(dim);
    inicializacionDistinta(T_anterior,n_evals);
    inicializarTensor(tensor,n_evals);

    subespacios = particionarSubespacios();

    while(t <= T_MAX && n_evals<=MAX_EVALUACIONES){
        vector<vector<double>> velocs_iteracion;

        for(int i=0;i<S;i++){
            vector<double> radios_energia_cinetica_sub_s = obtenerRadiosEnergiaCinetica(subespacios[i],i+1);
            vector<vector<double>> transformadas_Lorentz = obtenerVectoresTransformadaLorentz(radios_energia_cinetica_sub_s);
            vector<vector<double>> velocidades_sub_s = obtenerVelocidades(tensor,subespacios[i],transformadas_Lorentz);

            for(const auto& veloc : velocidades_sub_s){
                velocs_iteracion.push_back(veloc);
            }
        }

        vector<vector<double>> step_lengths = obtenerDistanciasARecorrer(velocs_iteracion,t);
        vector<vector<int>> deltas = obtenerDeltasGeodesicas(tensor,T_anterior);

        T_anterior = tensor;  

        if(n_evals<MAX_EVALUACIONES){
            // Mutación propuesta por artículo 
            tensor = obtenerSiguientePoblacion(tensor,step_lengths,deltas,n_evals);

            // Mejora propuesta por mi para diversidad
            if(n_evals<MAX_EVALUACIONES)
                mejoraPropuesta(tensor,n_evals);
            if(n_evals<MAX_EVALUACIONES)
                mutarPeorParticulaCadaSubespacio(tensor,subespacios,n_evals); 
        }

        escribirFitnessParaConvergencia(num_part,obtenerMejorParticula(tensor).first.fitness,t,"GRSA-Mej","convergencia/GRSA_Mej/");

        t++;

    }

    Particula mejorParticula = obtenerMejorParticula(tensor).first;

    return std::make_pair(mejorParticula.posicion,mejorParticula.fitness);
}

std::pair<Pesos,double> GRSA_Todo(int dim,int num_part){
    EspacioTiempo T_anterior, tensor;
    int t=1;
    vector<vector<int>> subespacios;
    int n_evals = 0;

    inicializarParametros(dim);
    inicializacionDistinta(T_anterior,n_evals);
    inicializarTensor(tensor,n_evals);
    int distance = T_MAX / 10;
    
    subespacios = particionarSubespacios();

    while(t <= T_MAX && n_evals<MAX_EVALUACIONES){
        vector<vector<double>> velocs_iteracion;
        for(int i=0;i<S;i++){
            vector<double> radios_energia_cinetica_sub_s = obtenerRadiosEnergiaCinetica(subespacios[i],i+1);
            vector<vector<double>> transformadas_Lorentz = obtenerVectoresTransformadaLorentz(radios_energia_cinetica_sub_s);
            vector<vector<double>> velocidades_sub_s = obtenerVelocidades(tensor,subespacios[i],transformadas_Lorentz);

            for(const auto& veloc : velocidades_sub_s){
                velocs_iteracion.push_back(veloc);
            }
        }
        vector<vector<double>> step_lengths = obtenerDistanciasARecorrer(velocs_iteracion,t);
        vector<vector<int>> deltas = obtenerDeltasGeodesicas(tensor,T_anterior);
        T_anterior = tensor;

        if(n_evals<MAX_EVALUACIONES){
            // Mutación propuesta por artículo  
            tensor = obtenerSiguientePoblacion(tensor,step_lengths,deltas,n_evals);
            // Mejora propuesta por mi para diversidad
            if(n_evals<MAX_EVALUACIONES)
                mejoraPropuesta(tensor,n_evals);
            if(n_evals<MAX_EVALUACIONES)
                mutarPeorParticulaCadaSubespacio(tensor,subespacios,n_evals); 
        } 

        if(t % distance == 0 && n_evals<MAX_EVALUACIONES){
            // Local search to the best particle per subspace
            vector<int> mejores_indices = obtenerPosicionesMejorParticulaPorSubespacion(tensor,subespacios);

            for(int i=0;i<mejores_indices.size() && n_evals<MAX_EVALUACIONES;i++){
                std::pair<Pesos,double> tras_BL = BL_dada_solucion_inicial(tensor.T[mejores_indices[i]],n_evals);
                tensor.T[mejores_indices[i]].posicion = tras_BL.first;
                if(tras_BL.second > tensor.T[mejores_indices[i]].fitness){
                    tensor.T[mejores_indices[i]].mejor_posicion = tras_BL.first;
                }
                tensor.T[mejores_indices[i]].fitness = tras_BL.second;
            }

        }
        
        escribirFitnessParaConvergencia(num_part,obtenerMejorParticula(tensor).first.fitness,t,"GRSA-Todo","convergencia/GRSA_Todo/");

        t++;
    }

    Particula mejorParticula = obtenerMejorParticula(tensor).first;

    return std::make_pair(mejorParticula.posicion,mejorParticula.fitness);
}

void mostrarResultadosPAlternativa(const string & caso){
    vector<int> dimensiones;
    dimensiones.push_back(10);
    dimensiones.push_back(30);
    dimensiones.push_back(50);
    
    cout<<"Algoritmo: "<<caso<<std::endl;
    cout<<"--------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;

    for(auto dim : dimensiones){
        cout<<"Dimensión del problema: "<<dim<<std::endl;

        for(int funcid = 1; funcid<=30;funcid++){

            std::pair<Pesos,double> sol;
            auto init=0.0,fin=0.0;

            init = std::clock();
            if(caso=="GRSA"){
                cec17_init("GRSA",funcid,dim);

                sol = GRSA(dim,funcid);
            }
            if(caso=="GRSA-BL"){
                cec17_init("GRSA_BL",funcid,dim);
                sol = GRSA_BL(dim,funcid);
            }
            if(caso=="GRSA-Mej"){
                cec17_init("GRSA_Mej",funcid,dim);
                sol = GRSA_Mej(dim,funcid);
            }
            if(caso=="GRSA-Todo"){
                cec17_init("GRSA_Todo",funcid,dim);
                sol = GRSA_Todo(dim,funcid);
            }
            fin = std::clock();
            
            double fitness = sol.second;
            double error = cec17_error(sol.second);
            double tiempo = (fin-init)/CLOCKS_PER_SEC;

            cout<<"Función "<<funcid<<". Fitness medio: "<<fitness<<". Error medio: "<<error<<". Tiempo(s): "<<tiempo<<std::endl;

        }        

        cout<<"--------------------------------------------------------------------------------------------------------------------------------------------"<<std::endl;

    }

}

/*void mostrarResultadosPAlternativa(const string & caso){
    string nombre;

    for(int i=0; i<NUM_CONJUNTOS_DATOS;i++){
        if(i==0)
            nombre = BCANCER;
        if(i==1)
            nombre = ECOLI;
        if(i==2)
            nombre = PARKINSON;

        // Se leen los ficheros
        vector<Dataset> particiones;
        Dataset part1 = lecturaFichero("../BIN/DATA/"+nombre+"1.arff");
        Dataset part2 = lecturaFichero("../BIN/DATA/"+nombre+"2.arff");
        Dataset part3 = lecturaFichero("../BIN/DATA/"+nombre+"3.arff");
        Dataset part4 = lecturaFichero("../BIN/DATA/"+nombre+"4.arff");
        Dataset part5 = lecturaFichero("../BIN/DATA/"+nombre+"5.arff");

        particiones.push_back(part1);
        particiones.push_back(part2);
        particiones.push_back(part3);
        particiones.push_back(part4);
        particiones.push_back(part5);

        // Normalizo ahora. Antes el error era que hacía una normalización local (máximos y mínimos locales de cada partición, no en general)
        normalizarDatos(particiones);

        cout << endl << endl;
        cout << "************************************ " << nombre << " ("<<caso<<") ************************************************" << endl;

        cout << endl << "....................................................................................................." << endl;
        cout << "::: Particion ::: Tasa de Clasificacion (%) ::: Tasa de Reduccion (%) ::: Fitness ::: Tiempo (s) :::" << endl;
        cout << "....................................................................................................." << endl;

        double tasa_clas_acumulada = 0.0;
        double tasa_red_acumulada = 0.0;
        double fitness_acumulada = 0.0;
        double tiempo_total = 0.0;

        // Para distintas particiones se ejecuta 1NN original
        for(int i=0; i < particiones.size(); i++){
            // Tomamos una partición de las NUM_PARTICIONES para test, guardando el resto para entrenamiento, que se unirán en un único dataset
            Dataset test = particiones[i];
            vector<Dataset> entrenam_particiones;
            for(int j=0; j < particiones.size();j++){
                if(j!=i){
                    entrenam_particiones.push_back(particiones[j]);
                }
            }
            Dataset entrenamiento = unirDatasets(entrenam_particiones);

            Pesos W;
            W.valores = std::vector<double>(test.muestras[0].caracteristicas.size());
            auto init=0.0,fin=0.0;
            double tasa_cla_i,tasa_red_i,fitness_i;

            if(caso == "GRSA"){
                init = std::clock();
                W = GRSA(entrenamiento,nombre,i+1).first;
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);
            }
            
            if(caso == "GRSA-BL"){
                init = std::clock();
                W = GRSA_BL(entrenamiento,nombre,i+1).first;
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);
            }

            if(caso == "GRSA-Mej"){
                init = std::clock();
                W = GRSA_Mej(entrenamiento,nombre,i+1).first;
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);
            }

            if(caso == "GRSA-Todo"){
                init = std::clock();
                W = GRSA_Todo(entrenamiento,nombre,i+1).first;
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);
            }

           

            double tiempo_i = (fin-init)/CLOCKS_PER_SEC;

            // Se acumula el total
            tasa_clas_acumulada+=tasa_cla_i;
            tasa_red_acumulada+=tasa_red_i;
            fitness_acumulada+=fitness_i;
            tiempo_total+=tiempo_i;

            // Para iteración i, se muestran los resultados
            cout << fixed << setprecision(5);
            cout << ":::" << setw(6) << (i+1) << setw(8) << ":::" << setw(15) << tasa_cla_i << setw(15) << ":::" << setw(13) << tasa_red_i;
            cout << setw(13) << ":::" << setw(7) << fitness_i << setw(5) << "::: " << setw(9) << tiempo_i << std::setw(7) << ":::" << endl;
            escribirCSV(tasa_cla_i,tasa_red_i,fitness_i,tiempo_i,false);
        }

        cout << ":::" << setw(8) << "MEDIA" << setw(6) << ":::" << setw(15) << (tasa_clas_acumulada/NUM_PARTICIONES) << setw(15) << ":::" << setw(13) << (tasa_red_acumulada/NUM_PARTICIONES);
        cout << setw(13) << ":::" << setw(7) << (fitness_acumulada/NUM_PARTICIONES) << setw(5) << "::: " << setw(9) << (tiempo_total/NUM_PARTICIONES) << std::setw(7) << ":::" << endl;  
        cout << "....................................................................................................." << endl << endl;
        escribirCSV(tasa_clas_acumulada/NUM_PARTICIONES,tasa_red_acumulada/NUM_PARTICIONES,fitness_acumulada/NUM_PARTICIONES,tiempo_total/NUM_PARTICIONES,true);
    }
}*/