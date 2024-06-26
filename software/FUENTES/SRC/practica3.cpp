// Alumno: Álvaro Rodríguez Gallardo
// DNI: 77034155W
// Correo: alvaro155w@correo.ugr.es
// Grupo 2 (Miércoles, 17.30-19.30)

#include "practica3.hpp"

//************************************************************************ FUNCIONES AUXILIARES GENERAL ***************************************************************************

std::pair<Pesos,double> BL_dada_sol_inicial(const Dataset & entrenamiento, const Pesos &  sol_inicial){
    // Primero se genera la solución inicial
    Pesos solucion_actual = sol_inicial;
    // Distribuciones que se van a usar
    std::normal_distribution<double> normal(0.0,std::sqrt(VARIANZA));
    //std::uniform_real_distribution<double> mutaciones_fuertes_distr(-0.25,0.25);
    double tasa_cla_sol_actual = tasaClasificacionLeaveOneOut(entrenamiento,solucion_actual);
    double tasa_red_sol_actual = tasaReduccion(solucion_actual);
    double fitness_actual = fitness(tasa_cla_sol_actual,tasa_red_sol_actual);

    int num_caracteristicas = entrenamiento.muestras[0].caracteristicas.size();
    //std::uniforme_int_distribution<int> distribucion_entera_uniforme(0,num_caracteristicas);

    double tasa_cla_vecino;
    double tasa_red_vecino;
    double fitness_vecino;
    vector<int> indexes(num_caracteristicas);
    bool hubo_mejora = true;
    int num_evaluaciones = 1;
    bool salir = false;     // Controlar caso extremo 

    mezclarIndices(num_caracteristicas,indexes,true);

    while(hubo_mejora && num_evaluaciones<=MAX_EVALS_BL_EN_BMB_ILS){
        // Cada vez que entra al bucle, se mezclan los índices. Así, se añade aleatoriedad en la generación, evitando repetir índice
        hubo_mejora = false;
        salir = false;
        for(int j=0;j<indexes.size() && !hubo_mejora && !salir;j++){
            Pesos vecino = solucion_actual;
            mutacionBL(vecino,indexes[j],normal);
           

            tasa_cla_vecino = tasaClasificacionLeaveOneOut(entrenamiento,vecino);
            tasa_red_vecino = tasaReduccion(vecino);
            fitness_vecino = fitness(tasa_cla_vecino,tasa_red_vecino);
            
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

void generarSolucionAleatoria(Pesos & solucion, int n_pesos){
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    for(int i=0; i<n_pesos; i++){
        solucion.valores.push_back(Random::get(distribution));
    }
}

//************************************************************************ FUNCIONES AUXILIARES BMB ***************************************************************************

bool compararFitness(const Solucion& s1, const Solucion& s2){
    return s1.fitness > s2.fitness; // Orden descendente
}

//************************************************************************ FUNCIONES AUXILIARES ILS ***************************************************************************

void mutacionFuerte(Pesos & pesos, int index, std::uniform_real_distribution<double> uniforme){
    pesos.valores[index] += Random::get(uniforme);

    // Se termina truncando si fuere necesario por si sale del rango [0,1]
    if(pesos.valores[index]>1.0){
        pesos.valores[index] = 1.0;
    }
    if(pesos.valores[index]<0.0){
        pesos.valores[index] = 0.0;
    }
}

//************************************************************************ FUNCIONES AUXILIARES ES ***************************************************************************

double obtenerTemperaturaInicial(double coste_sol_inicial){
    return MU*coste_sol_inicial/(-1.0 * log(PROBAB_ACEPTAR_SOLUC_MU_POR_1_PEOR));
}

double obtenerTemperaturaCauchy(double t_k, double beta){
    return t_k / (1+(beta*t_k));
}

double obtenerBeta(double t_inicial, int max_vecinos, double t_final){
    return (t_inicial-t_final)/((NUM_ITERACIONES_ES/max_vecinos)*t_inicial*t_final);
}

double limiteExponencial(double delta, double T){
    return exp(-delta/T);
}

//************************************************************************ ALGORITMOS ***************************************************************************

std::pair<Pesos,double> BMB(const Dataset & entrenamiento){
    vector<Solucion> soluciones;

    for(int i=0;i<N_SOLS_INICIALES_BMB_ILS;i++){
        Pesos sol_generada;
        generarSolucionAleatoria(sol_generada,entrenamiento.muestras[0].caracteristicas.size());
        std::pair<Pesos,double> sol_bl = BL_dada_sol_inicial(entrenamiento,sol_generada);
        Solucion soluc;
        soluc.pesos = sol_bl.first;
        soluc.fitness = sol_bl.second;

        soluciones.push_back(soluc);
    }

    std::sort(soluciones.begin(), soluciones.end(), compararFitness);
    Solucion mejor = soluciones.front();

    return std::make_pair(mejor.pesos,mejor.fitness);
}

std::pair<Pesos,double> ES(const Dataset & entrenamiento){
    Solucion mejor_solucion;
    generarSolucionAleatoria(mejor_solucion.pesos, entrenamiento.muestras[0].caracteristicas.size());
    mejor_solucion.fitness = fitness(tasaClasificacionLeaveOneOut(entrenamiento, mejor_solucion.pesos), tasaReduccion(mejor_solucion.pesos));

    double T = obtenerTemperaturaInicial(mejor_solucion.fitness);
    double T_FINAL = 0.001;
    const int MAX_VECINOS = LIMITE_MAX_VECINOS * mejor_solucion.pesos.valores.size();
    const int MAX_EXITOS = LIMITE_MAX_EXITOS * MAX_VECINOS;
    const double BETA = obtenerBeta(T, MAX_VECINOS, T_FINAL);

    // Distribuciones
    std::uniform_int_distribution<int> distrib_indices(0, mejor_solucion.pesos.valores.size() - 1);
    std::normal_distribution<double> normal(0.0, std::sqrt(VARIANZA));
    std::uniform_real_distribution<double> distribucion_cero_uno(0.0, 1.0);

    Pesos S = mejor_solucion.pesos;
    int num_exitos_enfriamiento_actual = 1;
    int num_iteraciones = 0;
    int n_vecinos = 0;

    while(T<T_FINAL){
        T_FINAL/=10;
    }

    while (T > T_FINAL && num_iteraciones < NUM_ITERACIONES_ES && num_exitos_enfriamiento_actual != 0) {
        n_vecinos = 0;
        num_exitos_enfriamiento_actual = 0;

        while (n_vecinos < MAX_VECINOS && num_exitos_enfriamiento_actual < MAX_EXITOS && num_iteraciones < NUM_ITERACIONES_ES) {
            Pesos vecino = mejor_solucion.pesos;
            int index_aleatorio = Random::get(distrib_indices);
            mutacionBL(vecino, index_aleatorio, normal);
            n_vecinos++;

            double fitness_vecino = fitness(tasaClasificacionLeaveOneOut(entrenamiento, vecino), tasaReduccion(vecino));
            double delta = fitness_vecino - mejor_solucion.fitness;

            if (delta > 0 || Random::get(distribucion_cero_uno) <= exp(delta / T)) {
                S = vecino;
                num_exitos_enfriamiento_actual++;
                if (fitness_vecino > mejor_solucion.fitness) {
                    mejor_solucion.pesos = S;
                    mejor_solucion.fitness = fitness_vecino;
                }
            }
            num_iteraciones++;
        }
        T = obtenerTemperaturaCauchy(T, BETA);
    }

    return std::make_pair(mejor_solucion.pesos, mejor_solucion.fitness);
}

std::pair<Pesos,double> ES_dada_sol_inicial(const Dataset& entrenamiento, const Pesos& sol_inicial){
    Solucion mejor_solucion;
    //generarSolucionAleatoria(mejor_solucion.pesos, entrenamiento.muestras[0].caracteristicas.size());
    mejor_solucion.pesos = sol_inicial;
    mejor_solucion.fitness = fitness(tasaClasificacionLeaveOneOut(entrenamiento, sol_inicial), tasaReduccion(sol_inicial));

    double T = obtenerTemperaturaInicial(mejor_solucion.fitness);
    double T_FINAL = 0.001;
    const int MAX_VECINOS = LIMITE_MAX_VECINOS * mejor_solucion.pesos.valores.size();
    const int MAX_EXITOS = LIMITE_MAX_EXITOS * MAX_VECINOS;
    const double BETA = obtenerBeta(T, MAX_VECINOS, T_FINAL);

    // Distribuciones
    std::uniform_int_distribution<int> distrib_indices(0, mejor_solucion.pesos.valores.size() - 1);
    std::normal_distribution<double> normal(0.0, std::sqrt(VARIANZA));
    std::uniform_real_distribution<double> distribucion_cero_uno(0.0, 1.0);

    Pesos S = mejor_solucion.pesos;
    int num_exitos_enfriamiento_actual = 1;
    int num_iteraciones = 0;
    int n_vecinos = 0;

    while(T<T_FINAL){
        T_FINAL/=10;
    }

    while (T > T_FINAL && num_iteraciones < MAX_EVALS_BL_EN_BMB_ILS && num_exitos_enfriamiento_actual != 0) {
        n_vecinos = 0;
        num_exitos_enfriamiento_actual = 0;

        while (n_vecinos < MAX_VECINOS && num_exitos_enfriamiento_actual < MAX_EXITOS && num_iteraciones < MAX_EVALS_BL_EN_BMB_ILS) {
            Pesos vecino = mejor_solucion.pesos;
            int index_aleatorio = Random::get(distrib_indices);
            mutacionBL(vecino, index_aleatorio, normal);
            n_vecinos++;

            double fitness_vecino = fitness(tasaClasificacionLeaveOneOut(entrenamiento, vecino), tasaReduccion(vecino));
            double delta = fitness_vecino - mejor_solucion.fitness;

            if (delta > 0 || Random::get(distribucion_cero_uno) <= exp(delta / T)) {
                S = vecino;
                num_exitos_enfriamiento_actual++;
                if (fitness_vecino > mejor_solucion.fitness) {
                    mejor_solucion.pesos = S;
                    mejor_solucion.fitness = fitness_vecino;
                }
            }
            num_iteraciones++;
        }
        T = obtenerTemperaturaCauchy(T, BETA);
    }

    return std::make_pair(mejor_solucion.pesos, mejor_solucion.fitness);
}

std::pair<Pesos,double> ILS(const Dataset & entrenamiento){
    Pesos sol_i;
    generarSolucionAleatoria(sol_i, entrenamiento.muestras[0].caracteristicas.size());
    std::pair<Pesos,double> sol_i_1_bl = BL_dada_sol_inicial(entrenamiento,sol_i);

    Solucion solucion_i;

    solucion_i.pesos = sol_i_1_bl.first;
    solucion_i.fitness = sol_i_1_bl.second;

    std::uniform_int_distribution<int> distribucion_indices_mutar(0,entrenamiento.muestras[0].caracteristicas.size()-1);
    std::uniform_real_distribution<double> mutac_fuerte(LIM_INF_DISTRIBUCION,LIM_SUP_DISTRIBUCION);

    for(int i=1;i<N_SOLS_INICIALES_BMB_ILS;i++){

        int num_indices_mutar = std::ceil(LIMITE_MUTACION_ILS * entrenamiento.muestras[0].caracteristicas.size());
        if(num_indices_mutar < N_CARACTS_DEBEN_CAMBIAR){
            num_indices_mutar = N_CARACTS_DEBEN_CAMBIAR;
        }

        // Mutación fuerte sin repetición
        std::unordered_set<int> indices_usados;
        for(int j=0;j<num_indices_mutar;j++){
            int indice_mutar;
            do {    // Sin repetición
                indice_mutar = Random::get(distribucion_indices_mutar);
            } while(indices_usados.find(indice_mutar) != indices_usados.end());

            mutacionFuerte(solucion_i.pesos,indice_mutar,mutac_fuerte);
            indices_usados.insert(indice_mutar);
        }

        solucion_i.fitness = fitness(tasaClasificacionLeaveOneOut(entrenamiento,solucion_i.pesos),tasaReduccion(solucion_i.pesos));

      
      //  Solucion solucion_i_1;
        std::pair<Pesos,double> sol_bl_i_1 = BL_dada_sol_inicial(entrenamiento,solucion_i.pesos);

        if(sol_bl_i_1.second > solucion_i.fitness){
            solucion_i.pesos = sol_bl_i_1.first;
            solucion_i.fitness = sol_bl_i_1.second;
        }
    }

    return std::make_pair(solucion_i.pesos,solucion_i.fitness);
}

std::pair<Pesos,double> ILS_ES(const Dataset & entrenamiento){
    Pesos sol_i;
    generarSolucionAleatoria(sol_i, entrenamiento.muestras[0].caracteristicas.size());
    std::pair<Pesos,double> sol_i_1_es = ES_dada_sol_inicial(entrenamiento,sol_i);

    Solucion solucion_i;

    solucion_i.pesos = sol_i_1_es.first;
    solucion_i.fitness = sol_i_1_es.second;

    std::uniform_int_distribution<int> distribucion_indices_mutar(0,entrenamiento.muestras[0].caracteristicas.size()-1);
    std::uniform_real_distribution<double> mutac_fuerte(LIM_INF_DISTRIBUCION,LIM_SUP_DISTRIBUCION);

    for(int i=1;i<N_SOLS_INICIALES_BMB_ILS;i++){

        int num_indices_mutar = std::ceil(LIMITE_MUTACION_ILS * entrenamiento.muestras[0].caracteristicas.size());
        if(num_indices_mutar < N_CARACTS_DEBEN_CAMBIAR){
            num_indices_mutar = N_CARACTS_DEBEN_CAMBIAR;
        }

        // Mutación fuerte sin repetición
        std::unordered_set<int> indices_usados;
        for(int j=0;j<num_indices_mutar;j++){
            int indice_mutar;
            do {    // Sin repetición
                indice_mutar = Random::get(distribucion_indices_mutar);
            } while(indices_usados.find(indice_mutar) != indices_usados.end());

            mutacionFuerte(solucion_i.pesos,indice_mutar,mutac_fuerte);
            indices_usados.insert(indice_mutar);
        }
        
        solucion_i.fitness = fitness(tasaClasificacionLeaveOneOut(entrenamiento,solucion_i.pesos),tasaReduccion(solucion_i.pesos));
       // Solucion solucion_i_1;

        std::pair<Pesos,double> sol_es_i_1 = ES_dada_sol_inicial(entrenamiento,solucion_i.pesos);

        if(sol_es_i_1.second > solucion_i.fitness){
            solucion_i.pesos = sol_es_i_1.first;
            solucion_i.fitness = sol_es_i_1.second;
        }
    }

    return std::make_pair(solucion_i.pesos,solucion_i.fitness);
}

void resultadosP3(const string & caso){
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
        cout << "************************************ " << nombre << " (" <<caso<<") ************************************************" << endl;

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

            if(caso=="BMB"){
                init = std::clock();
                W = BMB(entrenamiento).first;         
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);

            } else if(caso=="ILS"){
                init = std::clock();
                W = ILS(entrenamiento).first;
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);

            } else if(caso == "ES"){
                init = std::clock();
                W = ES(entrenamiento).first;
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);

            } else if(caso == "ILS-ES"){
                init = std::clock();
                W = ILS_ES(entrenamiento).first;
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);

            } else {
                std::cerr<<"\nNO ES UN ALGORITMO DE LA PRACTICA 3"<<std::endl;
                return;
            }
 
            cout<<std::endl;
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
}