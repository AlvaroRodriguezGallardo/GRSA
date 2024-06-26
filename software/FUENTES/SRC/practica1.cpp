// Alumno: Álvaro Rodríguez Gallardo
// DNI: 77034155W
// Correo: alvaro155w@correo.ugr.es
// Grupo 2 (Miércoles, 17.30-19.30)
#include "practica1.hpp"

//*******************************************************************************************************************************************************
//****************************************************************** 1-NN *******************************************************************************
//*******************************************************************************************************************************************************

string algoritmo1NN(const Muestra & m, const Dataset & entrenamiento, const Pesos & pesos){
    if(entrenamiento.muestras.empty()){
        cout<<"El conjunto de entrenamiento está vacío"<<std::endl;
        return "";
    }

    int index_min_distancia = 0;
    double min_distancia = distanciaEuclideaPonderada(m,entrenamiento.muestras[index_min_distancia],pesos);
    double distancia = min_distancia;

    for(int i=1; i<entrenamiento.muestras.size();i++){
        distancia = distanciaEuclideaPonderada(m,entrenamiento.muestras[i],pesos);
        if(distancia < min_distancia){
            index_min_distancia = i;
            min_distancia = distancia;
        }
    }

    return entrenamiento.muestras[index_min_distancia].clase;
}


//*******************************************************************************************************************************************************
//************************************************************** BÚSQUEDA LOCAL (BL) ********************************************************************
//*******************************************************************************************************************************************************

// CORREGIDA DE LA P1 para la P2
std::pair<Pesos,double> BusquedaLocal(const Dataset & entrenamiento){
    // Primero se genera la solución inicial
    Pesos solucion_actual;
    // Distribuciones que se van a usar
    std::uniform_real_distribution<double> distribution(0.0,1.0);
    //std::uniform_int_distribution<> distr_enteros(0,entrenamiento.muestras[0].caracteristicas.size()-1);
    std::normal_distribution<double> normal(0.0,std::sqrt(VARIANZA));
   // solucion.valores.resize(entrenamiento.muestras[0].caracteristicas.size());
    
    for(int i=0; i < entrenamiento.muestras[0].caracteristicas.size(); i++){
        solucion_actual.valores.push_back(Random::get(distribution));
    }

    double tasa_cla_sol_actual = tasaClasificacionLeaveOneOut(entrenamiento,solucion_actual);
    double tasa_red_sol_actual = tasaReduccion(solucion_actual);
    double fitness_actual = fitness(tasa_cla_sol_actual,tasa_red_sol_actual);
    int n_vecinos_generados_sin_mejora = 0;
    int num_caracteristicas = entrenamiento.muestras[0].caracteristicas.size();
    double tasa_cla_vecino;
    double tasa_red_vecino;
    double fitness_vecino;
    vector<int> indexes(num_caracteristicas);
    bool hubo_mejora = false;
    int num_iterac = 0;

    mezclarIndices(num_caracteristicas,indexes,true);

    while(num_iterac < MAXIMO_ITERACIONES && n_vecinos_generados_sin_mejora < LIMITE * num_caracteristicas){
        // Cada vez que entra al bucle, se mezclan los índices. Así, se añade aleatoriedad en la generación, evitando repetir índice
        hubo_mejora = false;
        for(int j=0;j<indexes.size() && !hubo_mejora && num_iterac<MAXIMO_ITERACIONES;j++){
            Pesos vecino = solucion_actual;
            mutacionBL(vecino,indexes[j],normal);

            tasa_cla_vecino = tasaClasificacionLeaveOneOut(entrenamiento,vecino);
            tasa_red_vecino = tasaReduccion(vecino);
            fitness_vecino = fitness(tasa_cla_vecino,tasa_red_vecino);
            
            num_iterac++;

            if(fitness_vecino > fitness_actual){
                fitness_actual = fitness_vecino;
                solucion_actual = vecino;
                hubo_mejora = true;
            } else {
                n_vecinos_generados_sin_mejora++;
            }
        }
        if (hubo_mejora){
            n_vecinos_generados_sin_mejora = 0;
            mezclarIndices(num_caracteristicas,indexes,false);
        }

    }

    return std::make_pair(solucion_actual,fitness_actual);

}

void mutacionBL(Pesos & pesos, int index, std::normal_distribution<double> normal){
    pesos.valores[index] += Random::get(normal);

    // Se termina truncando si fuere necesario por si sale del rango [0,1]
    if(pesos.valores[index]>1.0){
        pesos.valores[index] = 1.0;
    }
    if(pesos.valores[index]<0.0){
        pesos.valores[index] = 0.0;
    }
}

void mezclarIndices(int n_pesos, vector<int> & indexes, bool inicializar){
    if(inicializar){
        for(int i=0;i<n_pesos;i++){
            indexes[i] = i;
        }
    }
    Random::shuffle(indexes);
}

//*******************************************************************************************************************************************************
//******************************************************************* RELIEF ****************************************************************************
//*******************************************************************************************************************************************************

std::pair<Pesos,double> GreedyRelief(const Dataset & entrenamiento){
    Pesos solucion_actual;
    // Se inicializa el vector de pesos a 0 todos sus elementos
    for(int i=0;i<entrenamiento.muestras[0].caracteristicas.size();i++){
        solucion_actual.valores.push_back(0.0);
    }

    for(int i=0;i<entrenamiento.muestras.size();i++){
        Muestra enemigoMasCercano = buscarEnemigoMasCercano(entrenamiento,entrenamiento.muestras[i]);   //Bien
        Muestra amigoMasCercano = buscarAmigoMasCercano(entrenamiento,entrenamiento.muestras[i]);       // Bien

        if(amigoMasCercano.caracteristicas.size()>0){
            // Sumar diferencia de distancias
            for(int j=0;j<amigoMasCercano.caracteristicas.size();j++){
                solucion_actual.valores[j] = solucion_actual.valores[j] + std::abs(entrenamiento.muestras[i].caracteristicas[j]-enemigoMasCercano.caracteristicas[j]) - std::abs(entrenamiento.muestras[i].caracteristicas[j]-amigoMasCercano.caracteristicas[j]);
            }
        }
    }

    double w_max = *std::max_element(solucion_actual.valores.begin(),solucion_actual.valores.end());

    for(int j=0;j<solucion_actual.valores.size();j++){
        if(solucion_actual.valores[j] < 0.0){
            solucion_actual.valores[j] = 0.0;
        } else {
            solucion_actual.valores[j] = solucion_actual.valores[j]/w_max;
        }
    }
  
    return std::make_pair(solucion_actual,fitness(tasaClasificacionLeaveOneOut(entrenamiento,solucion_actual),tasaReduccion(solucion_actual)));
  }

Muestra buscarEnemigoMasCercano(const Dataset & entrenamiento, const Muestra & m){
    // Primero juntamos en un dataset aquellos candidatos a enemigo (clase distinta)
    Dataset candidatosEnemigo;

    for(int i=0;i<entrenamiento.muestras.size();i++){

        if(entrenamiento.muestras[i].clase != m.clase){
            candidatosEnemigo.muestras.push_back(entrenamiento.muestras[i]);
        }
        
        //candidatosEnemigo.atributos.push_back(entrenamiento.atributos[i]);  // No es necesaria esta línea para la función
    }

    // Cuando se tiene el dataset con los posibles enemigos, se busca el enemigo más cercano
    int index_enemigo_buscado = buscarMuestraAMenorDistancia(candidatosEnemigo,m);

    return candidatosEnemigo.muestras[index_enemigo_buscado];
}


Muestra buscarAmigoMasCercano(const Dataset & entrenamiento, const Muestra & m){
    // Primero juntamos en un dataset aquellos candidatos a amigos que no sean m
    Dataset candidatosAmigo;
    Muestra amigo;  // Si no hubiese candidatos a amigo por falta de datos, se identifica porque amigo.caracteristicas.size()==0
    int index_amigo_buscado = -1;

    for(int i=0;i<entrenamiento.muestras.size();i++){

        if(entrenamiento.muestras[i].clase == m.clase  && m.caracteristicas != entrenamiento.muestras[i].caracteristicas){
            candidatosAmigo.muestras.push_back(entrenamiento.muestras[i]);
        }      
        //    candidatosAmigo.atributos.push_back(entrenamiento.atributos[i]);    // No necesaria
    }

    // Se busca la muestra que sea amigo más cercano de m, si existe
    if(candidatosAmigo.muestras.size()>0){
        index_amigo_buscado = buscarMuestraAMenorDistancia(candidatosAmigo,m);
        amigo = candidatosAmigo.muestras[index_amigo_buscado];
    }

    return amigo;
}

//*******************************************************************************************************************************************************
//************************************************************** ESTADÍSTICOS ***************************************************************************
//*******************************************************************************************************************************************************

double tasaClasificacion(const Dataset & entrenamiento, const Dataset & test, const Pesos & pesos){
    int n_bien_clasificada = 0;
    string prediccion;

    for(int i=0; i<test.muestras.size();i++){
        prediccion = algoritmo1NN(test.muestras[i],entrenamiento,pesos);
        if(prediccion==test.muestras[i].clase){
            n_bien_clasificada++;
        }
    }
    return (100.0*static_cast<double>(n_bien_clasificada)) / static_cast<double>(test.muestras.size());
}

double tasaClasificacionLeaveOneOut(const Dataset &entrenamiento, const Pesos &pesos) {
    int num_instancias_bien_clas = 0;

    for (int i = 0; i < entrenamiento.muestras.size(); i++) {
        Muestra ejemploActual = entrenamiento.muestras[i];
        Dataset aux = entrenamiento;

        aux.muestras.erase(aux.muestras.begin() + i);

        if (algoritmo1NN(ejemploActual, aux, pesos) == ejemploActual.clase) {
            num_instancias_bien_clas++;
        }

    }

    return (100.0 * static_cast<double>(num_instancias_bien_clas)) / static_cast<double>(entrenamiento.muestras.size());

}


double tasaReduccion(const Pesos & pesos){
    const double UMBRAL = 0.1;
    int n_menor_umbral = 0;

    for(int i=0;i<pesos.valores.size();i++){
        if(pesos.valores[i]<UMBRAL){
            n_menor_umbral++;
        }
    }

    return (100.0*n_menor_umbral)/(1.0*pesos.valores.size());
}

double fitness(double tasa_cla, double tasa_red){
    return ALPHA*tasa_cla + (1.0-ALPHA)*tasa_red;
}


//*******************************************************************************************************************************************************
//******************************************************************* RESULTADOS ************************************************************************
//*******************************************************************************************************************************************************

void resultadosP1(string caso){
    string nombre;
    if (caso!="BL" && caso!="RELIEF")
        caso = "1-NN";

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

            if(caso=="BL"){
                init = std::clock();
                W = BusquedaLocal(entrenamiento).first;         
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);

            } else if(caso=="RELIEF"){
                init = std::clock();
                W = GreedyRelief(entrenamiento).first;
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);

            } else {
                // Inicializando vector de pesos
                for(int j=0; j<W.valores.size();j++){
                    W.valores[j] = 1.0;
                }

                // Ahora se toman los estadísticos, midiendo el tiempo
                
                init = std::clock();
                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);
                fin = std::clock();
            }
            cout<<std::endl;
            for(int i=0;i<W.valores.size();i++){
                cout<<W.valores[i]<<" ";
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
