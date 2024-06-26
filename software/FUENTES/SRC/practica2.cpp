// Alumno: Álvaro Rodríguez Gallardo
// DNI: 77034155W
// Correo: alvaro155w@correo.ugr.es
// Grupo 2 (Miércoles, 17.30-19.30)
#include "practica2.hpp"

//****************************************************** MÉTODOS CLASE CROMOSOMA *************************************************************************************

// 1 evaluación de fitness
Cromosoma::Cromosoma(const Dataset& entrenamiento,int N_pesos){
    std::uniform_real_distribution<double> distribution(0.0,1.0);

    for(int i=0;i<N_pesos;i++){
        this->pesos.valores.push_back(Random::get(distribution));
    }

    this->valor = fitness(tasaClasificacionLeaveOneOut(entrenamiento,this->pesos),tasaReduccion(this->pesos));
}

Cromosoma::Cromosoma(const Pesos& pesos, double valor){
    this->pesos = pesos;
    this->valor = valor;
}

bool Cromosoma::esMejorQue(const Cromosoma& crom) const{
    if(this->valor > crom.getValor()){
        return true;
    }

    return false;
}

Pesos Cromosoma::getPesos() const{
    return this->pesos;
}

double Cromosoma::getValor() const{
    return this->valor;
}

// 1 evaluación de fitness
void Cromosoma::mutacion(int indice,const Dataset& entrenamiento){
    std::normal_distribution<double> normal(0.0,std::sqrt(VARIANZA));

    this->pesos.valores[indice]+=Random::get(normal);

    if(this->pesos.valores[indice]>1.0){
        this->pesos.valores[indice] = 1.0;
    }
    if(this->pesos.valores[indice]<0.0){
        this->pesos.valores[indice] = 0.0;
    }
    this->valor = fitness(tasaClasificacionLeaveOneOut(entrenamiento,this->pesos),tasaReduccion(this->pesos));
}

//****************************************************** FUNCIONES AUXILIARES *************************************************************************************

void inicializacion(Poblacion & poblacion,const Dataset& entrenamiento){
    int N_pesos = entrenamiento.muestras[0].caracteristicas.size();

    for(int i=0; i<N_CROMOSOMAS; i++){
        Cromosoma crom(entrenamiento,N_pesos);
        poblacion.cromosomas.push_back(crom);
    }
}

//****************************************************** ALGORITMOS GENÉTICOS *************************************************************************************

// 2 evaluaciones de fitness
vector<Cromosoma> Cruce_BLX(const Cromosoma& c1, const Cromosoma& c2,const Dataset& entrenamiento){
    vector<Cromosoma> descendientes;
    double c_max, c_min;
    Pesos p_1 = c1.getPesos();
    Pesos p_2 = c2.getPesos();
    Pesos h1,h2;

    for(int i=0;i<p_1.valores.size();i++){
        c_min = min(p_1.valores[i],p_2.valores[i]);
        c_max = max(p_1.valores[i],p_2.valores[i]);

        double I = c_max-c_min;
        double minimo = c_min-I*ALPHA_BLX;
        double maximo = c_max+I*ALPHA_BLX;

        uniform_real_distribution<float> distrib(minimo, maximo);

        h1.valores.push_back(Random::get(distrib));
        h2.valores.push_back(Random::get(distrib));

        if(h1.valores[i] > 1)
            h1.valores[i] = 1.0;
        if(h1.valores[i] < 0)
            h1.valores[i] = 0.0;
        if(h2.valores[i] > 1)
            h2.valores[i] = 1.0;
        if(h2.valores[i] < 0)
            h2.valores[i] = 0.0;
    }

    double fitness1 = fitness(tasaClasificacionLeaveOneOut(entrenamiento,h1),tasaReduccion(h1));
    double fitness2 = fitness(tasaClasificacionLeaveOneOut(entrenamiento,h2),tasaReduccion(h2));

    descendientes.push_back(Cromosoma(h1,fitness1));
    descendientes.push_back(Cromosoma(h2,fitness2));

    assert(descendientes.size()==2);

    return descendientes;
}

// 2 evaluación de fitness
vector<Cromosoma> Cruce_Arit(const Cromosoma& c1, const Cromosoma& c2, const Dataset& entrenamiento){
    Pesos v_mutar1,v_mutar2;
    double fitness_mutar1, fitness_mutar2;
    Pesos p_1 = c1.getPesos();
    Pesos p_2 = c2.getPesos();
    std::uniform_real_distribution<double> distr(0.0,1.0);
    double alpha_arit = Random::get(distr);

    for(int i=0;i<p_1.valores.size();i++){
        v_mutar1.valores.push_back(alpha_arit*p_1.valores[i] + (1-alpha_arit)*p_2.valores[i]);
        v_mutar2.valores.push_back(alpha_arit*p_2.valores[i] + (1-alpha_arit)*p_1.valores[i]);
        
        if(v_mutar1.valores[i] > 1)
            v_mutar1.valores[i] = 1.0;
        if(v_mutar1.valores[i] < 0)
            v_mutar1.valores[i] = 0.0;
        if(v_mutar2.valores[i] > 1)
            v_mutar2.valores[i] = 1.0;
        if(v_mutar2.valores[i] < 0)
            v_mutar2.valores[i] = 0.0;
    }
    fitness_mutar1 = fitness(tasaClasificacionLeaveOneOut(entrenamiento,v_mutar1),tasaReduccion(v_mutar1));
    fitness_mutar2 = fitness(tasaClasificacionLeaveOneOut(entrenamiento,v_mutar2),tasaReduccion(v_mutar2));

    Cromosoma crom1(v_mutar1,fitness_mutar1);
    Cromosoma crom2(v_mutar2,fitness_mutar2);

    vector<Cromosoma> hijos;
    hijos.push_back(crom1);
    hijos.push_back(crom2);

    assert(hijos.size()==2);
    return hijos;
}

Cromosoma torneo(const Poblacion& poblacion){
    const int K = 3;    // A mayor K, mayor número de individuos compiten y puede discriminar más
    const int TAM_POBL = poblacion.cromosomas.size();
    uniform_int_distribution<int> distrib_entera(0, TAM_POBL - 1);
    Cromosoma ganador = poblacion.cromosomas[Random::get(distrib_entera)];

    // Escoger tres índices aleatorios de 0 a poblacion.cromosomas.size()
    for(int i=0; i < K;i++){
        Cromosoma candidato = poblacion.cromosomas[Random::get(distrib_entera)];
        if(candidato.esMejorQue(ganador)){
            ganador = candidato;
        }
    }
  
    return ganador;
}


Poblacion reemplazoGeneracional(const Poblacion& poblacion, const Poblacion& hijos){
    Poblacion PoblReemplazo = hijos;
    std::pair<Cromosoma,int> mejorPadre = mejorCromosoma(poblacion);
    std::pair<Cromosoma,int> peorHijo = peorCromosoma(hijos);

    if(!estaDentro(hijos,mejorPadre.first) && mejorPadre.first.esMejorQue(peorHijo.first)){
        PoblReemplazo.cromosomas[peorHijo.second] = mejorPadre.first;
    }

    return PoblReemplazo;
}

void reemplazoEstacionario(Poblacion& poblacion, const vector<Cromosoma>& hijos){
    assert(hijos.size()==2);

    Cromosoma hijo1 = hijos[0];
    Cromosoma hijo2 = hijos[1];

    std::pair<Cromosoma,int> peorCrom1 = peorCromosoma(poblacion);
    std::pair<Cromosoma,int> peorCrom2 = std::make_pair(poblacion.cromosomas[0],0);

    if(peorCrom1.second == 0){
        peorCrom2 = {poblacion.cromosomas[1],1};
    }
    // Buscamos el segundo peor cromosoma de la población (distinto a peorCrom1)
    for(int i=0;i<N_CROMOSOMAS;i++){
        if(i != peorCrom1.second){
            if(peorCrom2.first.esMejorQue(poblacion.cromosomas[i])){
                peorCrom2.first = poblacion.cromosomas[i];
                peorCrom2.second = i;
            }
        }
    }

    std::vector<Cromosoma> croms = {peorCrom1.first, peorCrom2.first, hijo1, hijo2};

    std::sort(croms.begin(), croms.end(), [](const Cromosoma& a, const Cromosoma& b) {
        return a.esMejorQue(b);
    });


    // Ahora croms contiene solo los dos mejores cromosomas ordenados de mejor a peor
    std::sort(croms.begin(), croms.end(), [](const Cromosoma& a, const Cromosoma& b) {
        return a.getValor() > b.getValor(); // Ordena por getValor de mayor a menor
    });
 
    poblacion.cromosomas[peorCrom1.second] = croms[0];
    poblacion.cromosomas[peorCrom2.second] = croms[1];

}

std::pair<Cromosoma,int> mejorCromosoma(const Poblacion& poblacion){
    Cromosoma mejor = poblacion.cromosomas[0];
    int ind_mejor = 0;

    for(int i=1;i<poblacion.cromosomas.size();i++){
        if(poblacion.cromosomas[i].esMejorQue(mejor)){
            mejor = poblacion.cromosomas[i];
            ind_mejor = i;
        }
    }

    return std::make_pair(mejor,ind_mejor);
}

std::pair<Cromosoma,int> peorCromosoma(const Poblacion& poblacion){
    Cromosoma peor = poblacion.cromosomas[0];
    int ind_peor = 0;

    for(int i=0;i<poblacion.cromosomas.size();i++){
        if(peor.esMejorQue(poblacion.cromosomas[i])){
            peor = poblacion.cromosomas[i];
            ind_peor = i;
        }
    }

    return std::make_pair(peor,ind_peor);
}

bool estaDentro(const Poblacion & poblacion, const Cromosoma& crom){
    for(int i=0;i<poblacion.cromosomas.size();i++){
        if(crom.getPesos().valores == poblacion.cromosomas[i].getPesos().valores){
            return true;
        }
    }

    return false;
}

std::pair<Pesos,double> AGG_BLX(const Dataset& entrenamiento, string nombre, int n_part){
    Poblacion poblacion;
    int evaluaciones = 0;
    int num_generaciones = 1;

    inicializacion(poblacion,entrenamiento);
    evaluaciones+=poblacion.cromosomas.size();

    int num_cruces_esperanza_mat = std::round(PROB_CRUCE_AGG * (poblacion.cromosomas.size()/2));
    if(num_cruces_esperanza_mat<1){
        num_cruces_esperanza_mat = 1;
    }
    int numGenes = poblacion.cromosomas.size()*poblacion.cromosomas[0].getPesos().valores.size();
    int num_mutaciones_esp_mat = std::round(PROB_MUTAR_INDIVIDUO_GENETICOS*numGenes);
    
    if(num_mutaciones_esp_mat < 1){
        num_mutaciones_esp_mat = 1;
    }

    uniform_int_distribution<int> distrib_indices_uniforme(0,poblacion.cromosomas[0].getPesos().valores.size()-1);
    uniform_int_distribution<int> distrib_cromos_uniforme(0,poblacion.cromosomas.size()-1);

    while(evaluaciones < MAX_EVALUACIONES){
        Poblacion intermedio;
        if (nombre != BCANCER)
            escribirFitnessParaConvergencia(n_part,mejorCromosoma(poblacion).first.getValor(),num_generaciones,"AGG_BLX",nombre);

        // Torneo N veces en AGG
        for(int i=0;i<N_CROMOSOMAS;i++){
            intermedio.cromosomas.push_back(torneo(poblacion));
        }

        // Cruces según BL-0.3 a partir de cierta probabilidad para la cantidad de
        // cruces esperados
        // Se hace shuffle, pues si la población es parecida a la anterior, se va a quedar estancado
        Random::shuffle(intermedio.cromosomas);

        for(int i=0;i<2*num_cruces_esperanza_mat;i+=2){
            vector<Cromosoma> desc_BLX = Cruce_BLX(intermedio.cromosomas[i],intermedio.cromosomas[i+1],entrenamiento);
            intermedio.cromosomas[i] = desc_BLX[0];
            intermedio.cromosomas[i+1] = desc_BLX[1];
            evaluaciones+=2;
        }


        // Mutación de población intermedia para cromosomas hijos
        // 1- Seleccionar un cromosoma al azar de intermedio
        // 2- Mutar una posición aleatoria
        // En cada mutación se ha hecho evaluación de fitness
        for(int i=0;i<num_mutaciones_esp_mat;i++){      
            int index_cromosoma_mutar = Random::get(distrib_cromos_uniforme);
            int index_gen_mutar = Random::get(distrib_indices_uniforme);
            intermedio.cromosomas[index_cromosoma_mutar].mutacion(index_gen_mutar,entrenamiento);
            evaluaciones++;
        }

        // Reemplazo con elitismo y reasignación
        poblacion = reemplazoGeneracional(poblacion,intermedio);
    //    cout<<"Tamaño población: "<<poblacion.cromosomas.size()<<std::endl;
        num_generaciones++;
    }

    Cromosoma mejorCrom = mejorCromosoma(poblacion).first;

    return std::make_pair(mejorCrom.getPesos(),mejorCrom.getValor()); 
}

std::pair<Pesos,double> AGG_Arit(const Dataset& entrenamiento, string nombre, int n_part){
    Poblacion poblacion;
    int evaluaciones = 0;
    int num_generaciones = 1;

    inicializacion(poblacion,entrenamiento);

    int num_cruces_esperanza_mat = std::round(PROB_CRUCE_AGG * (poblacion.cromosomas.size()/2));
    int numGenes = poblacion.cromosomas.size()*poblacion.cromosomas[0].getPesos().valores.size();
    int num_mutaciones_esp_mat = std::round(PROB_MUTAR_INDIVIDUO_GENETICOS*numGenes);
    uniform_int_distribution<int> distrib_indices_uniforme(0,poblacion.cromosomas[0].getPesos().valores.size()-1);
    uniform_int_distribution<int> distrib_cromos_uniforme(0,poblacion.cromosomas.size()-1);

    evaluaciones+=poblacion.cromosomas.size();

    while(evaluaciones < MAX_EVALUACIONES){
        Poblacion intermedio;
        if (nombre != BCANCER)
            escribirFitnessParaConvergencia(n_part,mejorCromosoma(poblacion).first.getValor(),num_generaciones,"AGG_Arit",nombre);

        for(int i=0;i<N_CROMOSOMAS;i++){
            intermedio.cromosomas.push_back(torneo(poblacion));
        }

        // Cruces según BL-Arit a partir de cierta probabilidad para la cantidad de
        // cruces esperados
        // Se hace shuffle, pues si la población es parecida a la anterior, los hijos serán parecidos
        // de los cruces pueden ser parecidos a los de la población y no dar heterogeneidad
        Random::shuffle(intermedio.cromosomas);
       // cout<<"El número esperado de cruces es "<<2*num_cruces_esperanza_mat<<std::endl;
        for(int i=0;i<2*num_cruces_esperanza_mat;i+=2){
            vector<Cromosoma> desc_Arit = Cruce_Arit(intermedio.cromosomas[i],intermedio.cromosomas[N_CROMOSOMAS-i-1],entrenamiento);
            intermedio.cromosomas[i] = desc_Arit[0];
            intermedio.cromosomas[i+1] = desc_Arit[1];
            evaluaciones+=2;
        }
        //cout<<"Tras el cruce intermedio tiene tamaño de "<<intermedio.cromosomas.size()<<std::endl;

        // Mutación de población intermedia para cromosomas hijos
        // 1- Seleccionar un cromosoma al azar de intermedio
        // 2- Mutar una posición aleatoria
        // En cada mutación se ha hecho evaluación de fitness
        //cout<<"El número esperado de mutaciones es "<<num_mutaciones_esp_mat<<std::endl;
        for(int i=0;i<num_mutaciones_esp_mat;i++){
            int index_cromosoma_mutar = Random::get(distrib_cromos_uniforme);
            int index_gen_mutar = Random::get(distrib_indices_uniforme);

            intermedio.cromosomas[index_cromosoma_mutar].mutacion(index_gen_mutar,entrenamiento);
            evaluaciones++;
        }
        //cout<<"Hasta ahora hay estas evaluaciones "<<evaluaciones<<std::endl;
        // Reemplazo con elitismo y reasignación
        poblacion = reemplazoGeneracional(poblacion,intermedio);
   
        num_generaciones++;
    }

    Cromosoma mejorCrom = mejorCromosoma(poblacion).first;

    return std::make_pair(mejorCrom.getPesos(),mejorCrom.getValor()); 
}

std::pair<Pesos,double> AGE_BLX(const Dataset& entrenamiento, string nombre, int n_part){
    Poblacion poblacion;
    int evaluaciones = 0;
    int num_generaciones = 1;
    const int MAX_GANADORES = 2;
    inicializacion(poblacion,entrenamiento);
    evaluaciones+=poblacion.cromosomas.size();

    //int num_cruces_esperanza_mat = std::round(PROB_CRUCE_AGE * (poblacion.cromosomas.size()/2));
    //int numGenes = poblacion.cromosomas.size()*poblacion.cromosomas[0].getPesos().size();
    //int num_mutaciones_esp_mat = std::round(PROB_MUTAR_INDIVIDUO_GENETICOS*numGenes);
    uniform_int_distribution<int> distrib_indices_uniforme(0,poblacion.cromosomas[0].getPesos().valores.size()-1);
   // uniform_int_distribution<int> distrib_cromos_uniforme(0,poblacion.cromosomas.size());
    uniform_int_distribution<int> distr_para_cruce(0,1);

    while(evaluaciones < MAX_EVALUACIONES){
        Poblacion intermedio;
        if (nombre != BCANCER && ((num_generaciones-1) % 15 == 0))
            escribirFitnessParaConvergencia(n_part,mejorCromosoma(poblacion).first.getValor(),num_generaciones,"AGE_BLX",nombre);

        // Selección, ya se han hecho las evaluaciones correspondientes
        for(int i=0;i<MAX_GANADORES;i++){
            intermedio.cromosomas.push_back(torneo(poblacion));
        }
        // Cruce, dos evaluaciones de fitness por llamada
        vector<Cromosoma> cruces = Cruce_BLX(intermedio.cromosomas[0],intermedio.cromosomas[1],entrenamiento);
        evaluaciones+=2;

        int numGenes = cruces.size()*cruces[0].getPesos().valores.size();
        int num_mutaciones_esp_mat = std::round(PROB_MUTAR_INDIVIDUO_GENETICOS*numGenes);

        Random::shuffle(cruces);

        for(int i=0;i<num_mutaciones_esp_mat;i++){
            cruces[Random::get(distr_para_cruce)].mutacion(Random::get(distrib_indices_uniforme),entrenamiento);
            evaluaciones++;
        }

        // Reemplazar según AGE, ambos hijos luchan por incorporarse a la población
        reemplazoEstacionario(poblacion,cruces);
        num_generaciones++;
    }

    Cromosoma mejorCrom = mejorCromosoma(poblacion).first;

    return std::make_pair(mejorCrom.getPesos(),mejorCrom.getValor());
}

std::pair<Pesos,double> AGE_Arit(const Dataset& entrenamiento, string nombre, int n_part){
    Poblacion poblacion;
    int evaluaciones = 0;
    int num_generaciones = 1;
    const int MAX_GANADORES = 2;
    inicializacion(poblacion,entrenamiento);
    evaluaciones+=poblacion.cromosomas.size();

    //int num_cruces_esperanza_mat = std::round(PROB_CRUCE_AGE * (poblacion.cromosomas.size()/2));
    //int numGenes = poblacion.cromosomas.size()*poblacion.cromosomas[0].getPesos().size();
    //int num_mutaciones_esp_mat = std::round(PROB_MUTAR_INDIVIDUO_GENETICOS*numGenes);
    uniform_int_distribution<int> distrib_indices_uniforme(0,poblacion.cromosomas[0].getPesos().valores.size()-1);
   // uniform_int_distribution<int> distrib_cromos_uniforme(0,poblacion.cromosomas.size());
    uniform_int_distribution<int> distr_para_cruce(0,1);
    
    while(evaluaciones < MAX_EVALUACIONES){
        Poblacion intermedio;
        if(nombre != BCANCER  && ((num_generaciones-1) % 15 == 0))
            escribirFitnessParaConvergencia(n_part,mejorCromosoma(poblacion).first.getValor(),num_generaciones,"AGE_Arit",nombre);

        // Selección, ya se han hecho las evaluaciones correspondientes
        for(int i=0;i<MAX_GANADORES;i++){
            intermedio.cromosomas.push_back(torneo(poblacion));
        }

        // Cruce, dos evaluaciones de fitness por llamada
        vector<Cromosoma> cruces = Cruce_Arit(intermedio.cromosomas[0],intermedio.cromosomas[1],entrenamiento);
        evaluaciones+=2;

        int numGenes = cruces.size()*cruces[0].getPesos().valores.size();
        int num_mutaciones_esp_mat = std::round(PROB_MUTAR_INDIVIDUO_GENETICOS*numGenes);

        Random::shuffle(cruces);
        for(int i=0;i<num_mutaciones_esp_mat;i++){
            cruces[Random::get(distr_para_cruce)].mutacion(Random::get(distrib_indices_uniforme),entrenamiento);
            evaluaciones++;
        }

        // Reemplazar según AGE, ambos hijos luchan por incorporarse a la población
        reemplazoEstacionario(poblacion,cruces);
        num_generaciones++;
    }

    Cromosoma mejorCrom = mejorCromosoma(poblacion).first;

    return std::make_pair(mejorCrom.getPesos(),mejorCrom.getValor());
}

//****************************************************** ALGORITMOS MEMÉTICOS *************************************************************************************

vector<int> devolverSubconjuntoCromosomas(const Poblacion& poblacion, int tam_subconjunto, bool elMejor){
    vector<int> subconjunto;

    if(elMejor){
        const std::vector<Cromosoma>& misCromos = poblacion.cromosomas;

        std::vector<std::pair<Cromosoma, int>> cromosomasConIndices;
        for (size_t i = 0; i < N_CROMOSOMAS; ++i) {
            cromosomasConIndices.push_back({misCromos[i], i});
        }

        // Ordenar los cromosomas por getValor de mayor a menor
        std::sort(cromosomasConIndices.begin(), cromosomasConIndices.end(), 
            [](const std::pair<Cromosoma, int>& a, const std::pair<Cromosoma, int>& b) {
                return a.first.getValor() > b.first.getValor();
            });

        for(int i=0;i<tam_subconjunto;i++){
            subconjunto.push_back(cromosomasConIndices[i].second);
        }

    } else {
        std::random_device rd;
        std::mt19937 gen(rd());
        // Como se escoge o no según probabilidad de 0.1, la distribución de Bernoulli
        // es la que mejor modela el escenario
        std::bernoulli_distribution distrib_bern(PROB_MEMETICOS_BL);

        for(size_t i=0;i<N_CROMOSOMAS;i++){
            if(distrib_bern(gen)){
                subconjunto.push_back(i);
            }
        }

    }

    return subconjunto;

}

int busquedaLocalPractica2(const Dataset& entrenamiento, Cromosoma& crom){
    Pesos v_pesos = crom.getPesos();
    int evaluaciones = 0;
    int n_vecinos_evaluados = 0;
    const int TOPE = LIMITE_BL*v_pesos.valores.size();
    const int MAX_ITERACS = 15000;
    int num_caracteristicas = entrenamiento.muestras[0].caracteristicas.size();
    vector<int> indexes(num_caracteristicas);
    double tasa_cla_vecino,tasa_red_vecino,fitness_vecino;
    //bool hubo_mejora = false;
    std::normal_distribution<double> normal(0.0,std::sqrt(VARIANZA));
    double fitness_actual = crom.getValor();

    mezclarIndices(num_caracteristicas,indexes,true);

    while(n_vecinos_evaluados < LIMITE * num_caracteristicas){
        // Cada vez que entra al bucle, se mezclan los índices. Así, se añade aleatoriedad en la generación, evitando repetir índice
     //   hubo_mejora = false;

        for(int j=0;j<indexes.size()/* && !hubo_mejora*/;j++){
            Pesos vecino = v_pesos;
            mutacionBL(vecino,indexes[j],normal);

            tasa_cla_vecino = tasaClasificacionLeaveOneOut(entrenamiento,vecino);
            tasa_red_vecino = tasaReduccion(vecino);
            fitness_vecino = fitness(tasa_cla_vecino,tasa_red_vecino);
            
            evaluaciones++;
            n_vecinos_evaluados++;

            if(fitness_vecino > fitness_actual){
                fitness_actual = fitness_vecino;
                v_pesos = vecino;
            //    hubo_mejora = true;
            }
        }

        if (/*hubo_mejora || */n_vecinos_evaluados % num_caracteristicas==0){
            mezclarIndices(num_caracteristicas,indexes,false);
        //    hubo_mejora = false;
        }

    }

    crom = Cromosoma(v_pesos,fitness_actual);
    return evaluaciones; 

}

std::pair<Pesos,double> AM_All(const Dataset& entrenamiento, string nombre, int n_part){
    Poblacion poblacion;
    int evaluaciones = 0;
    int num_generaciones = 1;

    inicializacion(poblacion,entrenamiento);

    int num_cruces_esperanza_mat = std::round(PROB_CRUCE_MEMETICOS * (poblacion.cromosomas.size()/2));
    int numGenes = poblacion.cromosomas.size()*poblacion.cromosomas[0].getPesos().valores.size();
    int num_mutaciones_esp_mat = std::round(PROB_MUTAR_INDIVIDUO_MEMETICOS*numGenes);
    uniform_int_distribution<int> distrib_indices_uniforme(0,poblacion.cromosomas[0].getPesos().valores.size()-1);
    uniform_int_distribution<int> distrib_cromos_uniforme(0,poblacion.cromosomas.size()-1);

    evaluaciones+=poblacion.cromosomas.size();

    while(evaluaciones < MAX_EVALUACIONES){
        Poblacion intermedio;
        if(nombre != BCANCER)
            escribirFitnessParaConvergencia(n_part,mejorCromosoma(poblacion).first.getValor(),num_generaciones,"AM_All",nombre);

        for(int i=0;i<N_CROMOSOMAS;i++){
            intermedio.cromosomas.push_back(torneo(poblacion));
        }

        // Cruces según BL-Arit a partir de cierta probabilidad para la cantidad de
        // cruces esperados
        // Se hace shuffle, pues si la población es parecida a la anterior, los hijos serán parecidos
        // de los cruces pueden ser parecidos a los de la población y no dar heterogeneidad
        Random::shuffle(intermedio.cromosomas);
       // cout<<"El número esperado de cruces es "<<2*num_cruces_esperanza_mat<<std::endl;
        for(int i=0;i<2*num_cruces_esperanza_mat;i+=2){
            vector<Cromosoma> desc_Arit = Cruce_Arit(intermedio.cromosomas[i],intermedio.cromosomas[N_CROMOSOMAS-i-1],entrenamiento);
            intermedio.cromosomas[i] = desc_Arit[0];
            intermedio.cromosomas[i+1] = desc_Arit[1];
            evaluaciones+=2;
        }
        //cout<<"Tras el cruce intermedio tiene tamaño de "<<intermedio.cromosomas.size()<<std::endl;

        // Mutación de población intermedia para cromosomas hijos
        // 1- Seleccionar un cromosoma al azar de intermedio
        // 2- Mutar una posición aleatoria
        // En cada mutación se ha hecho evaluación de fitness
        //cout<<"El número esperado de mutaciones es "<<num_mutaciones_esp_mat<<std::endl;
        for(int i=0;i<num_mutaciones_esp_mat;i++){
            int index_cromosoma_mutar = Random::get(distrib_cromos_uniforme);
            int index_gen_mutar = Random::get(distrib_indices_uniforme);

            intermedio.cromosomas[index_cromosoma_mutar].mutacion(index_gen_mutar,entrenamiento);
            evaluaciones++;
        }
        //cout<<"Hasta ahora hay estas evaluaciones "<<evaluaciones<<std::endl;
        // Reemplazo con elitismo y reasignación
        poblacion = reemplazoGeneracional(poblacion,intermedio);
   
        num_generaciones++;
    }
       

    Cromosoma mejorCrom = mejorCromosoma(poblacion).first;

    return std::make_pair(mejorCrom.getPesos(),mejorCrom.getValor());
}

std::pair<Pesos,double> AM_Rand(const Dataset& entrenamiento, string nombre, int n_part){
    Poblacion poblacion;
    int evaluaciones = 0;
    int num_generaciones = 1;

    inicializacion(poblacion,entrenamiento);

    int num_cruces_esperanza_mat = std::round(PROB_CRUCE_MEMETICOS * (poblacion.cromosomas.size()/2));
    int numGenes = poblacion.cromosomas.size()*poblacion.cromosomas[0].getPesos().valores.size();
    int num_mutaciones_esp_mat = std::round(PROB_MUTAR_INDIVIDUO_MEMETICOS*numGenes);
    uniform_int_distribution<int> distrib_indices_uniforme(0,poblacion.cromosomas[0].getPesos().valores.size()-1);
    uniform_int_distribution<int> distrib_cromos_uniforme(0,poblacion.cromosomas.size()-1);

    evaluaciones+=poblacion.cromosomas.size();

    while(evaluaciones < MAX_EVALUACIONES){
        Poblacion intermedio;
        if(nombre != BCANCER)
            escribirFitnessParaConvergencia(n_part,mejorCromosoma(poblacion).first.getValor(),num_generaciones,"AM_Rand",nombre);

        for(int i=0;i<N_CROMOSOMAS;i++){
            intermedio.cromosomas.push_back(torneo(poblacion));
        }

        // Cruces según BL-Arit a partir de cierta probabilidad para la cantidad de
        // cruces esperados
        // Se hace shuffle, pues si la población es parecida a la anterior, los hijos serán parecidos
        // de los cruces pueden ser parecidos a los de la población y no dar heterogeneidad
        Random::shuffle(intermedio.cromosomas);
       // cout<<"El número esperado de cruces es "<<2*num_cruces_esperanza_mat<<std::endl;
        for(int i=0;i<2*num_cruces_esperanza_mat;i+=2){
            vector<Cromosoma> desc_Arit = Cruce_Arit(intermedio.cromosomas[i],intermedio.cromosomas[N_CROMOSOMAS-i-1],entrenamiento);
            intermedio.cromosomas[i] = desc_Arit[0];
            intermedio.cromosomas[i+1] = desc_Arit[1];
            evaluaciones+=2;
        }
        //cout<<"Tras el cruce intermedio tiene tamaño de "<<intermedio.cromosomas.size()<<std::endl;

        // Mutación de población intermedia para cromosomas hijos
        // 1- Seleccionar un cromosoma al azar de intermedio
        // 2- Mutar una posición aleatoria
        // En cada mutación se ha hecho evaluación de fitness
        //cout<<"El número esperado de mutaciones es "<<num_mutaciones_esp_mat<<std::endl;
        for(int i=0;i<num_mutaciones_esp_mat;i++){
            int index_cromosoma_mutar = Random::get(distrib_cromos_uniforme);
            int index_gen_mutar = Random::get(distrib_indices_uniforme);

            intermedio.cromosomas[index_cromosoma_mutar].mutacion(index_gen_mutar,entrenamiento);
            evaluaciones++;
        }
        //cout<<"Hasta ahora hay estas evaluaciones "<<evaluaciones<<std::endl;
        // Reemplazo con elitismo y reasignación
        poblacion = reemplazoGeneracional(poblacion,intermedio);
   
        num_generaciones++;
    }

    Cromosoma mejorCrom = mejorCromosoma(poblacion).first;

    return std::make_pair(mejorCrom.getPesos(),mejorCrom.getValor());
}

std::pair<Pesos,double> AM_Best(const Dataset& entrenamiento, string nombre, int n_part){
    Poblacion poblacion;
    int evaluaciones = 0;
    int num_generaciones = 1;

    inicializacion(poblacion,entrenamiento);

    int num_cruces_esperanza_mat = std::round(PROB_CRUCE_MEMETICOS * (poblacion.cromosomas.size()/2));
    int numGenes = poblacion.cromosomas.size()*poblacion.cromosomas[0].getPesos().valores.size();
    int num_mutaciones_esp_mat = std::round(PROB_MUTAR_INDIVIDUO_MEMETICOS*numGenes);
    uniform_int_distribution<int> distrib_indices_uniforme(0,poblacion.cromosomas[0].getPesos().valores.size()-1);
    uniform_int_distribution<int> distrib_cromos_uniforme(0,poblacion.cromosomas.size()-1);

    evaluaciones+=poblacion.cromosomas.size();

    while(evaluaciones < MAX_EVALUACIONES){
        Poblacion intermedio;
        if(nombre != BCANCER)
            escribirFitnessParaConvergencia(n_part,mejorCromosoma(poblacion).first.getValor(),num_generaciones,"AM_Best",nombre);

        for(int i=0;i<N_CROMOSOMAS;i++){
            intermedio.cromosomas.push_back(torneo(poblacion));
        }

        // Cruces según BL-Arit a partir de cierta probabilidad para la cantidad de
        // cruces esperados
        // Se hace shuffle, pues si la población es parecida a la anterior, los hijos serán parecidos
        // de los cruces pueden ser parecidos a los de la población y no dar heterogeneidad
        Random::shuffle(intermedio.cromosomas);
       // cout<<"El número esperado de cruces es "<<2*num_cruces_esperanza_mat<<std::endl;
        for(int i=0;i<2*num_cruces_esperanza_mat;i+=2){
            vector<Cromosoma> desc_Arit = Cruce_Arit(intermedio.cromosomas[i],intermedio.cromosomas[N_CROMOSOMAS-i-1],entrenamiento);
            intermedio.cromosomas[i] = desc_Arit[0];
            intermedio.cromosomas[i+1] = desc_Arit[1];
            evaluaciones+=2;
        }
        //cout<<"Tras el cruce intermedio tiene tamaño de "<<intermedio.cromosomas.size()<<std::endl;

        // Mutación de población intermedia para cromosomas hijos
        // 1- Seleccionar un cromosoma al azar de intermedio
        // 2- Mutar una posición aleatoria
        // En cada mutación se ha hecho evaluación de fitness
        //cout<<"El número esperado de mutaciones es "<<num_mutaciones_esp_mat<<std::endl;
        for(int i=0;i<num_mutaciones_esp_mat;i++){
            int index_cromosoma_mutar = Random::get(distrib_cromos_uniforme);
            int index_gen_mutar = Random::get(distrib_indices_uniforme);

            intermedio.cromosomas[index_cromosoma_mutar].mutacion(index_gen_mutar,entrenamiento);
            evaluaciones++;
        }
        //cout<<"Hasta ahora hay estas evaluaciones "<<evaluaciones<<std::endl;
        // Reemplazo con elitismo y reasignación
        poblacion = reemplazoGeneracional(poblacion,intermedio);
   
        num_generaciones++;
    }

    Cromosoma mejorCrom = mejorCromosoma(poblacion).first;

    return std::make_pair(mejorCrom.getPesos(),mejorCrom.getValor());
}

//************************************************************* RESULTADOS *************************************************************************************

void resultadosP2(const string& caso){
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

            if(caso=="AGG-BX"){
                init = std::clock();
                W = AGG_BLX(entrenamiento,nombre,i+1).first;         
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);

            } else if(caso=="AGG-CA"){
                init = std::clock();
                W = AGG_Arit(entrenamiento,nombre,i+1).first;
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);

            } else if(caso == "AGE-BX"){
                init = std::clock();
                W = AGE_BLX(entrenamiento,nombre,i+1).first;
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);

            } else if(caso == "AGE-CA"){
                init = std::clock();
                W = AGE_Arit(entrenamiento,nombre,i+1).first;
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);

            } else if(caso == "AM-1"){
                init = std::clock();
                W = AM_All(entrenamiento,nombre,i+1).first;
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);
            } else if(caso == "AM-01"){
                init = std::clock();
                W = AM_Rand(entrenamiento,nombre,i+1).first;
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);
            } else if(caso == "AM-mej"){
                init = std::clock();
                W = AM_Best(entrenamiento,nombre,i+1).first;
                fin = std::clock();

                tasa_cla_i = tasaClasificacion(entrenamiento,test,W);
                tasa_red_i = tasaReduccion(W);
                fitness_i = fitness(tasa_cla_i,tasa_red_i);
            } else {
                std::cerr<<"\nNO ES UN ALGORITMO DE LA PRACTICA 2"<<std::endl;
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