// Alumno: Álvaro Rodríguez Gallardo
// DNI: 77034155W
// Correo: alvaro155w@correo.ugr.es
// Grupo 2 (Miércoles, 17.30-19.30)
#include "aux.h"

using namespace std;


//*******************************************************************************************************************************************************
//******************************************************* PREPROCESAMIENTO DE DATOS *********************************************************************
//*******************************************************************************************************************************************************

// 'fichero' es el PATH al fichero .arff
// len(características)+1 == len(atributos)
Dataset lecturaFichero(const string fichero) {
    Dataset ds;

    ifstream archivo(fichero);
    string linea;
    bool leyendoDatos = false;
    if(!archivo.is_open()){
        cout << "No se ha podido abrir el archivo "<<fichero<<std::endl;
        return ds;
    }

    while(getline(archivo,linea)) {
        if(!leyendoDatos && linea.find("@attribute") != string::npos) {
            Atributo attr;
            stringstream ss(linea);
            string temp;
            ss >> temp; // Descarta "@attribute"
            ss >> attr.nombre >> attr.tipo;
            ds.atributos.push_back(attr);

        } else if(!leyendoDatos && linea.find("@data") != string::npos) {
            leyendoDatos = true; // Encontró @data, ahora comenzará a leer las muestras
            continue;
        } else if(leyendoDatos) {
            stringstream ss(linea);
            vector<string> vals;
            string valor;

            while(getline(ss, valor, ',')) {
                vals.push_back(valor);
            }

            string clase = vals.back();

            vals.pop_back();
            vector<double> numeros;
            for(auto& v : vals) {
                numeros.push_back(stod(v));
            }

            Muestra m;
            m.caracteristicas = numeros;
            m.clase = clase;

            ds.muestras.push_back(m);
        }
    }

    archivo.close();

    return ds;
}

void normalizarDatos(vector<Dataset> & vect_ds) {
    if (vect_ds.empty())
        return;
    vector<double> vals_min, vals_max;
    int num_caracts = vect_ds[0].muestras[0].caracteristicas.size();
    double min_def;
    double max_def;

    // Calculo máximos y mínimos GLOBALES por característica
    for(int i=0;i<num_caracts;i++){
        min_def = vect_ds[0].muestras[0].caracteristicas[i];
        max_def = vect_ds[0].muestras[0].caracteristicas[i];

        for(int j=0; j<vect_ds.size();j++){
            Dataset aux = vect_ds[j];
            for(int k=0;k<aux.muestras.size();k++){
                if(aux.muestras[k].caracteristicas[i] < min_def){
                    min_def = aux.muestras[k].caracteristicas[i];
                }
                if(aux.muestras[k].caracteristicas[i] > max_def){
                    max_def = aux.muestras[k].caracteristicas[i];
                }
            }
        }
        vals_min.push_back(min_def);
        vals_max.push_back(max_def);
    }
    
    for(int i=0;i<num_caracts;i++){
        for(int j=0;j<vect_ds.size();j++){
            for(int k=0;k<vect_ds[j].muestras.size();k++){
                vect_ds[j].muestras[k].caracteristicas[i] = (vect_ds[j].muestras[k].caracteristicas[i] - vals_min[i]) / (vals_max[i]-vals_min[i]);
            }
        }
    }
}


Dataset unirDatasets(const vector<Dataset>& vect_ds) {
    Dataset ds;

    if (!vect_ds.empty()) {
        // Asumiendo que los atributos son los mismos en todas las particiones y vect_ds no está vacío
        ds.atributos = vect_ds[0].atributos;

        for (const auto& aux : vect_ds) {
            ds.muestras.insert(ds.muestras.end(), aux.muestras.begin(), aux.muestras.end());
        }
    }
  //  mostrarDataset(ds);
    return ds;
}


//*******************************************************************************************************************************************************
//************************************************************* DISTANCIAS ******************************************************************************
//*******************************************************************************************************************************************************

double distanciaEuclidea(const Muestra & m1, const Muestra & m2) {
    if(m1.caracteristicas.size()!=m2.caracteristicas.size()){
        cout<<"Las dos muestras no son iguales" << std::endl;
        return -1;
    }

    double sum=0.0;

    for(int i=0; i<m1.caracteristicas.size();i++) {
        sum+=(m1.caracteristicas[i]-m2.caracteristicas[i])*(m1.caracteristicas[i]-m2.caracteristicas[i]);
    }

    return std::sqrt(sum);
}

double distanciaEuclideaPonderada(const Muestra & m1, const Muestra & m2, const Pesos & pesos) {
    if(m1.caracteristicas.size()!=m2.caracteristicas.size()){
        cout<<"Las dos muestras no son iguales" << std::endl;
        return -1;
    }

    if(m1.caracteristicas.size()!=pesos.valores.size()){
        cout<<"El vector de pesos no tiene mismo tamaño que las muestras"<<std::endl;
        return -1;
    }

    double sum=0.0;

    for(int i=0; i<m1.caracteristicas.size();i++) {
        if(pesos.valores[i] > 0.1){
            sum+=pesos.valores[i] * (m1.caracteristicas[i] - m2.caracteristicas[i]) * (m1.caracteristicas[i] - m2.caracteristicas[i]);
        }
    }

    return std::sqrt(sum);
}

int buscarMuestraAMenorDistancia(const Dataset & ds, const Muestra & m){
    int index_buscado = 0;
    double distancia_m_buscado = distanciaEuclidea(ds.muestras[index_buscado],m);

    for(int i=1; i<ds.muestras.size();i++){
        double distancia_m_probable = distanciaEuclidea(ds.muestras[i],m);

        if(distancia_m_probable < distancia_m_buscado){
            index_buscado = i;
            distancia_m_buscado = distancia_m_probable;
        }
    }

    return index_buscado;
}

//*******************************************************************************************************************************************************
//******************************************************************* AUXILIARES ************************************************************************
//*******************************************************************************************************************************************************

void mostrarDataset(const Dataset & ds){
    cout<<std::endl;
    cout<<"Los atributos del dataset y su tipo son:";
    cout<<std::endl;
    for(int i=0;i<ds.atributos.size();i++){
        cout<<"("<<ds.atributos[i].nombre<<","<<ds.atributos[i].tipo<<"), ";
    }
    cout<<"Los valores de las muestras son, junto a su clase, los siguientes";
    cout<<std::endl;
    for(int i=0;i<ds.muestras.size();i++){
        cout<<"Muestra "<<i+1<<" de tamaño "<<ds.muestras[i].caracteristicas.size()<<":"<<std::endl;
        for(int j=0;j<ds.muestras[i].caracteristicas.size();j++){
            cout<<ds.muestras[i].caracteristicas[j]<<" ";
            if(j==ds.muestras[i].caracteristicas.size()-1){
                cout<<std::endl;
                cout<<"Su clase es "<<ds.muestras[i].clase;
                cout<<std::endl;
            }
        }
    }
    cout<<std::endl;
    cout<<std::endl;
}

void escribirCSV(double tasa_cla,double tasa_red,double fitness,double tiempo,bool media){
    // Abre el archivo 'resultados.csv' en modo de escritura al final del archivo
    std::ofstream archivo("resultados.csv", std::ios::app);

    // Verifica si el archivo se ha abierto correctamente
    if (!archivo.is_open()) {
        std::cerr << "Error al abrir el archivo!" << std::endl;
        return;
    }

    // Configura el archivo para usar punto como separador decimal
    archivo.imbue(std::locale("en_US.UTF-8"));

    // Escribir los datos en el archivo con dos decimales y en notación científica si es necesario
    // Función lambda para escribir un valor notación científica si fuere necesario
    auto escribirValor = [&](double valor) {
        if (std::abs(valor) < 1e-3 && valor != 0.0) {
            archivo << std::scientific << std::setprecision(2) << valor << std::fixed;
        } else {
            archivo << std::fixed << std::setprecision(2) << valor;
        }
    };
    escribirValor(tasa_cla); archivo << ';';
    escribirValor(tasa_red); archivo << ';';
    escribirValor(fitness); archivo << ';';
    escribirValor(tiempo);

    // Si media es true, añade una marca o indicación de que los datos son promedios
    if (media) {
        archivo << ";medias";
    }

    archivo << std::endl;  // Inserta un salto de línea al final de cada registro

    // Cierra el archivo
    archivo.close();
}

void escribirFitnessParaConvergencia(int num_part, double fitness, int n_generacion, string algoritmo, string nom_carpeta){
   /* std::string filename = "fitness_Partic_" + std::to_string(num_part) + "_algoritmo_" + algoritmo + "_dataset_" + nombre_dataset + ".csv";
    
    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::app);
    
    if (!file.is_open()) {
        std::cerr << "Error al abrir el archivo: " << filename << std::endl;
        return;
    }
    
    // Comprobar si el archivo está vacío para escribir los encabezados
    file.seekp(0, std::ios::end); // Mover al final del archivo
    if (file.tellp() == 0) {
        // Escribir encabezados de las columnas si el archivo está vacío
        file << "Generación,Fitness\n";
    }
    
    // Escribir los datos proporcionados en el archivo
    file << n_generacion << "," << fitness << "\n";
    
    file.close();*/

    std::string filename = nom_carpeta + "/fitness_FuncID_" + std::to_string(num_part) + "_algoritmo_" + algoritmo + ".csv";

    //std::filesystem::create_directories(nom_carpeta);

    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::app);

    if (!file.is_open()) {
        std::cerr << "Error al abrir el archivo: " << filename << std::endl;
        return;
    }

    // Comprobar si el archivo está vacío para escribir los encabezados
    file.seekp(0, std::ios::end); // Mover al final del archivo
    if (file.tellp() == 0) {
        // Escribir encabezados de las columnas si el archivo está vacío
        file << "Generación,Fitness\n";
    }

    // Escribir los datos proporcionados en el archivo
    file << n_generacion << "," << fitness << "\n";

    file.close();
}