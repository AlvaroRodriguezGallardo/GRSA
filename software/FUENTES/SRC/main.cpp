// Alumno: Álvaro Rodríguez Gallardo
// DNI: 77034155W
// Correo: alvaro155w@correo.ugr.es
// Grupo 2 (Miércoles, 17.30-19.30)
#include <iostream>
#include <vector>
#include "practicaAlternativa.hpp"

using namespace std;

// ******************************************************************************************************************************************************
// *************************************************** FUNCIÓN PRINCIPAL ********************************************************************************
// ******************************************************************************************************************************************************

int main(int argc, char **argv){
    long int semilla;

    if(argc <= 1){
        cout << "Introduzca una semilla"<<endl;
        return -1;
    } else {
        semilla = atoi(argv[1]);
        Random::seed(semilla);
        cout <<"Semilla usada: " <<semilla <<endl;
    
        mostrarResultadosPAlternativa("GRSA-Todo");
        mostrarResultadosPAlternativa("GRSA-Mej");
        mostrarResultadosPAlternativa("GRSA-BL");
        mostrarResultadosPAlternativa("GRSA");

    }

    return 0;
}

