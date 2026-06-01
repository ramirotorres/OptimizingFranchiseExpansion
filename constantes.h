#ifndef CONSTANTES_H
#define CONSTANTES_H

#include <iostream>
#include <sstream>
#include <gurobi_c++.h>
#include <vector>
#include <list>
#include <queue>
#include <set>
#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <numeric>
using namespace std;

///PARAMETROS GUROBI

const double TIEMPO_IP = 1000;
const double GAP=0;

const string NOMBRE_INSTANCIA="instancia1.txt";
const string REPORTE="salidas/Reporte.csv";
const double INF = 1e+10;
///COLORES

const string RESET="\033[0m";
const string  ROJO="\033[31m";      /* Red */
const string VERDE="\033[32m";     /* Green */
const string  AZUL="\033[94m";     /* Blue */
const string  BLANCO="\033[37m";      /* White */


#endif // CONSTANTES_H
