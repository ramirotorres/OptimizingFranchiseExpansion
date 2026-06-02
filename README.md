
# Associated data and code

The files in this code repository are updated and commented versions of the code used for the following work:

           Optimizing Franchise Expansion: Insights from Emerging Markets

## Main files in this repository
The code is written in C++ programming language interface. The ```Gurobi``` package is needed for the linear models and run the the C++ programs. 
The IPs formulation provide the optimal solution, meanwhile heuristics can lead to suboptimal decisions because they explore only a tiny percentage of possible solutions. 


* [`ReadMe.md`](ReadMe.md): This file.
  
C++ files:
* [`Modelo_IP.cpp`](Modelo_IP.cpp): MIP model (M<sub>1</sub>).
* [`Modelo_Heuristica.cpp`](Modelo_IP.cpp): Solution Approach (Heuristic Approach) (M<sub>1</sub>).
* [`constantes.h`](constantes.h): Contains the Gurobi parameters and C++ libraries.

This repository contains input data for IP and Hueristc Approach:

* [`coordenadas_brasil.txt`](Coordenates): A geographical representation of potential locations for opening new franchises, based on their latitude and longitude coordinates. The total number of areas considered is 1,775. Moreover, each geographic area is complemented by the estimation of consumer demand during the period of the franchisee contract. Each row represents:
id	Location_Name	latitude	longitude	Demand_year_1  Demand_year_2  Demand_year_3  Demand_year_4  Demand_year_5
  (M<sub>1</sub>).

* [`frontera_brasil.txt`](frontera_brasil.txt): The file contains the the locations with geographical border and their distance (M<sub>1</sub>) 

* [`Instance_Brasil.txt`](Instance_Brasil.txt): The file contains the following data (M<sub>2</sub>).

m T
F r G
R[1] R[2],....,R[T]
S[1] S[2],....,S[T]
latitude longitude \beta

Thus, the parameters used in the aplication
m=1775 
T=5
F=60000 
r=0.05 
G=1260000
R=(0.06 , 0.06 , 0.06 , 0.06 , 0.06)
S=(25200 , 25200 , 25200 , 25200 , 25200)
