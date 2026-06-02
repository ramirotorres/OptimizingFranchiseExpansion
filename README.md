
# Associated data and code

The files in this code repository are updated and commented versions of the code used for the following work:

           ## Optimizing Franchise Expansion: Insights from Emerging Markets

## Main files in this repository
The code is written in C++ programming language interface. The ```Gurobi``` package is needed for the linear models and run the the C++ programs. 
The IPs formulation provide the optimal solution, meanwhile heuristics can lead to suboptimal decisions because they explore only a tiny percentage of possible solutions. 


* [`ReadMe.md`](ReadMe.md): This file.
  
C++ files:
* [`Modelo_IP.cpp`](Modelo_IP.cpp): MIP model.
* [`Modelo_Heuristica.cpp`](Modelo_Heuristica.cpp): Solution Approach (Heuristic Approach) .
* [`constantes.h`](constantes.h): Contains the Gurobi parameters and C++ libraries.

## Input data

This repository contains input data for IP and Heuristc Approach:

* [`coordenadas_brasil.txt`](coordenadas_brasil.txt): A geographical representation of potential locations for opening new franchises, based on their latitude and longitude coordinates. The total number of areas considered is $m=1,775$ and $T=5$. Moreover, each geographic area is complemented by the estimation of consumer demand ($D_ik$) during the period of the franchisee contract. This, each row represents:
id	Location_Name	latitude	longitude	Demand_year_1  Demand_year_2  Demand_year_3  Demand_year_4  Demand_year_5
 
* [`frontera_brasil.txt`](frontera_brasil.txt): The file contains the borders among areas (parameter $\delta_{ij}\in {0,1}$). If some oulet has no no border, the closest oulet is considered as a border. 

* [`Instance_Brasil.txt`](Instance_Brasil.txt): The file contains the following dat:

             m   T

             F   r   G

             R[1]   R[2],....,R[T]

             S[1]   S[2],....,S[T]

             latitude      longitude     \beta_i

* [Parameters](Instance_Brasil.txt): Thus, the parameters used in the aplication: $m=1775 , T=5, F=60000 , r=0.05  , G=1260000$ and vectors $R=(0.06 , 0.06 , 0.06 , 0.06 , 0.06)$, $S=(25200 , 25200 , 25200 , 25200 , 25200)$
