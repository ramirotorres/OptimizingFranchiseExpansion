
# Associated data and code

The files in this code repository are updated and commented versions of the code used for the following work:

            Optimizing Franchise Expansion: Insights from Emerging Markets

## Main files in this repository
The code is written using the C++ programming language. The Gurobi package is required to solve the linear models and execute the C++ programs. Integer Programming (IP) formulations provide optimal solutions, whereas heuristic methods may lead to suboptimal decisions because they explore only a small fraction of the possible solution space.


* [`ReadMe.md`](ReadMe.md): This file.
  
C++ files:
* [`Modelo_IP.cpp`](Modelo_IP.cpp): MIP model.
* [`Modelo_Heuristica.cpp`](Modelo_Heuristica.cpp): Solution Approach (Heuristic Approach) .
* [`constantes.h`](constantes.h): Contains the Gurobi parameters and C++ libraries.

## Input data

This repository contains the input data for the IP and Heuristic approaches.

* [`coordenadas_brasil.txt`](coordenadas_brasil.txt): A geographical representation of potential locations for opening new franchises is provided, based on their latitude and longitude coordinates. The total number of areas considered is $m = 1,775$, and the time horizon $T = 5$. Moreover, each geographical area includes an estimate of consumer demand ($D_{ik}$) during the franchise contract period. Specifically, each row represents:
 
          id	Location_Name	latitude	longitude	Demand_year_1  Demand_year_2  Demand_year_3  Demand_year_4  Demand_year_5
 
* [`frontera_brasil.txt`](frontera_brasil.txt): The file contains the borders between areas (parameter $\delta_{ij} \in {0,1}$). If an outlet has no neighboring areas, the closest outlet is considered a neighbor.

* [`Instance_Brasil.txt`](Instance_Brasil.txt): The file contains the following data:

             m   T

             F   r   G

             R[1]   R[2],....,R[T]

             S[1]   S[2],....,S[T]

             latitude      longitude     $\beta_i$

* [Parameters](Instance_Brasil.txt): Thus, the parameters used in this application are: $m=1775 , T=5, F=60000 , r=0.05  , G=1260000$ and vectors $R=(0.06 , 0.06 , 0.06 , 0.06 , 0.06)$, $S=(25200 , 25200 , 25200 , 25200 , 25200)$
