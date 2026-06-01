#include "constantes.h"


// PARAMETERS OF THE MODELS
//AVERAGE DISTANCE
double d_prom=24.19;
double mu=3*d_prom;
//EFECT OF THE DISTANCE PN THE DEMAND
double alfa=0.01;
//MAXIMUM NUMBER OF OUTLETS 
int L=5;

template<class T> 
void escribir_matriz(vector<vector<T> > A)
{
    for(int i=0;i<A.size();i++)
    {
        for(int j=0;j<A[i].size();j++)
        {
            cout<<A[i][j]<<" ";
        }
        cout<<endl;
    }
}

class ciudad
{
private:
    int id;
    string Nombre;
    pair<double,double> coordenada;//latitud y longitud en polares
public:
    vector<int> demanda;
    list<int> frontera;
    int Bi;

ciudad(int k=0)
{
    demanda.resize(k);
    id=-1;
    Nombre="SN";
    coordenada=make_pair(0,0);
    Bi=0;
}

string Nombre_Ciudad()
{
    return Nombre;
}
int  Obt_id()
{
    return id;
}

double latitud()
{
    return coordenada.first;
}

double longitud()
{
    return coordenada.second;
}

friend ostream& operator<<(ostream& os, ciudad c)
{
    os<<"("<<c.id<<","<<c.Nombre<<","<<c.coordenada.first<<","<<c.coordenada.second<<")"
      <<endl<<"pob: ";

    for(int i=0;i<c.demanda.size();i++)
        os<<c.demanda[i]<<" , ";

    return os;
}

friend istream& operator>>(istream& is, ciudad& c)
{
    double x,y;
    is>>c.id>>c.Nombre>>x>>y;
    c.coordenada=make_pair(x,y);
    
    for(int i=0;i<c.demanda.size();i++)
    {
        is>>c.demanda[i];
        c.demanda[i]=24*c.demanda[i];
    }

    return is;
}
};

class Modelo
{
private:
  // Puntero al Solver
GRBEnv* env_modelo;
GRBModel*  modelo_ip;
    // Vector de (punteros a las) variables
    
  vector< vector<GRBVar> >  _x;//x_ij
  vector< GRBVar >  _y;//y_i
 
    vector<vector<int> > _Ady;
public:

  // Constructor y destructor
  Modelo();
  ~Modelo();

void incluir_parametros();
void crear_vars(vector<ciudad> C, double F,double r, vector<double> Rk, vector<double> Sk);
void crear_restricciones(vector<ciudad> C, double G,double F,double r, vector<double> Rk, vector<double> Sk);
ostream& presentar_sol(ostream& os, int &nl, vector<ciudad> C);

void write_MIP(string archivo)
{
    modelo_ip->write(archivo);
}
void resolver()
{
    modelo_ip->optimize();
}

double Lower()
{
    return modelo_ip->get(GRB_DoubleAttr_ObjBound);
}

int ObtN_Vars()
{
    return modelo_ip->get(GRB_IntAttr_NumVars);
}
int ObtN_Restr()
{
    return modelo_ip->get(GRB_IntAttr_NumConstrs);
}

double Obt_Obj()
{
    return modelo_ip->get(GRB_DoubleAttr_ObjVal);
}

double Obt_tiempo()
{
    return modelo_ip->get(GRB_DoubleAttr_Runtime);
}

double Obt_Brecha()
{
    return 100*modelo_ip->get(GRB_DoubleAttr_MIPGap);
}
int status()
{
    return modelo_ip->get(GRB_IntAttr_Status);
}
int Obt_N_Nodos()
{
    return modelo_ip->get(GRB_DoubleAttr_NodeCount);
}

};

double distancia_coord(ciudad c1,ciudad c2);
/* void leer_instancia(string ARCHIVO,int&m, int&t, double &F, double& r,
double& alfa, double& G,  
vector<double> &Rk,vector<double>& Sk,vector<int>& Bi,
vector<vector<int> > &Dik, vector<vector<double> > &dij,vector<vector<int> > &fij);
*/
int main()
{
    cout<<endl<<endl
        <<"***********************"<<endl
        <<"MODELO DE FRANQUICIAS: "<<endl
        <<"***********************"<<endl;   
//READING DATA
//CITIES 
vector<ciudad> Ciudades;
ifstream fa("coordenadas_brasil.txt");
if(fa.is_open()==true)
{
    int n,k;
    fa>>n>>k;
    Ciudades.resize(n);
    for(int i=0;i<n;i++)
    {
        ciudad c(k);
        fa>>c;
        Ciudades[i]=c;
    }
}
else
{
    cout<<"Problemas de lectura"<<endl;
}
fa.close();

//Se ingresan como áreas de frontera si son menores que el promedio
//AREAS ARE CONSIDERED IN THE BORDER IF THE DISTANCE IS LESS THAN THE AVERAGE
cout<<"... Calculando fronteras"<<endl;
for(int ci=0;ci<Ciudades.size();ci++)
{
    for(int cj=ci+1;cj<Ciudades.size();cj++)
    {
        if(distancia_coord(Ciudades[ci],Ciudades[cj])<=mu)
        {
            Ciudades[ci].frontera.push_back(cj);
            Ciudades[cj].frontera.push_back(ci);
         }
    }
}
//READING THE FILE WITH BORDERS AREAS. 
//BORDER AREAS INCLUDED IN THE PREVIOUS PHASE ARE NOT CONSIDERED  
//leyendo archivo de fronteras
//Del archivo de frontera se ingresan aquellas que la distancia es mayor al promedio
//se buscan y si existen sus coordenadas se las ingresa.

fa.open("frontera_brasil.txt");
if(fa.is_open()==true)
{

    double la1,la2,lg1,lg2,d;
    int n;
    fa>>n;
    for(int i=1;i<=n;i++)
    {
        fa>>la1>>lg1>>la2>>lg2>>d;
        if(d>mu)
        {
            int id1=-1,id2=-1;
            for(int ci=0;ci<Ciudades.size();ci++)
            {
                if((Ciudades[ci].latitud()==la1)and(Ciudades[ci].longitud()==lg1))
                    id1=ci;
                if((Ciudades[ci].latitud()==la2)and(Ciudades[ci].longitud()==lg2))
                    id2=ci;
            }
            if((id1!=-1)and(id2!=-1))
            {
                Ciudades[id1].frontera.push_back(id2);
                Ciudades[id2].frontera.push_back(id1);
            }
            else
                cout<<"**** id no encontrados..."<<endl;
        }
    }
}
else
{
    cout<<"Problemas de lectura frontera..."<<endl;
}
fa.close();

//EVERY OUTLET HAS TO HAVE AT LEAST ONE BORDER CITY. 
//IF SOME OUTLET HAS NO BORDER, THE CLOSEST OUTLET IS INCLUDED
//Toda ciudad debe tener frontera. Si la zona ci no tiene frontera, entonces se ingresa 
//la zona más cercana a ci
for(int ci=0;ci<Ciudades.size();ci++)
{
    Ciudades[ci].frontera.sort();
    Ciudades[ci].frontera.unique();
    if(Ciudades[ci].frontera.size()==0)
    {
        int cmin=-1;
        double dmin=+INF , dt;
        for(int cj=0;cj<Ciudades.size();cj++)
        {
            if(ci!=cj)
            {
                dt=distancia_coord(Ciudades[ci],Ciudades[cj]);
                if((dt>0)and(dt<dmin))
                {
                    cmin=cj;
                    dmin=dt;
                }                
            }
        }
        Ciudades[cmin].frontera.push_back(ci);
        Ciudades[ci].frontera.push_back(cmin);
    }
}

//LOADING DATA OF THE FRANCHISER: SPOLETO
cout<<"Cargando datos franquicia..."<<endl;

 double F,r,G;
 int m,t;
 vector<double> Rk, Sk;
 vector<int> Bi;

ifstream ff("Instancia_SPOLETO.txt");
if(ff.is_open()==true)
{
    ff>>m>>t>>F>>r>>G;
    Rk.resize(t);
    Sk.resize(t);

    for(int i=0;i<Rk.size();i++)
        ff>>Rk[i];

    for(int i=0;i<Sk.size();i++)
        ff>>Sk[i];
    
    double lt, lg; 
    int bi;
    while(ff.eof()!=true)
    {
        ff>>lt>>lg>>bi;
        bool encontrado=false;
        for(int ci=0;ci<Ciudades.size();ci++)
        {
            if((Ciudades[ci].latitud()==lt)and(Ciudades[ci].longitud()==lg))
            {
                encontrado=true;
                Ciudades[ci].Bi+=bi;//Se carga Bi=Número de locales del mismo tipo de producto
                break;
            }
        }
        if(encontrado==false)
            cout<<"Coordenadas NO encontradas "<<bi<<endl;
    }
}
else
{
    cout<<"Error de apertura datos franquicia"<<endl;
}
ff.close();

cout<<"PARÁMETROS:"<<endl
    <<"F="<<F<<endl
    <<"r="<<r<<endl
    <<"alfa="<<alfa<<endl
    <<"G="<<G<<endl
    <<"L="<<L<<endl
    <<"R(k)  ,  S(k)"<<endl;
for(int i=0;i<Sk.size();i++)
    cout<<Rk[i]<<" , "<<Sk[i]<<endl;

//CONSTRUCTING THE IP MODEL
int n_locales=0;
ofstream f_sol("Resultados.csv",ios::app);

// for(L=3;L<=8;L++)
//for(alfa=0.002;alfa<=0.041;alfa=alfa+0.002)
    {
    Modelo MF;
    MF.incluir_parametros();
    MF.crear_vars(Ciudades,F,r,Rk,Sk);
    MF.crear_restricciones(Ciudades, G, F, r,  Rk,  Sk);
    //MF.write_MIP("MFranquicias.lp");
    MF.resolver();
    /*
    ofstream f("sol_grafico.txt",ios::app);
    MF.presentar_sol(f,n_locales,Ciudades);
    f.close();
    */
    MF.presentar_sol(cout,n_locales,Ciudades);
//SAVING THE SOLUTION
    f_sol<<"("<<L<<" / "<<mu<<" / "<<alfa<<");"<<MF.Obt_Obj()<<";"<<n_locales<<";"
         <<MF.Obt_tiempo()<<";"<<MF.ObtN_Vars()<<";"<<MF.ObtN_Restr()<<endl;
    }	
    f_sol.close();
    return 0;
}

//https://www.ehowenespanol.com/calcular-distancia-puntos-latitud-longitud-como_452715/
//Fórmula de Haversine
double distancia_coord(ciudad c1,ciudad c2)
{
    pair<double,double> delta;
    delta=make_pair(c2.latitud()-c1.latitud(),c2.longitud()-c1.longitud());
    //a = [sin ² (Δlat / 2) + cos (lat1)] x cos (lat2) x sin ² (Δlong / 2)
    double a,c,d,R=6371;
    a=pow(sin(delta.first/2.0),2.0)+cos(c2.latitud())*cos(c1.latitud())*pow(sin(delta.second/2.0),2.0);
    //c = 2 x arctan (√ a / √ (1-a))
    c=2*atan(sqrt(a)/sqrt(1-a));
    //d = R xc, en donde "R" representa el radio de la Tierra, o 6371
    d=R*c;
    return d;
}

Modelo::Modelo()
{
    env_modelo = new GRBEnv();
    modelo_ip = new GRBModel(*env_modelo);
    modelo_ip->set(GRB_StringAttr_ModelName, "Modelo de franquicias");
    modelo_ip->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
}

Modelo::~Modelo()
{
    for(int i=0;i<_x.size();i++)
        _x[i].clear();

    _x.clear();
    _y.clear();
    
    delete modelo_ip;
    delete env_modelo;
}

void Modelo::incluir_parametros()
{
    modelo_ip->getEnv().resetParams();
    
    modelo_ip->getEnv().set(GRB_IntParam_MIPFocus, 1);
    modelo_ip->getEnv().set(GRB_IntParam_Symmetry, 2);
    modelo_ip->getEnv().set(GRB_DoubleParam_NodefileStart, 0.5);

    //modelo_ip->getEnv().set(GRB_IntParam_NoRelHeuristic,1);
    //modelo_ip->getEnv().set(GRB_DoubleParam_Heuristics,0.1);
    
    modelo_ip->getEnv().set(GRB_IntParam_Seed, 1);
    modelo_ip->getEnv().set(GRB_DoubleParam_TimeLimit, TIEMPO_IP);
    modelo_ip->getEnv().set(GRB_DoubleParam_MIPGap, GAP);
    
    modelo_ip->update();
}


void Modelo::crear_vars(vector<ciudad> C, double F,double r, vector<double> Rk, vector<double> Sk)
{

    _Ady.resize(C.size());
    for(int i=0;i<_Ady.size();i++)
        _Ady[i].resize(C.size(),false);

    for(int i=0;i<C.size();i++)
    {
        list<int>::iterator itf;
        for(itf=C[i].frontera.begin();itf!=C[i].frontera.end();itf++)
        {
            _Ady[i][*itf]=true;
            _Ady[*itf][i]=true;
        }
    }


    int m, t;
    m=C.size();
    t=C[0].demanda.size();

    _y.resize(m);
    _x.resize(m);
    for(int i=0;i<_x.size();i++)
        _x[i].resize(m);
        
    for(int i=0;i<_y.size();i++)
    {
        double costo=0;
        for(int k=0;k<Rk.size();k++)
            costo=costo+(pow(1+r,-(k+1))*(-Sk[k]+(Rk[k]*C[i].demanda[k]/(C[i].Bi+1))));
        
        costo+=F;

        ostringstream nombre_y;
        nombre_y << "y_"<<i;
        _y[i] = modelo_ip->addVar(0.0, 1.0, costo , GRB_BINARY, nombre_y.str());
        
    }


    for(int i=0;i<_x.size();i++)
    {
        list<int>::iterator itj;
        for(itj=C[i].frontera.begin();itj!=C[i].frontera.end();itj++)
        {
                double costo=0,dij;
                dij=distancia_coord(C[i],C[*itj]);
                for(int k=0;k<Rk.size();k++)
                    costo=costo+pow(1+r,-(k+1))*((Rk[k]*C[*itj].demanda[k])/(C[*itj].Bi+1))*exp(-alfa*dij);
                
                if(costo<0.1)
                    costo=0;
                
                ostringstream nombre_x;
                nombre_x << "x_"<<i<<"_"<<*itj;
                _x[i][*itj] = modelo_ip->addVar(0.0, 1.0, costo , GRB_BINARY, nombre_x.str());
                
        }
    }

    modelo_ip->update();
}

void Modelo::crear_restricciones(vector<ciudad> C, double G,double F,double r,  vector<double> Rk, vector<double> Sk)
{
 /*  
cout<<"Restricción extra..."<<endl;
for(int i=0;i<m;i++)
{
    ostringstream restr0;
    restr0 << "Rest_0_" <<i;
    GRBLinExpr eq_restr0 = 0;

    double costo_x=0;

     for(int j=0;j<m;j++)
        if(i!=j)
        {
            for(int k=0;k<Rk.size();k++)
            {
                costo_x = (Dik[j][k]/(Bi[j]+1)*exp(-alfa_exp*dij[i][j]));
            }
            eq_restr0+=costo_x*_x[i][j];
        }
    modelo_ip->addConstr(eq_restr0<=ceil(0.1*Dik[i][0]/(Bi[i]+1)) , restr0.str());
}
*/

int m, t;
m=C.size();
t=C[0].demanda.size();
    cout<<"Restricción 1..."<<endl;

    for(int j=0;j<_x.size();j++)
    {
        ostringstream restr1;
        restr1 << "Rest_1_" <<j;
        GRBLinExpr eq_restr1 = 0;
        for(int i=0;i<_x[j].size();i++)
        {
            if((j!=i)and(_Ady[i][j]==1))
            {
                eq_restr1+=_x[i][j]; 
            }
        }
        modelo_ip->addConstr(eq_restr1+_y[j]<=1 , restr1.str());      
    }

    cout<<"Restricción 2..."<<endl;

    for(int i=0;i<_x.size();i++)
    {
        ostringstream restr2;
        restr2 << "Rest_2_" <<i;
        GRBLinExpr eq_restr2 = 0;
        for(int j=0;j<_x[i].size();j++)
        {
            if((j!=i)and(_Ady[i][j]>0))
            {
                eq_restr2+=_x[i][j]; 
                //modelo_ip->addConstr(_x[i][j]<=_y[i]); 
            }
        }
        modelo_ip->addConstr(eq_restr2<=_y.size()*_y[i] , restr2.str());      
    }



cout<<"Restricción 3..."<<endl;
ostringstream restr3;
restr3 << "Rest_3_"<<L;
GRBLinExpr eq_restr3 = 0;
for(int i=0;i<_y.size();i++)
    eq_restr3+=_y[i];

modelo_ip->addConstr(eq_restr3<= L , restr3.str()); 
 

cout<<"Restricción 4..."<<endl;

for(int i=0;i<_y.size();i++)
{
    list<int>::iterator itj;
    for(itj=C[i].frontera.begin();itj!=C[i].frontera.end();itj++)
        if(distancia_coord(C[i],C[*itj])<2*mu)
        {
            ostringstream restr4;
            restr4 << "Rest_4_" <<i<<"_"<<*itj;
            modelo_ip->addConstr(_y[i]+_y[*itj]<= 1, restr4.str());
        }
}

cout<<"Restricción 5..."<<endl;

for(int i=0;i<m;i++)
{
    ostringstream restr5;
    restr5 << "Rest_5_" <<i;
    GRBLinExpr eq_restr5 = 0;

    double costo_y=0, costo_x=0;
    for(int k=0;k<Rk.size();k++)
    {
        costo_y=costo_y+(pow(1+r,-(k+1))*((Rk[k]*C[i].demanda[k]/(C[i].Bi+1))));
        //costo_y=costo_y+(pow(1+r,k+1)*((1-alfa)*(Rk[k]*C[i].demanda[k]/(C[i].Bi+1))));
    }
    costo_y-=F;
    
    for(int k=0;k<t;k++)
    {
        costo_y=costo_y-(pow(1+r,-(k+1))*G);
    }
    
    eq_restr5+=costo_y*_y[i];


     for(int j=0;j<m;j++)
        if((i!=j)and(_Ady[i][j]==1))
        {
            double dij=distancia_coord(C[i],C[j]);
            for(int k=0;k<Rk.size();k++)
            {
                //double aux= pow(1+r,k+1)*(1-alfa)*Rk[k];
                double aux= pow(1+r,-(k+1))*Rk[k];
                costo_x = aux*(C[i].demanda[k]/(C[i].Bi+1)*exp(-alfa*dij));
            }
            eq_restr5+=costo_x*_x[i][j];
        }
    modelo_ip->addConstr(eq_restr5>=0  , restr5.str());
}

    modelo_ip->update();
}



ostream& Modelo::presentar_sol(ostream& os, int &nl, vector<ciudad> C)
{
    int st = modelo_ip->get(GRB_IntAttr_Status);
    
    if (modelo_ip->get(GRB_DoubleAttr_MIPGap)<1e+20)
    {
        
        os << "status=" << st << endl
        << "F objetivo=" << modelo_ip->get(GRB_DoubleAttr_ObjVal) << endl
        << "Variables: " << modelo_ip->get(GRB_IntAttr_NumVars) << endl
        << "Restricciones: " << modelo_ip->get(GRB_IntAttr_NumConstrs) << endl
        << "Tiempo=" << modelo_ip->get(GRB_DoubleAttr_Runtime) << endl
        << "Gap=" << 100 * modelo_ip->get(GRB_DoubleAttr_MIPGap)
        << endl<<endl;
        
        nl=0;
        double costo_i=0,costo_ij=0; 
        for(int i=0;i<_y.size();i++)
        {
            if(_y[i].get(GRB_DoubleAttr_X)>0.1)
            {
                nl++;
                costo_i=costo_i+_y[i].get(GRB_DoubleAttr_Obj);
                os<<endl<<ROJO<<"Punto de venta  "<<endl<<C[i]<<":"<<endl;
                for(int j=0;j<_x[i].size();j++)
                {
                    if((i!=j)and(_Ady[i][j]==1))
                    {
                        double val = 0.0;
                        val = _x[i][j].get(GRB_DoubleAttr_X);
                        if (val > 0)
                        {
                            costo_ij=costo_ij+_x[i][j].get(GRB_DoubleAttr_Obj);
                            os  <<AZUL<< _x[i][j].get(GRB_StringAttr_VarName) << "="
                                << val <<"-->"<<C[j]<<endl;
                                
                        }
                    }
                }
            }
        }
        os<<"Costo_i  = "<<costo_i<<endl
          <<"Costo_ij = "<<costo_ij<<endl
          <<"CT="<<costo_i+costo_ij<<endl;
    }
    
    return os;
}
/*
void leer_instancia(string ARCHIVO,int&m, int&t, double &F, double& r,
double& alfa, double& G,  
vector<double> &Rk,vector<double>& Sk,vector<int>& Bi,
vector<vector<int> > &Dik, vector<vector<double> > &dij,vector<vector<int> > &fij)
{

    ifstream fo(ARCHIVO.c_str());
    if(fo.is_open()==true)
    {
//LECTURA DE PARÁMETROS
        fo>>m>>t>>F>>r>>alfa>>G;
        Rk.resize(t);
        Sk.resize(t);
        Bi.resize(m);

        Dik.resize(m);
        dij.resize(m);
        fij.resize(m);
        for(int i=0;i<Dik.size();i++)
        {
            Dik[i].resize(t);
            dij[i].resize(m);
            fij[i].resize(m);
        }
        for(int i=0;i<Rk.size();i++)
            fo>>Rk[i];

        for(int i=0;i<Sk.size();i++)
            fo>>Sk[i];
        
        for(int i=0;i<Bi.size();i++)
            fo>>Bi[i];
        
        for(int i=0;i<Dik.size();i++)
        {
            for(int j=0;j<Dik[i].size();j++)
            {
                 fo>>Dik[i][j];
            }
        }

        for(int i=0;i<dij.size();i++)
        {
            for(int j=0;j<dij[i].size();j++)
            {
                 fo>>dij[i][j];
            }
        } 

        for(int i=0;i<fij.size();i++)
        {
            for(int j=0;j<fij[i].size();j++)
            {
                 fo>>fij[i][j];
            }
        } 
    }
    fo.close();
}
*/