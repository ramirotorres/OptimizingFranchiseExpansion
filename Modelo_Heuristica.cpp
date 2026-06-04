#include "constantes.h"


// PARAMETERS OF THE MODELS
//AVERAGE DISTANCE
double d_prom=24.19;
double mu=4*d_prom;
//EFECT OF THE DISTANCE PN THE DEMAND
double alfa=0.01;
//MAXIMUM NUMBER OF OUTLETS
int L=3;

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

ostream& operator<<(ostream& os, vector<int> x)
{
    for(int j=0;j<x.size();j++)
        os<<x[j]<<" ";
        
    os<<endl;
    return os;
}

class ciudad
{
private:
    int id;
    string Nombre;
    pair<double,double> coordenada;
public:
    vector<double> demanda;
    list<int> frontera;
    long int Bi;

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
  //  Solver
GRBEnv* env_modelo;
GRBModel*  modelo_ip;
    // Vector of variables

  vector< vector<GRBVar> >  _x;//x_ij
  vector< GRBVar >  _y;//y_i

public:
  // Constructor and destructor
  Modelo();
  ~Modelo();

void incluir_parametros();
void crear_vars(vector<ciudad> C, double F,double r, vector<double> Rk, vector<double> Sk);
void crear_restricciones(vector<ciudad> C, double G,double F,double r, vector<double> Rk, vector<double> Sk);

void fijar_vars_Y(vector<int> S);
void Cargar_Sol_Inicial(vector<int> S);

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

double Obj_var_Y(int i)
{
    return _y[i].get(GRB_DoubleAttr_Obj);
}
};

class Modelo_Reducido
{
private:
  //  Solver
GRBEnv* env_modelo;
GRBModel*  modelo_ip;
    // Vector of variables

  vector< vector<GRBVar> >  _x;//x_ij
  vector< GRBVar >  _y;//y_i

public:

  // Constructor y destructor
  Modelo_Reducido();
  ~Modelo_Reducido();

void incluir_parametros();
void crear_vars(vector<ciudad> C, double F,double r, vector<double> Rk, vector<double> Sk);
void crear_restricciones(vector<ciudad> C, double G,double F,double r, vector<double> Rk, vector<double> Sk);

void fijar_vars_Y(vector<int> S);
void limpiar_vars_Y();

ostream& presentar_sol(ostream& os, vector<ciudad> C);

void write_MIP(string archivo)
{
    modelo_ip->write(archivo);
}
void resolver()
{
    modelo_ip->optimize();
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

int status()
{
    return modelo_ip->get(GRB_IntAttr_Status);
}

};

double distancia_coord(ciudad c1,ciudad c2);

double Utilidad_yi(ciudad Ci, double F,double r, vector<double> Rk, vector<double> Sk);
double Coef_R5_yi(ciudad Ci, double G,double F,double r,  vector<double> Rk);
double Utilidad_xij(ciudad Ci,ciudad Cj, double r, vector<double> Rk);

int variante_1(vector<ciudad> Ciudades,vector<int> S,vector<int> C_factible);
int variante_2(vector<ciudad> Ciudades,vector<int> S,vector<int> C_factible);
void Heuristica(vector<ciudad> Ciudades,double F,double r,vector<double> Rk,vector<double> Sk, double G, int L, int MAX_ITR,double &zf, vector<int> &Sf);


int main()
{
    cout<<endl<<endl
        <<"***********************"<<endl
        <<" FRANCHISE HERISTIC: "<<endl
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
    cout<<"Error:  Cities file..."<<endl;
}
fa.close();

//Se ingresan como áreas de frontera si son menores que el promedio
//AREAS ARE CONSIDERED IN THE BORDER IF THE DISTANCE IS LESS THAN THE AVERAGE
cout<<"... Computing borders "<<endl;
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
    cout<<"ERROR: border file..."<<endl;
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

//LOADING DATA OF THE FRANCHISER: 
cout<<"Loading franchise data ..."<<endl;

 double F,r,G;
 int m,t;
 vector<double> Rk, Sk;
 vector<int> Bi;

ifstream ff("Instance_Brasil.txt");
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
    cout<<"Error: Franchise file..."<<endl;
}
ff.close();

cout<<"PARAMETERS:"<<endl
    <<"F="<<F<<endl
    <<"r="<<r<<endl
    <<"alfa="<<alfa<<endl
    <<"G="<<G<<endl
    <<"L="<<L<<endl
    <<"R(k)  ,  S(k)"<<endl;
for(int i=0;i<Sk.size();i++)
    cout<<Rk[i]<<" , "<<Sk[i]<<endl;


//CONSTRUCTING THE IP MODEL y HEURISTICA
int n_locales=0;
ofstream f_sol("Resultados.txt",ios::app);

for(L=3;L<=8;L++)
//for(alfa=0.002;alfa<=0.041;alfa=alfa+0.002)
    {
    //HEURISTIC

    vector<int> Sf;
    double zf=0;
    int MAX_ITR=50;
    double th1,th2,th;
    th1=clock();

    Heuristica(Ciudades,F,r,Rk,Sk,G,L,MAX_ITR,zf,Sf);
    th2=clock();
    th=(th2-th1)/CLOCKS_PER_SEC;


    // LP MODEL
    
    Modelo MF;
    MF.incluir_parametros();
    MF.crear_vars(Ciudades,F,r,Rk,Sk);
    MF.crear_restricciones(Ciudades, G, F, r,  Rk,  Sk);

    MF.Cargar_Sol_Inicial(Sf);
    MF.resolver();
    /*
    ofstream f("sol_grafico.txt",ios::app);
    MF.presentar_sol(f,n_locales,Ciudades);
    f.close();
    */
    MF.presentar_sol(cout,n_locales,Ciudades);
//SAVING THE SOLUTION
    f_sol<<"("<<L<<" / "<<mu<<" / "<<alfa<<");"<<MF.Obt_Obj()<<";"<<n_locales<<";"
         <<MF.Obt_tiempo()<<";"<<MF.ObtN_Vars()<<";"<<MF.ObtN_Restr()<<";;"
         <<zf<<";"<<Sf.size()<<";"<<th<<";"<<MAX_ITR<<endl;
    }
    f_sol.close();
    return 0;
}

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

double Utilidad_yi(ciudad Ci, double F,double r, vector<double> Rk, vector<double> Sk)
{
    double costo=0;
    for(int k=0;k<Rk.size();k++)
        costo=costo+(pow(1+r,-(k+1))*(-Sk[k]+(Rk[k]*Ci.demanda[k]/(Ci.Bi+1))));

    costo+=F;
    return costo;
}

double Coef_R5_yi(ciudad Ci, double G,double F,double r,  vector<double> Rk)
{
    double costo_y=0;
    for(int k=0;k<Rk.size();k++)
    {
        costo_y=costo_y+(pow(1+r,-(k+1))*((Rk[k]*Ci.demanda[k]/(Ci.Bi+1))));
    }
    costo_y-=F;

    for(int k=0;k<Rk.size();k++)
    {
        costo_y=costo_y-(pow(1+r,-(k+1))*G);
    }

    return costo_y;
}

double Utilidad_xij(ciudad Ci,ciudad Cj, double r, vector<double> Rk)
{
    double costo=0,dij;
    dij=distancia_coord(Ci,Cj);
    for(int k=0;k<Rk.size();k++)
        costo=costo+pow(1+r,-(k+1))*((Rk[k]*Cj.demanda[k])/(Cj.Bi+1))*exp(-alfa*dij);

    return costo;
}

int variante_1(vector<ciudad> Ciudades , vector<int> S , vector<int> C_factible)
{
    int j;
    double d_total=0;
    vector<double> C_demanda(C_factible.size(),0);
    for(int c=0;c<C_factible.size();c++)
    {
        C_demanda[c]=0.0;
        for(int cd=0;cd<Ciudades[ C_factible[c] ].demanda.size();cd++)
             C_demanda[c]+=Ciudades[ C_factible[c] ].demanda[cd];
        d_total+=C_demanda[c];
    }

    bool d_bool;
    int cont=0;
    do
    {
        cont++;
        double prob=rand()%10001/10000.0;
        double prob_acum=0.0;
        j=C_factible.back();
        for(int c=0;c<C_demanda.size();c++)
        {
                prob_acum+=(C_demanda[c]/d_total);
                if(prob_acum>prob)
                {
                    j=C_factible[c];
                    break;
                }
        }

        d_bool=true;
        if(S.size()>0)
        {
            for(int cs=0;cs<S.size();cs++)
            {
                if(distancia_coord(Ciudades[S[cs]],Ciudades[j])<=mu)
                {
                    d_bool=false;
                    break;
                }
            }
        }
    }
    while((d_bool==false)and(cont<C_factible.size()));

    return j;
}

int variante_2(vector<ciudad> Ciudades,vector<int> S,vector<int> C_factible)
{
    int j;

    double d_total=0;
    vector<double> C_demanda(C_factible.size(),0);
    for(int c=0;c<C_factible.size();c++)
    {
         C_demanda[c]=0.0;
        for(int cd=0;cd<Ciudades[ C_factible[c] ].demanda.size();cd++)
             C_demanda[c]+=Ciudades[ C_factible[c] ].demanda[cd];
        d_total+=C_demanda[c]; 
    }

    if(S.size()==0)
        j=variante_1(Ciudades,S, C_factible);
    else
    {
        for(int i=0;i<S.size();i++)
            for(int j=0;j<C_factible.size();j++)
                if(S[i]==C_factible[j])
                    C_demanda[ j ]=-1;
            
            
        bool d_bool;
        int cont=0;
        do
        {
            cont++;
            double d_max=0;
            int ind=-1;
            for(int c=0;c<C_demanda.size();c++)
            {
                    if(C_demanda[c]>d_max)
                    {
                        j=C_factible[c];
                        d_max=C_demanda[c];
                        ind=c;
                    }
            }

            C_demanda[ind]=-1;
            d_bool=true;

            for(int cs=0;cs<S.size();cs++)
            {
                if(distancia_coord(Ciudades[S[cs]],Ciudades[j])<=mu)
                    d_bool=false;

            }
        }
        while((d_bool==false)and(cont<C_factible.size()));
    }
    return j;
}


void Heuristica(vector<ciudad> Ciudades,double F,double r,vector<double> Rk,vector<double> Sk, double G, int L, int MAX_ITR,double &zf, vector<int> &Sf)
{

//Validando cuales ciudades satisfacen la condición de ganancia minima.

    vector<int> Y_factible;

    double costo_yi,costo_xij;
    for(int i=0;i<Ciudades.size();i++)
    {
        costo_yi= Coef_R5_yi(Ciudades[i], G, F, r, Rk);

        costo_xij=0;
        list<int>::const_iterator itf;
        for(itf=Ciudades[i].frontera.begin();itf!=Ciudades[i].frontera.end();itf++)
        {
            costo_xij = costo_xij+Utilidad_xij(Ciudades[i],Ciudades[*itf], r,  Rk);
        }

        if(costo_xij+costo_yi>0)
        {
            Y_factible.push_back(i);
        }
    }
cout<<"Y_factible="<<Y_factible.size()<<endl;
//Aplicado heurística
    vector<int> Y_fijas;//indice de variables Y fijadas
    int j;
    double z;//obj de la iteración
    zf=0;//Solución inicial
    Sf.clear();//Solución inicial

    if(Y_factible.size()>L)
    {
        Modelo_Reducido MH;
        MH.incluir_parametros();
        MH.crear_vars(Ciudades,F,r,Rk,Sk);
        MH.crear_restricciones(Ciudades, G, F, r,  Rk,  Sk);
        //MH.write_MIP("Red.lp");

        for (int p=1; p<=MAX_ITR;p++)
        {
            cout<<"---->Iteración="<<p<<endl;
            Y_fijas.clear();
            for(int sd=1;sd<=L;sd++)
            {
                //j=variante_1(Ciudades,Y_fijas,Y_factible);
                j=variante_2(Ciudades,Y_fijas,Y_factible);
                Y_fijas.push_back(j);
            
                cout<<"...Fijando S="<<Y_fijas<<endl;
                MH.fijar_vars_Y(Y_fijas);
                MH.resolver();
                if(MH.status()==2)
                {
                    if(MH.Obt_Obj()>zf)
                    {
                        zf=MH.Obt_Obj();
                        Sf=Y_fijas;
                        cout<<".... Obj actualizado ="<<zf<<"-->"<<Sf<<endl;
                    }
                }
                else
                    Y_fijas.pop_back();

                MH.limpiar_vars_Y();
                
            }
        }
    }
    else
    {
        Sf=Y_factible;
    }

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
    modelo_ip->getEnv().set(GRB_DoubleParam_NodefileStart, 0.5);


    modelo_ip->getEnv().set(GRB_IntParam_Seed, 1);
    modelo_ip->getEnv().set(GRB_DoubleParam_TimeLimit, TIEMPO_IP);
    modelo_ip->getEnv().set(GRB_DoubleParam_MIPGap, GAP);

    modelo_ip->update();
}



void Modelo::crear_vars(vector<ciudad> C, double F,double r, vector<double> Rk, vector<double> Sk)
{

    _y.resize( C.size() );
    _x.resize( C.size() );
    for(int i=0;i<_x.size();i++)
        _x[i].resize( C.size() );

cout<<"... Vars Y"<<endl;
    for(int i=0;i<_y.size();i++)
    {
        double costo_yi=Utilidad_yi(C[i],F,r,Rk,Sk);
       
        ostringstream nombre_y;
        nombre_y << "y_"<<i;
        _y[i] = modelo_ip->addVar(0.0, 1.0, costo_yi , GRB_BINARY, nombre_y.str());
    }

cout<<"... Vars X"<<endl;
    for(int i=0;i<_x.size();i++)
    {
        list<int>::iterator itj;
        for(itj=C[i].frontera.begin();itj!=C[i].frontera.end();itj++)
        {
            double costo_xij=Utilidad_xij(C[i],C[*itj],r,Rk);
          
            ostringstream nombre_x;
            nombre_x << "x_"<<i<<"_"<<*itj;
            _x[i][*itj] = modelo_ip->addVar(0.0, 1.0, costo_xij , GRB_CONTINUOUS, nombre_x.str());

        }
    }

    modelo_ip->update();
}

void Modelo::crear_restricciones(vector<ciudad> C, double G,double F,double r,  vector<double> Rk, vector<double> Sk)
{

cout<<"Restricción 1..."<<endl;

    for(int j=0;j<_x.size();j++)
    {
        ostringstream restr1;
        restr1 << "Rest_1_" <<j;
        GRBLinExpr eq_restr1 = 0;
        list<int>::iterator iti;
        for(iti=C[j].frontera.begin();iti!=C[j].frontera.end();iti++)
        {
            eq_restr1+=_x[*iti][j];
        }
        modelo_ip->addConstr(eq_restr1+_y[j]<=1 , restr1.str());
    }

cout<<"Restricción 2..."<<endl;

    for(int i=0;i<_x.size();i++)
    {
        ostringstream restr2;
        restr2 << "Rest_2_" <<i;
        GRBLinExpr eq_restr2 = 0;
        list<int>::iterator itj;
        for(itj=C[i].frontera.begin();itj!=C[i].frontera.end();itj++)
        {
            eq_restr2+=_x[i][*itj];
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

    for(int i=0;i<C.size();i++)
    {
        ostringstream restr5;
        restr5 << "Rest_5_" <<i;
        GRBLinExpr eq_restr5 = 0;

        double costo_x=0;

        eq_restr5+=Coef_R5_yi(C[i],G, F, r, Rk)*_y[i];

        list<int>::const_iterator itf;
        for(itf=C[i].frontera.begin();itf!=C[i].frontera.end();itf++)
            {
                eq_restr5+=Utilidad_xij(C[i],C[*itf],r,Rk)*_x[i][*itf];
            }
        modelo_ip->addConstr(eq_restr5>=0  , restr5.str());
    }

    modelo_ip->update();
}


void Modelo::fijar_vars_Y(vector<int> S)
{
    
    for(int i=0;i<_y.size();i++)
    {
        _y[i].set(GRB_DoubleAttr_UB,0.0);
    }
    
    for(int i=0;i<S.size();i++)
    {
        _y[S[i]].set(GRB_DoubleAttr_LB,1.0);
        _y[S[i]].set(GRB_DoubleAttr_UB,1.0);
    }
    
    modelo_ip->update();
}

void Modelo::Cargar_Sol_Inicial(vector<int> S)
{
    cout<<"... Cargando Solución Incial Y="<<S<<endl;
    
    for(int i=0;i<S.size();i++)
    {
        _y[S[i]].set(GRB_DoubleAttr_Start,1.0);
    }
    
   // modelo_ip->update();
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
                
                list<int>::iterator itj;
                for(itj=C[i].frontera.begin();itj!=C[i].frontera.end();itj++)
                {
                    double val = 0.0;
                    val = _x[i][*itj].get(GRB_DoubleAttr_X);
                    if (val > 0)
                    {
                        costo_ij=costo_ij+_x[i][*itj].get(GRB_DoubleAttr_Obj);
                        os  <<AZUL<< _x[i][*itj].get(GRB_StringAttr_VarName) << "="
                            << val <<"-->"<<C[*itj]<<endl;
                    }
                }
            }
        }
        os<<"Costo_i  = "<<costo_i<<endl
          <<"Costo_ij = "<<costo_ij<<endl
          <<"CT="<<costo_i+costo_ij<<endl
          <<"CT_IP="<<modelo_ip->get(GRB_DoubleAttr_ObjVal)<<endl;
    }

    return os;
}

//***************************************************************************
//***************** MODELO REDUCIDO *****************************************
//***************************************************************************

Modelo_Reducido::Modelo_Reducido()
{
    env_modelo = new GRBEnv();
    modelo_ip = new GRBModel(*env_modelo);
    modelo_ip->set(GRB_StringAttr_ModelName, "Modelo de franquicias Reducido");
    modelo_ip->set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
}

Modelo_Reducido::~Modelo_Reducido()
{
    for(int i=0;i<_x.size();i++)
        _x[i].clear();

    _x.clear();
    _y.clear();

    delete modelo_ip;
    delete env_modelo;
}

void Modelo_Reducido::incluir_parametros()
{
    modelo_ip->getEnv().resetParams();
    modelo_ip->getEnv().set(GRB_DoubleParam_NodefileStart, 0.5);
    modelo_ip->getEnv().set(GRB_IntParam_Seed, 1);
    modelo_ip->getEnv().set(GRB_DoubleParam_TimeLimit, TIEMPO_IP);
    modelo_ip->getEnv().set(GRB_DoubleParam_MIPGap, GAP);

    modelo_ip->update();
}


void Modelo_Reducido::crear_vars(vector<ciudad> C, double F,double r, vector<double> Rk, vector<double> Sk)
{

    int m, t;
    m=C.size();
    t=C[0].demanda.size();

    _y.resize(m);
    _x.resize(m);
    for(int i=0;i<_x.size();i++)
        _x[i].resize(m);

    for(int i=0;i<_y.size();i++)
    {
        double costo_yi=Utilidad_yi(C[i],F, r, Rk, Sk);

        ostringstream nombre_y;
        nombre_y << "y_"<<i;
        _y[i] = modelo_ip->addVar(0.0, 1.0, costo_yi , GRB_CONTINUOUS , nombre_y.str());

    }


    for(int i=0;i<_x.size();i++)
    {
        list<int>::iterator itj;
        for(itj=C[i].frontera.begin();itj!=C[i].frontera.end();itj++)
        {
                double costo_xij=Utilidad_xij(C[i],C[*itj],r, Rk);
                
                ostringstream nombre_x;
                nombre_x << "x_"<<i<<"_"<<*itj;
                _x[i][*itj] = modelo_ip->addVar(0.0, 1.0, costo_xij , GRB_CONTINUOUS, nombre_x.str());
        }
    }

    modelo_ip->update();
}


void Modelo_Reducido::crear_restricciones(vector<ciudad> C, double G,double F,double r,  vector<double> Rk, vector<double> Sk)
{


cout<<"Restricción 1..."<<endl;

    for(int j=0;j<_x.size();j++)
    {
        ostringstream restr1;
        restr1 << "Rest_1_" <<j;
        GRBLinExpr eq_restr1 = 0;
        list<int>::iterator iti;
        for(iti=C[j].frontera.begin();iti!=C[j].frontera.end();iti++)
        {
            eq_restr1+=_x[*iti][j];
        }
        modelo_ip->addConstr(eq_restr1+_y[j]<=1 , restr1.str());
    }

cout<<"Restricción 2..."<<endl;

    for(int i=0;i<_x.size();i++)
    {
        ostringstream restr2;
        restr2 << "Rest_2_" <<i;
        GRBLinExpr eq_restr2 = 0;
        list<int>::iterator itj;
        for(itj=C[i].frontera.begin();itj!=C[i].frontera.end();itj++)
        {
            eq_restr2+=_x[i][*itj];
        }
        modelo_ip->addConstr(eq_restr2<=_y.size()*_y[i] , restr2.str());
    }

 


cout<<"Restricción 5..."<<endl;

    for(int i=0;i<C.size();i++)
    {
        ostringstream restr5;
        restr5 << "Rest_5_" <<i;
        GRBLinExpr eq_restr5 = 0;

        eq_restr5+=Coef_R5_yi(C[i], G, F, r, Rk)*_y[i];

        list<int>::const_iterator itf;
        for(itf=C[i].frontera.begin();itf!=C[i].frontera.end();itf++)
            {
                double costo_xij=Utilidad_xij(C[i],C[*itf],r,Rk);
                eq_restr5+=costo_xij*_x[i][*itf];
            }
        modelo_ip->addConstr(eq_restr5>=0  , restr5.str());
    }

    modelo_ip->update();
}


void Modelo_Reducido::fijar_vars_Y(vector<int> S)
{
    
    for(int i=0;i<_y.size();i++)
    {
        _y[i].set(GRB_DoubleAttr_UB,0.0);
    }
    
    for(int i=0;i<S.size();i++)
    {
        _y[S[i]].set(GRB_DoubleAttr_LB,1.0);
        _y[S[i]].set(GRB_DoubleAttr_UB,1.0);
    }
    
    modelo_ip->update();
}


void Modelo_Reducido::limpiar_vars_Y()
{
    for(int i=0;i<_y.size();i++)
    {
        _y[i].set(GRB_DoubleAttr_LB,0.0);
        _y[i].set(GRB_DoubleAttr_UB,1.0);
    }
    modelo_ip->update();
}


ostream& Modelo_Reducido::presentar_sol(ostream& os, vector<ciudad> C)
{
    int st = modelo_ip->get(GRB_IntAttr_Status);

    //if (modelo_ip->get(GRB_DoubleAttr_MIPGap)<1e+20)
    {

        os << "status=" << st << endl
        << "F objetivo=" << modelo_ip->get(GRB_DoubleAttr_ObjVal) << endl
        << "Variables: " << modelo_ip->get(GRB_IntAttr_NumVars) << endl
        << "Restricciones: " << modelo_ip->get(GRB_IntAttr_NumConstrs) << endl
        << "Tiempo=" << modelo_ip->get(GRB_DoubleAttr_Runtime) << endl
        << endl<<endl;

        double costo_i=0,costo_ij=0;
        for(int i=0;i<_y.size();i++)
        {
            if(_y[i].get(GRB_DoubleAttr_X)>0.1)
            {
                costo_i=costo_i+_y[i].get(GRB_DoubleAttr_Obj);
                os<<endl<<ROJO<<"Punto de venta  "<<endl<<C[i]<<":"<<endl;
                list<int>::const_iterator itf;
                for(itf=C[i].frontera.begin();itf!=C[i].frontera.end();itf++)
                {
                        double val = 0.0;
                        val = _x[i][*itf].get(GRB_DoubleAttr_X);
                        if (val > 0)
                        {
                            costo_ij=costo_ij+_x[i][*itf].get(GRB_DoubleAttr_Obj);
                            os  <<AZUL<< _x[i][*itf].get(GRB_StringAttr_VarName) << "="
                                << val <<"-->"<<C[*itf]<<endl;

                        }
                }
            }
        }
        os<<"Costo_i  = "<<costo_i<<endl
          <<"Costo_ij = "<<costo_ij<<endl
          <<"CT="<<costo_i+costo_ij<<endl
          <<"CT_IP="<<modelo_ip->get(GRB_DoubleAttr_ObjVal)<<endl;
    }

    return os;
}


