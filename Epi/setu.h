/*
 *   Interface for a class implimenting sets of natural numbers
 *
 */
#ifndef    _SETU_H
#define    _SETU_H

#include <vector>

using namespace std;

//The is the number of new vertices created when a set enlarges
#define SETINCR 10

#define dna_len (int)128
#define init_bits (int)32
// TODO: Update below?
// Comment
#define max_vars (int)32
#define max_len (int)80

//fitness proportional selector used in simulations
int rselect(double *v, double ttl, int N);

class set {

 public:

  set();                 //default constructor
  set(int *m, int z);     //construct with a list of elements
  set(const set &);       //copy constructor
  ~set();                //destructor

  //utilities
  void create(int *m, int z);   //create a set with a given membership and size
  void destroy();              //destroy a set
  void setempty();             //mark as empty for mass allocation
  void copy(set &);            //copy another set
  void copyO(set &, int q);     //copy another set with an offset
  void enlarge();              //increment max
  int add(int z);              //add a member, returns true if a not ALREADY
  int remo(int z);             //remove a member, returns true if in there
  void clear();                //make the set empty

  int infected(int n, double alpha);//SIR utility routine

  //information
  int size();                //what is the size of the set
  int memb(int z);           //is z a member?  0=no 1=yes
  int memz(int z);           //zth member

  //operations
  void unyun(set &A, set &B);   //union
  void inter(set &A, set &B);   //intersection
  void setdf(set &A, set &B);   //set difference
  void symdf(set &A, set &B);   //symmetric difference
  double sumAt(double *ft);    //sum indecies of ft in the set
  int FPS(double *ft);         //fitness proportional selection
  int RNB(int v);              //Random neighbor

  //input-output
  void write(ostream &aus);     //write set
  void writememb(ostream &aus); //write members on one line
  void
  writememb(ostream &aus, char sep); //write members on one line with seperator
  void read(istream &inp);      //read set
  void readmemb(istream &inp);  //read members on one line

 private:

  int max;   //maximum number of objects in the set
  int n;     //number of objects in the set
  int *mem;  //members, sorted into increasing order


};

/***********************SET BASED GRAPH UNIT*************************/

//This is the largest number of different qualities or weights that
//may be saved in a graph
#define MAXQ 256

class graph {

 public:

  //constructors and destructors
  graph();           //initialize an empty structure
  explicit graph(int max);    //initialize to maximum of M vertices
  ~graph();          //delete s structure

  //utilities
  void create(int max);     //create with max maximum vertices
  void destroy();           //deallocate everything
  void Enlarge(int newmax); //increase maximum number of vertices to newmax
  void clearE();            //change graph to empty

  //RQ means request
  int RQquality(int v);      //index of quality or -1 if none left  init to v
  int RQquality(int *Q);     //as above, but initialize to Q  Sive is V
  void RecordQ(int num, int *Q);   //assign Q to the quality in question
  void RetrieveQ(int num, int *Q); //get the values in the quality -> Q
  void RecordQ(int num, int dex, int val);    //Q[num][dex]=val
  int RetrieveQ(int num, int dex);           //return Q[num][dex]

  //RQ means request
  int RQweight(double v);   //index of quality or -1 if none left  init to v
  int RQweight(double *W);  //as above but initialize to W.  Size is V.        
  void RecordW(int num, double *W);   //assign W to the weight in question
  void RetrieveW(int num, double *W); //get the values in the quality -> W
  void RecordW(int num, int dex, double val); //W[num][dex]=val
  double RetrieveW(int num, int dex);        //return Q[num][dex]

  int infected(int n, double alpha); //SIR utility routine
  void varInfected(bitset<dna_len> &immunity,
                   vector<pair<int, bitset<dna_len>>> &strains,
                   double alpha, int &str_id);
  static void create_new_variant(bitset<dna_len> &orig,
                                 bitset<dna_len> &variant, vector<int> rv,
                                 int low_bnd, int up_bnd);

  //initializers
  void empty(int n);            //empty graph
  void Kn(int n);               //complete 
  void Knm(int n, int m);        //complete bipartite graph
  void Cn(int n);               //cycle
  void Pn(int n, int m);         //Petersen n,m
  void Hn(int dim);             //Hypercube
  void RNGnm(int n, int m);      //Ring with +/-m neighbors
  void UTAM(int *ed);           //initialize from an upper tri. adj. matrix
  void WalkO(int *wk, int wl);   //overlaying walk representation
  void WalkT(int *wk, int wl);   //toggling walk representation
  void EdgeLst(int *el, int ne); //Edge list
  void ADTS(int **el, int L);    //implement an add delete toggle swap
  void HADTS(int **el, int L);   //implement a hop add delete toggle swap
  void copy(graph &other);      //copy another graph

  //Random graph generators

  //Barabasi-Albert graph with n vertices; m edges added at each step
  //This is a random graph model that iteratively adds a degree m
  //vertex to the growing graph, choosing its neighbors in proportion
  //to their current degree.
  void BA(int n, int m);


  //Power law clustering graph with n vertices; m edges added at each
  //step; p is the chance of adding an extra edge to create a triangle
  //(modification of BA graph to enhance clustering coeff of vertices)
  void PCG(int n, int m, double prob);

  //This is the "each edge is a coin flip" with probability p of making the edge
  void ER(int n, double p);//Erdo-Renyi graph with n vertices; prob p

  //This sort of random graph has n nodes connected to k nearest
  //neighbours with edges rewired with probability p
  void WS(int n, int k, double p);//Watts-Strogatz graph

  //This sort of random graph has n nodes connected to k nearest
  //neighbours with edges added with probability p
  void NWS(int n, int k, double p); //Newman-Watts-Strogatz graph

  //modifiers -- local operations are configured for evolutionary algorithm
  void add(int a, int b);             //force an edge to add
  void ladd(int v, int n1, int n2);    //local force edge
  void del(int a, int b);             //force an edge to be gone
  void ldel(int v, int n1, int n2);    //local remove edge
  void toggle(int a, int b);          //toggle an edge
  void loggle(int v, int n1, int n2);  //local toggle
  void simplexify(int a);            //simplexify at a
  void hop(int v, int n1, int n2);     //perform a hoperation
  void edgeswap(int a, int b, int k);  //decode and perform an edge swap
  //with degree bound k
  void Absorb(graph &other);         //add a copy of other to yourself
  void Prism();                      //create the prism of a graph

  //information
  int size();                 //number of vertices
  int edges();                //number of edges
  int edgeP(int a, int b);     //is a-b an edge?
  void dfrom(int z, int *ds);  //distances from z
  void dfrom0(int z, int *ds); //distances from z in color zero
  int ecc(int z);             //eccentricity of a vertex
  int diameter();             //eccentricity of a vertex
  int radius();               //eccentricity of a vertex
  int connectedP();           //is the graph connected?
  void eccSeq(int *ecs);      //compute the eccentricity sequence
  int nbrmod(int v, int n);    //compute the n%degreeth neighbor of v
  int degree(int v);          //report the degree of v
  double meandegree();        //report the mean degree of the graph
  int Nbrs(int v, int *nb);    //report the neighboors of v
  int MaxCol();               //report the maximal color

  void DiffChar(int v, double omega, double *dc); //diffusion character at v

  //Quality and weight handling

  //Simulation methods
  /* Run and SIR epidemic with a given patient zero; return maximum number
   * of people infected, length of epidemic, total number infected.
   * The parameter alpha is the probability of passing an infection
   */
  void SIR(int p0, int &max, int &len, int &ttl, double alpha, vector<int> &prof); //SIR method
  void SIRProfile(int p0, int &max, int &len, int &ttl, double alpha, double *prof); //Sir Method, with profile
  void varSIR(int p0, int &vcnt, vector<int> vprofs[], bitset<dna_len> variants[], int *vorigs, pair<int, int> vtimes[], double alpha, int lB, int uB, double var_prob);

  /*SIRr is the same as SIR except that patient zero is assigned at random*/
  void SIRr(int &max, int &len, int &ttl, double alpha); //SIR method

  /***********************************************************************/
  /* The following SIR method superceeds the others, they are included
     for upward compatibility.

     S=0;
     I=1;
     R=2;

  */


  //Alpha is the chance of transmission, advance is the I->R
  //probability if q is -1 then it is assigned the first time the
  //routine is called.  It is the quality that stores agent staus.
  //This routine returns the number of sick people The status of
  //patients are stored in the quality
  int SIRupdate(int &q, double alpha, double advance);


  //Alpha is the chance of transmission, creating an E, infec is the
  //chance of E advancing to I, asym is the chance of E going direct
  //to I -- that is an asymptomtic person, advance is the I->R
  //probability if q is -1 then it is assigned the first time the
  //routine is called.  It is the quality that stores agent staus.
  //This routine returns the number of sick people The status of
  //patients are stored in the quality
  int SEIRupdate(int &q, double alpha,   //pairwise infection
                 double infec,   //E->I transition
                 double asym,    //E->R transition
                 double advance  //I->R transition
  );//one time step


  //Alpha is the chance of transmission, creating an E, infec is the
  //chance of E advancing to Ep, visible is the chance of the Ep-I
  //transition, asym is the chance of Ep going to R directly, skipping
  //I, and finally advance is the chance of the I->R transition 
  //---
  //If q is -1 then it is assigned the first time the routine is
  //called.  It is the quality that stores agent staus.  This routine
  //returns the number of sick people The status of patients are
  //stored in the quality
  int SEEpIRupdate(int &q, double alpha,   //pairwise infection
                   double infec,   //E->Ep transition
                   double visible, //Ep->I transition
                   double asym,    //Ep->R transition
                   double advance  //I->R transition
  );//one time step




  //This version adds the longer term asymptomatic class A
  int SEEIARupdate(int &q, double alpha,   //pairwise infection
                   double infec,   //E->Ep transition
                   double infecP,  //Ep->I transition
                   double cksee,   //Ep->A transition
                   double gamma,   //A,I -> R transition
                   double epsilon  //R->S transition
  );//one time step





  /*
   *   Run an attack until disconnection simulation on the graph.  The
   *   vector vulnerability must be a probability distribution on the
   *   vertices.  Vertices are chosen at random until the graph
   *   becomes disconnected and the number of vertices until
   *   disconnection is returned.  Getting pr right is the users
   *   responsibility.  The routine uses the color array to mark dead
   *   vertices.
   */
  int attack(double *pr);  //probabalistic attack method

  /************ANT METHODS GO HERE*************************************/
  //Start an ant at each vertex and run it for reps steps and return
  //the number of ands entering each vertex
  void Ant(int reps, int *cnt);

  //color methods
  void setC2(int vl);         //set colors to vl
  void GDC(int *ord);         //run the greedy coloring algorithm/order ord
  int AUC(double *gn, int tg); //run Austrailian coloring with target tg

  //genetics
  int FPN(int v, double *ft); //perform fitness perprotional neighbor selection
  int RNB(int v);            //Random neighbor

  //I-O
  void write(ostream &aus);   //write the graph
  void read(istream &inp);    //read the graph
  void readadj(istream &inp, int numV); //read a python adjacency graph
  void writeC(ostream &aus);  //write the colors
  void dotout(ostream &aus);  //write the dot format
  void
  adjout(ostream &aus, char sep);  //write the adjacency list with seperator
  void writeadj(ostream &aus); //write the graph as a python adjacency graph


 private:

  int M;    //maximum number of vertices
  int V;    //number of vertices
  int E;    //number of edges
  set *nbr; //neighbor lists are stored as sets
  int *clr; //colors
  int nqual;              //number of qualities
  int *quality[MAXQ];     //pointers to quality vectors
  int nwgt;               //number of weights
  double *weights[MAXQ];  //pointers to weight vectors

};

#endif /* _SETU_H */
