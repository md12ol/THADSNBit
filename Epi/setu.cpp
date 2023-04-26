#include <iostream>
#include <fstream>
#include <random>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <random>
#include <chrono>
#include <bitset>

using namespace std;

#include "setu.h"

bool infected(unsigned long size, double alpha);
//fitness proportional selector used in simulations
int rselect(double *v, double ttl, int N) {

    double dart;   //selection variable
    int i;         //loop index

    dart = drand48() * ttl - v[0];  //throw dart
    i = 0;  //stat at the front
    while ((i < N) && (dart > 0))
        dart -= v[++i];  //figure out where the dart landed
    if (i >= N)i = N - 1; //stupid failsafe
    return (i); //say where the dart landed

}

set::set() {//default constructor

    max = n = 0;  //max==0 is the clue that the structure is unallocated
    mem = 0;

}

set::set(int *m, int z) {//construct with a list of elements

    create(m, z);

}

set::set(const set &other) {//copy constructor

    if (max == 0) {
        max = n = 0;
        mem = 0;
        return;
    }

    max = other.max;
    n = other.n;
    mem = new int[max];
    for (int i = 0; i < n; i++)mem[i] = other.mem[i];

}

set::~set() {//destructor

    destroy();

}

//utilities
void set::create(int *m, int z) {//create a set with a given membership and size

    for (int i = 0; i < z; i++)add(m[i]);  //put the elements in

}

void set::destroy() {//destroy a set

    if (max == 0)return;  //don't try to destroy empty structures

    delete[] mem;
    mem = 0;
    max = n = 0;

}

void set::setempty() {//mark as empty for mass allocation

    max = n = 0;
    mem = 0;

}

void set::copy(set &other) {//copy another set

    destroy();
    if (other.max == 0)return;
    max = other.max;
    n = other.n;
    mem = new int[max];
    for (int i = 0; i < n; i++)mem[i] = other.mem[i];

}

void set::copyO(set &other, int q) {//copy another set

    destroy();
    if (other.max == 0)return;
    max = other.max;
    n = other.n;
    mem = new int[max];
    for (int i = 0; i < n; i++)mem[i] = other.mem[i] + q;

}

void set::enlarge() {//increment max

    int i;
    int *nw;

    nw = new int[max + SETINCR];       //create new larger memory
    for (i = 0; i < n; i++)nw[i] = mem[i];  //transfer data
    delete[] mem;                 //delete old memory
    mem = nw;                        //install new memory
    max += SETINCR;                  //record new size

}

int set::add(int z) {//add a member, returns true if a not ALREADY

    int i, j;

    //cout << "ADD " << z << endl;
    //write(cout);

    if (max == 0) {//empty set, create everything
        max = SETINCR;
        n = 1;
        mem = new int[max];
        mem[0] = z;
        return (1);
    } else if (n == 0) {
        mem[0] = z;
        n = 1;
        return (1);
    } else {//existing set, see if the element is new
        for (i = 0; i < n; i++)if (mem[i] == z)return (0);
        if (n == max)enlarge();  //create more space if it is needed
        if (mem[n - 1] < z) {
            //cout << "add end" << endl;
            mem[n++] = z;
        } else {//ripple insert
            i = 0;
            while (z > mem[i])i++;
            for (j = n; j > i; j--)mem[j] = mem[j - 1];
            mem[i] = z;
            n++;
        }
    }
    return (1);
}

int set::remo(int z) {//remove a member, returns true if in there

    int i, j;

    //cout << "REMO " << z << endl;

    for (i = 0; i < n; i++)
        if (mem[i] == z) {//found it
            for (j = i; j < n; j++)mem[j] = mem[j + 1];//ripple delte
            n--;  //reduce set size
            return (1); //report succesful deletion
        }

    return (0);  //report vertex not found

}

void set::clear() {//clear the set of members

    n = 0;  //max the set empty

}

//information

int set::size() {//what is the size of the set

    return (n);

}

int set::memb(int z) {//is z a member? 0=no 1=yes

    int i;

    for (i = 0; i < n; i++) {
        if (mem[i] == z)return (1);
        if (mem[i] > z)return (0);
    }
    return (0);

}

int set::memz(int z) {//zth member

    if ((z < 0) || (z >= n))return (0);  //crock failsafe

    return (mem[z]);

}
//operations
void set::unyun(set &A, set &B) {//union

    int i;

    copy(A);
    for (i = 0; i < B.n; i++)add(B.mem[i]);

}

void set::inter(set &A, set &B) {//intersection

    int a, b;

    n = 0;  //erase current content
    a = 0;
    b = 0;
    do {
        if ((a >= A.n) || (b >= B.n))return; //done
        if (A.mem[a] == B.mem[b]) {
            add(A.mem[a]);
            a++;
            b++;
        } else if (A.mem[a] < B.mem[b])a++; else b++;
    } while (1);

}

void set::setdf(set &A, set &B) {//set difference

    int i;

    copy(A);
    for (i = 0; i < B.n; i++)remo(B.mem[i]);

}

void set::symdf(set &A, set &B) {//symmetric difference

    int i;

    for (i = 0; i < A.n; i++)if (B.memb(A.mem[i]) == 0)add(A.mem[i]);
    for (i = 0; i < B.n; i++)if (A.memb(B.mem[i]) == 0)add(B.mem[i]);

}

double set::sumAt(double *ft) {//sum indecies of ft in the set

    double ttl;   //summing register
    int i;        //loop index

    ttl = 0.0;                           //zero the register
    for (i = 0; i < n; i++)ttl += ft[mem[i]];   //sum the indexed members
    return (ttl);                       //return value

}

//ft should be strictly positive - it places a probability measure on the
//members of the set.
int set::FPS(double *ft) {//fitness proportional selection

    double ttl, dart;  //dart and dartboard
    int i;            //index

    if (n == 0)return (0);  //dumbass filter

    ttl = sumAt(ft);           //get the size of the dartboard
    if (ttl > 0.01) {
        dart = drand48() * ttl;      //throw the dart
        i = 0;
        dart -= ft[mem[0]];    //hit first?
        while ((dart > 0) && (i < n)) {  //iterated:
            i++;
            dart -= ft[mem[i]];  //hit next?
        }
        if (i == n)i = n - 1;           //failsafe
        return (mem[i]);               //return result
    } else return (mem[lrand48() % n]); //zero total?  select at radnom
}

int graph::RNB(int v) {//Random neighbor

    if ((v < 0) || (v >= V))return (-1);  //return "nope"
    if (degree(v) == 0)return (-1);   //return "nope"
    return (nbrmod(v, lrand48() % degree(v))); //return a random neighbor

}

//input-output
void set::write(ostream &aus) {//write set

    aus << n << " " << max << endl;
    for (int i = 0; i < n; i++)cout << mem[i] << endl;

}

void set::writememb(ostream &aus) { //write members on one line

    if (n == 0) {
        aus << endl;
        return; //nothing to write
    }

    aus << mem[0];
    for (int i = 1; i < n; i++)aus << "\t" << mem[i];
    aus << endl;

}

void set::writememb(ostream &aus, char sep) { //write members on one line

    if (n == 0) {
        aus << endl;
        return; //nothing to write
    }

    aus << mem[0];
    for (int i = 1; i < n; i++)aus << sep << mem[i];
    aus << endl;

}

void set::read(istream &inp) {//read set

    char buf[60];
    int k;

    destroy();
    inp.getline(buf, 59);
    n = atoi(buf);
    k = 0;
    while (buf[k] != ' ')k++;
    while (buf[k] == ' ')k++;
    max = atoi(buf + k);
    mem = new int[max];
    for (k = 0; k < n; k++) {
        inp.getline(buf, 59);
        mem[k] = atoi(buf);
    }
}

void set::readmemb(istream &inp) {//read members on one line

    char buf[10000];

    int k, l;

    destroy();
    inp.getline(buf, 999);
    l = strlen(buf);
    if (l > 0) {
        n = 0;
        for (k = 0; k < l; k++)if (buf[k] == ' ')n++;
        n++;
        max = n;
        mem = new int[max];
        mem[0] = atoi(buf);
        k = 0;
        l = 1;
        while (l < n) {
            while (buf[k] != ' ')k++;
            while (buf[k] == ' ')k++;
            mem[l++] = atoi(buf + k);
        }
    }
}

/******************************Graph CODE*******************************/



graph::graph() {//initialize an empty structure

    M = V = E = 0;  //M==0 isthe empty structure clue
    nbr = 0;    //nil nieghbor list
    clr = 0;    //nil color buffer
    nqual = nwgt = 0;  //record that no weights or qualities exist yet
    for (int i = 0; i < MAXQ; i++) {
        quality[i] = 0;
        weights[i] = 0;
    }  //zero out weights,quals

}

graph::graph(int max) {//initialize to maximum of M vertices

    M = 0;        //prevent pre-natal calls to destroy
    create(max);  //call the creation method

}

graph::~graph() {//delete s structure

    destroy();  //deallocate everything

}

//utilities
void graph::create(int max) {//create with max maximum vertices

    //cout << "Create " << max << endl;

    if (M != 0)destroy();  //clear the graph if it is not empty

    M = max;          //set the maximum number of vertices
    V = E = 0;          //indicate the graph starts empty
    nbr = new set[M]; //create set variables
    clr = 0;          //no colors
    nqual = nwgt = 0;   //record that no weights or qualities exist yet
    for (int i = 0; i < MAXQ; i++) {
        quality[i] = 0;
        weights[i] = 0;
    }  //zero out weights,quals

}

void graph::destroy() {//deallocate everything

    int i;  //loop index

    if (M > 0) {
        delete[] nbr;
        if (clr != 0)delete[] clr;
    }

    //clear everything as well
    for (i = 0; i < nqual; i++) {
        delete[] quality[i];
        quality[i] - 0;
    } //clear quals
    nqual = 0;                                                 //record no quals
    for (i = 0; i < nwgt; i++) {
        delete[] weights[i];
        weights[i] - 0;
    }  //clear weights
    nwgt =
        0;                                                  //record no weights
    M = V = E = 0;
    nbr = 0;
    clr = 0;

}

void graph::Enlarge(int newmax) { //increase maximum vertices to newmax

    set *nw, sw;  //new set list for new neighbors
    int *nc;     //expanded color list
    int i;       //loop index

    if (newmax <= M)return;   //dumbass request
    nw = new set[newmax];    //create new neighbor list
    if (clr != 0)nc = new int[newmax];  //create new color list if needed
    for (i = 0; i < newmax; i++)nw[i].setempty(); //mark as empty
    for (i = 0; i < V; i++) {
        nw[i].copy(nbr[i]);  //copy current neighbors
        if (clr != 0)nc[i] = clr[i];
    }
    for (i = 0; i < M; i++)nbr[i].destroy();  //deallocate nightbors
    delete[] nbr;  //deallocate neighbor list
    if (clr != 0) {//if color needs updating
        delete[] clr; //delete the old color
        clr = nc; //update
    }
    nbr = nw;    //assign new neighbor list
    M = newmax;  //update maximum

}

void graph::clearE() {//change graph to empty

    empty(V);

}

//Request a quality
int graph::RQquality(int v) {//index of quality or -1 if none left

    if (nqual == MAXQ)return (-1);    //Sorry -- no qualities left
    quality[nqual] = new int[V];    //create an integer register for each vertex
    for (int i = 0; i < V; i++)quality[nqual][i] = v; //initialize to v
    nqual++;                      //record that quality was allocated
    return (nqual - 1);              //return index of quality

}

int graph::RQquality(int *Q) {//as above, but initialize to Q -- size==V

    if (nqual == MAXQ)return (-1);    //Sorry -- no qualities left
    quality[nqual] = new int[V];    //create an integer register for each vertex
    for (int i = 0; i < V; i++)quality[nqual][i] = Q[i]; //initialize to Q
    nqual++;                      //record that quality was allocated
    return (nqual - 1);              //return index of quality

}

void graph::RecordQ(int num, int *Q) {//assign Q to the quality in question

    if ((num < 0) || (num >= nqual))
        return;  //try to initialize a non-existant quality
    for (int i = 0; i < V; i++)quality[num][i] = Q[i]; //set to Q

}

void graph::RetrieveQ(int num, int *Q) {//get the values in the quality -> Q

    if ((num < 0) || (num >= nqual))
        return;  //try to retrieve a non-existant quality
    for (int i = 0; i < V; i++)Q[i] = quality[num][i]; //dump values into Q

}

//request a weight
int graph::RQweight(double v) {//index of quality or -1 if none left

    if (nwgt == MAXQ)return (-1);     //Sorry -- no qualities left
    weights[nwgt] = new double[V];  //create an double register for each vertex
    for (int i = 0; i < V; i++)weights[nwgt][i] = v; //initialize to v
    nwgt++;                       //record that quality was allocated
    return (nwgt - 1);               //return index of quality

}

int graph::RQweight(double *W) {//as above but initialize to W -- size==V

    if (nwgt == MAXQ)return (-1);     //Sorry -- no qualities left
    weights[nwgt] = new double[V];  //create an double register for each vertex
    for (int i = 0; i < V; i++)weights[nwgt][i] = W[i]; //initialize to W
    nwgt++;                       //record that quality was allocated
    return (nwgt - 1);               //return index of quality

}

void graph::RecordW(int num, double *W) {//assign W to the weights in question

    if ((num < 0) || (num >= nwgt))
        return;  //try to retrieve a non-existant weight
    for (int i = 0; i < V; i++)weights[num][i] = W[i]; //set to W

}

void graph::RetrieveW(int num, double *W) {//get the values in the quality -> W

    if ((num < 0) || (num >= nwgt))
        return;  //try to retrieve a non-existant weight
    for (int i = 0; i < V; i++)W[i] = weights[num][i]; //transfer weights to W

}

//are you infected with n contacts of strength alpha
int graph::infected(int n, double alpha) {//SIR utility routine

    double beta;

    beta = 1 - exp(n * log(1 - alpha));
    if (drand48() < beta)return (1);
    return (0);

}

//initializers
void graph::empty(int n) {//empty graph

    int i;

    if ((M == 0) || (M < n))create(n);  //make sure storage is available
    V = n;  //vertices for an empty graph
    E = 0;  //edges for an empty graph
    for (i = 0; i < V; i++)nbr[i].clear();  //empty out the neghbor lists

}

void graph::Kn(int n) {//complete

    int i, j;

    if ((M == 0) || (M < n))create(n);  //make sure storage is available
    V = n;
    E = n * (n - 1) / 2;             //create values for V and E

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)if (j != i)nbr[i].add(j);
    }

}

void graph::Knm(int n, int m) {//complete bipartite

    int i, j;

    if ((M == 0) || (M < n + m))create(n + m);  //make sure storage is available
    V = n + m;
    E = n * m;                    //create values for V and E

    for (i = 0; i < n; i++)
        for (j = 0; j < m; j++) {//loop over relevant pairs
            nbr[i].add(n + j);
            nbr[n + j].add(i);
        }

}

void graph::Cn(int n) {//cycle

    int i;
    int ls[2];

    if ((M == 0) || (M < n))create(n);  //make sure storage is available
    V = n;
    E = n;                     //create values for V and E
    for (i = 0; i < n; i++) {
        ls[0] = (i + 1) % n;
        ls[1] = (i + n - 1) % n;
        nbr[i].create(ls, 2);
    }

}

void graph::Pn(int n, int m) {//Petersen n,m

    int i;

    if ((M == 0) || (M < 2 * n))create(2 * n); //make sure storage is available

    if (m > n)m %= n;  //failsafe the arithmetic

    V = 2 * n;
    E = 3 * n;  //addign the vertex and edge values
    for (i = 0; i < n; i++) {
        nbr[i].add((i + 1) % n);        //outer cycle left
        nbr[i].add((i - 1 + n) % n);      //outer cycle right
        nbr[i].add(i + n);            //spoke
        nbr[i + n].add(i);            //other end of spoke
        nbr[i + n].add((i + m) % n + n);    //inner cycle left
        nbr[i + n].add((i - m + n) % n + n);  //inner cycle right
    }

}

void graph::Hn(int dim) {//Hypercube

    int i, j, b;       //index variables, size buffer
    int bits[20];    //bit array
    int nb;          //neighbor buffer


    b = 1;                         //initialize size buffer
    for (i = 0; i < dim; i++) {//compute size and bit array
        bits[i] = b; //save bit
        b *= 2;      //compute size
    }
    if ((M == 0) || (M < b))create(b);  //make sure storage is available
    V = b;        //record number of vertices
    E = b * dim / 2;  //compute and record number of vertices
    for (i = 0; i < b; i++) {//loop over vertices
        for (j = 0; j < dim; j++) {//loop over neighbors
            nb = (i ^ bits[j]); //compute neighbor with bitwise xor
            nbr[i].add(nb);
        }
    }
}

void graph::RNGnm(int n, int m) {//Ring with +/-m neighbors

    int i, j;

    if ((M == 0) || (M < n))create(n); //make sure storage is available

    if (m > n)m %= n;  //failsafe the arithmetic

    V = n;
    E = m * n;  //adding the vertex and edge values
    for (i = 0; i < n; i++) {
        for (j = 1; j <= m; j++) {
            nbr[i].add((i + j) % n);
            nbr[i].add((i - j + n) % n);
        }
    }
}

/* the UTAM method ASSUMES that ed has length V(V-1)/2*/
void graph::UTAM(int *ed) {//initialize from an upper triangular adj. matrix

    int i, j, k;  //loop index variables

    clearE();
    V = M;
    k = 0;
    for (i = 0; i < V - 1; i++)
        for (j = i + 1; j < V; j++) {
            if (ed[k++] == 1) {
                nbr[i].add(j);
                nbr[j].add(i);
                E++;
            }
        }
}


//wk* holds the walk, wl is the size of the array */
void graph::WalkO(int *wk, int wl) {//overlaying walk representation

    int i;

    clearE();
    for (i = 0; i < wl - 1; i++)
        if (edgeP(wk[i], wk[i + 1]) == 0)
            toggle(wk[i],
                   wk[i + 1]);

}
/* holds the walk, wl is the size of the array */
void graph::WalkT(int *wk, int wl) {//toggling walk representation

    int i;

    clearE();
    for (i = 0; i < wl - 1; i++)toggle(wk[i], wk[i + 1]);

}

//el* holds the edges, ne is the size of the array */
void graph::EdgeLst(int *el, int ne) {//Edge list

    int i;

    clearE();
    for (i = 0; i < ne - 1; i += 2)toggle(el[i], el[i + 1]);

}
/********************************************************************/
/*
 * The following routine has been superceded by running individual
 * commands in the calling software because we need to manage
 * operation probabilities.
 *
 */

//This routine takes a vector of triples of integers The first is the
//command A(dd)=0 D(elete)=1 T(oggle)=2 S(wap)=3 the second specifies
//the first argument of the command the third specifies the second.
//L is the length of the gene
void graph::ADTS(int **cs, int L) {//implement an add delete toggle swap

    int c;     //command index


    if (V < 2)return; //nothing to do
    for (c = 0; c < L; c++) {//loop over commands
        switch (cs[c][0]) {
            case 0: //Add
                add(cs[c][1] % V, cs[c][2] % V);
                //cout << "ADD" << endl;
                break;
            case 1: //Delete
                del(cs[c][1] % V, cs[c][2] % V);
                //cout << "DEL" << endl;
                break;
            case 2: //Toggle
                toggle(cs[c][1] % V, cs[c][2] % V);
                //cout << "TOG" << endl;
                break;
            case 3: //edge swap with degree bound two
                edgeswap(cs[c][1], cs[c][2], 2);
                //cout << "SWP" << endl;
                break;
        }
    }
}

//This routine takes a vector of triples of integers The first is the
//command H(op)=0 A(dd)=1 D(elete)=2 T(oggle)=3 S(wap)=4 the second specifies
//the first argument of the command the third specifies the second.
//L is the length of the gene
void graph::HADTS(int **cs, int L) {//implement an add delete toggle swap

    int c;     //command index


    if (V < 2)return; //nothing to do
    for (c = 0; c < L; c++) {//loop over commands
        switch (cs[c][0] % 5) {
            case 0: //Hop
                hop(cs[c][1] % V,
                    cs[c][2] % V,
                    cs[c][2] / V);  //assume large integer rep
                break;
            case 1: //Add
                add(cs[c][1] % V, cs[c][2] % V);
                //cout << "ADD" << endl;
                break;
            case 2: //Delete
                del(cs[c][1] % V, cs[c][2] % V);
                //cout << "DEL" << endl;
                break;
            case 3: //Toggle
                toggle(cs[c][1] % V, cs[c][2] % V);
                //cout << "TOG" << endl;
                break;
            case 4: //edge swap with degree bound two
                edgeswap(cs[c][1], cs[c][2], 2);
                //cout << "SWP" << endl;
                break;
        }
    }
}

void graph::copy(graph &other) {//copy another graph

    int i, m;    //loop index variable

    if (other.V == 0) {//clear the graph to "copy" and empty graph
        clearE();
        return;
    }

    if (other.V > M) {//make surthe graphi s big enough
        Enlarge(other.V);
    }

    V = other.V;
    E = other.E;

    for (i = 0; i < other.V; i++) {//copy the other vertex
        nbr[i].copy(other.nbr[i]);
    }

}

//Random graph generators -- new section Nov 2 2020
void graph::BA(int n, int m) {//Barabasi-Albert graph with n vertices;
    //m edges added at each step

    int i, j;         //loop indices
    vector<int> p;   //vertices repeated to enhance probability of being chosen
    set temp;        //vertices already chosen
    int dart, hold;   //select randomly
    int dx[n];       //shuffled vertices


    if ((M == 0) || (M < n))create(n); //make sure storage is available

    V = n;
    E = 1; //assign the vertex and edge values
    for (i = 0; i < n; i++) dx[i] = i; //initialize shuffle
    for (i = 0; i < n - 2; i++) {//shuffle
        dart = lrand48() % (n - i) + i;
        hold = dx[dart];
        dx[dart] = dx[i];
        dx[i] = hold;
    }

    //for(i=0;i<n;i++) cout << dx[i] << " "; cout << endl;
    p.push_back(0); //initialize prob vector

    for (i = 1; i < n; i++) { //add vertices one by one
        temp.setempty();
        for (j = 0; j < min(m, i); j++) { //add m edges to existing vertices
            do { //choose from prob vector
                dart = lrand48() % p.size();
            } while (temp.memb(p[dart]));
            add(dx[i], dx[p[dart]]); //add an edge
            E++; //increment count
            temp.add(p[dart]); //add new edge to temp
        }
        for (j = 0; j < temp.size(); j++)
            p.push_back(temp.memz(j)); //add new connections
        p.push_back(i); //add new vertex to prob vector
        //for(j=0;j<p.size();j++) cout << p[j] << " "; cout << endl;
    }
}

//Power law clustering graph with n vertices; m edges added at each
//step; p is the chance of adding an extra edge to create a triangle
//(modification of BA graph to enhance clustering coeff of vertices)
void graph::PCG(int n, int m, double prob) {

    int i, j, k;
    vector<int> p; //vertices repeated to enhance probability of being chosen
    vector<int> pp; //holding vector
    set temp; //vertices already chosen
    set neigh; //neighbours of target vertex
    int dart, hold; //select randomly
    int dx[n], ix[n]; //shuffled vertices
    int nn; //neigbour vertex
    int target; //vertex to connect to

    if ((M == 0) || (M < n))create(n); //make sure storage is available

    V = n;
    E = 0; //assign the vertex and edge values

//shuffle vertices to ensure randomness
    for (i = 0; i < n; i++) dx[i] = i; //initialize shuffle
    for (i = 0; i < n - 2; i++) {//shuffle
        dart = lrand48() % (n - i) + i;
        hold = dx[dart];
        dx[dart] = dx[i];
        dx[i] = hold;
    }
//store dx inverse
    for (i = 0; i < n; i++) ix[dx[i]] = i;

//for(i=0;i<n;i++) cout << dx[i] << " "; cout << endl;
    p.push_back(0); //initialize prob vector

    for (i = 1; i < n; i++) { //add vertices one by one
        temp.setempty(); //keep track of added edges
        pp.clear();
        std::copy(p.begin(),
                  p.end(),
                  back_inserter(pp)); //copy prob vector to scratch vector
        for (j = 0; j < min(m, i); j++) { //add m edges to existing vertices
            for (k = 0; k < pp.size(); k++)
                if (edgeP(dx[i], dx[pp[k]]))
                    pp.erase(pp.begin() + k);//delete existing edges from pp
            if (pp.size() > 0) { //if there are possible edges to add
                dart = lrand48() % pp.size();  //choose from prob vector
                target = pp[dart];
                if (i > m && drand48() < prob) {//make a triangle
                    //check if possible
                    if (degree(dx[target]) > 0) {
                        neigh.setempty();
                        for (k = 0; k < degree(dx[target]); k++)
                            neigh.add(nbrmod(dx[target], k));
                        do {
                            dart = lrand48() % neigh.size();
                            nn = neigh.memz(dart);
                            neigh.remo(nn);
                        } while (edgeP(nn, dx[i]) && neigh.size() > 0);

                        if (!edgeP(nn,
                                   dx[i])) { //add edge that creates triangle
                            add(nn, dx[i]);
                            temp.add(i); //add vertex to temp
                            temp.add(ix[nn]); //add new edge to temp
                        }
                    }
                }
                add(dx[i], dx[target]); //add an edge to target
                temp.add(target); //add target to temp
                temp.add(i); //add new vertex to temp
            }
        }
        for (j = 0; j < temp.size(); j++)
            p.push_back(temp.memz(j)); //add new vertices to prob vector
    }
}

void graph::ER(int n, double p) {//Erdo-Renyi graph with n vertices; prob p
//makes a connected graph -- will not work well if p*(n-1)<1

    int i, j; //loop index

    if ((M == 0) || (M < n))create(n); //make sure storage is available
    V = n;
    E = 0; //start with n vertices and no edges

    do {
        for (i = 0; i < n - 1; i++)
            for (j = i + 1; j < n; j++) {//cycle through all possible edges
                if (drand48() < p) {//add edge with prob p
                    add(i, j);
                }
            }
    } while (!connectedP());
}

//This sort of random graph has n nodes connected to k nearest
//neighbours with edges rewired with probability p
void graph::WS(int n, int k, double p) {//Watts-Strogatz graph

    int i, j; //loop index
    int x; //random node
    int nb[n][k]; //neighbours
    int deg; //degree of node

    deg = k - k % 2; //subtract one if odd
    RNGnm(n, deg / 2);

    for (i = 0; i < n; i++)Nbrs(i, nb[i]); //save neighbours

    do {
        for (i = 0; i < n; i++) {
            for (j = 0; j < deg; j++) {
                if (drand48() < p) {//swap with prob p
                    do {//choose another node
                        x = lrand48() % n;
                    } while (x == i || edgeP(i, x));
                    add(i, x);
                    del(i, nb[i][j]);
                }
            }
        }
    } while (!connectedP());

}


//This sort of random graph has n nodes connected to k nearest
//neighbours with edges added with probability p
void graph::NWS(int n, int k, double p) {//Newman-Watts-Strogatz graph

    int i, j; //loop index
    int x; //random node
    int deg; //degree of node

    deg = k - k % 2; //subtract one if odd
    RNGnm(n, deg / 2);

    for (i = 0; i < n; i++) {
        for (j = 0; j < deg; j++) {
            if (drand48() < p) {//add edge with prob p
                do {//choose another node
                    x = lrand48() % n;
                } while (x == i || edgeP(i, x));
                add(i, x);
            }
        }
    }
}


//modifiers
void graph::add(int a, int b) {//force an edge to add

    if ((a < 0) || (a >= V))a = ((a % V) + V) % V;  //failsafe vertex identity a
    if ((b < 0) || (b >= V))b = ((b % V) + V) % V;  //failsafe vertex identity b
    if (a == b)return;  //enforce simplicity
    if (nbr[a].memb(b) == 0)toggle(a, b);

}

void graph::ladd(int v, int n1, int n2) {//hop an edge

    int nb1, nb2;  //identity of new and old neighbors
    int d1, d2;    //degrees of neighbors

    if (M == 0)return;  //empty graph, nothing to do
    //cout << "Not empty" << endl;
    if ((v < 0) || (v >= V))v = (v % V + V) % V; //force correct vertex number
    //cout << "Vertex is " << v << endl;
    d1 = degree(v);   //retrieve the degree
    //cout << "Its degree is " << d1 << endl;
    if (d1 < 1)return; //cannot hop with no neighbors
    if ((n1 < 0) || (n1 >= d1))
        n1 = (n1 % d1 + d1) % d1;  //force possible neighbor
    nb1 = nbrmod(v, n1);  //retrieve nieghbor
    //cout << "Neighbor is " << nb1 << endl;
    d2 = degree(nb1);    //get degree of neighbor
    if (d2 < 2)return;    //no place to hop
    //cout << "Its degree is " << d2 << endl;
    if ((n2 < 0) || (n2 >= d2))
        n2 = (n2 % d2 + d2) % d2;  //force possible neighbor, again
    nb2 = nbrmod(nb1, n2);  //retrieve nieghbor-squared
    //cout << v << " " << nb1 << " " << nb2 << endl;
    if (edgeP(v, nb2) == 1)return;  //no hop possible - its a triangle
    if (v == nb2)return; //trying to add a loop
    //cout << "Inserting " << v << " " << nb2 << endl;
    add(v, nb2); //add the new edge

}

void graph::del(int a, int b) {//force an edge to be gone

    if ((a < 0) || (a >= V))a = ((a % V) + V) % V;  //failsafe vertex identity a
    if ((b < 0) || (b >= V))b = ((b % V) + V) % V;  //failsafe vertex identity b
    if (a == b)return;  //enforce simplicity
    if (nbr[a].memb(b) == 1)toggle(a, b);

}

void graph::ldel(int v, int n1, int n2) {//hop an edge

    int nb1, nb2;  //identity of new and old neighbors
    int d1, d2;    //degrees of neighbors

    if (M == 0)return;  //empty graph, nothing to do
    //cout << "Not empty" << endl;
    if ((v < 0) || (v >= V))v = (v % V + V) % V; //force correct vertex number
    //cout << "Vertex is " << v << endl;
    d1 = degree(v);   //retrieve the degree
    //cout << "Its degree is " << d1 << endl;
    if (d1 < 1)return; //cannot hop with no neighbors
    if ((n1 < 0) || (n1 >= d1))
        n1 = (n1 % d1 + d1) % d1;  //force possible neighbor
    nb1 = nbrmod(v, n1);  //retrieve nieghbor
    //cout << "Neighbor is " << nb1 << endl;
    d2 = degree(nb1);    //get degree of neighbor
    if (d2 < 2)return;    //no place to hop
    //cout << "Its degree is " << d2 << endl;
    if ((n2 < 0) || (n2 >= d2))
        n2 = (n2 % d2 + d2) % d2;  //force possible neighbor, again
    nb2 = nbrmod(nb1, n2);  //retrieve nieghbor-squared
    //cout << v << " " << nb1 << " " << nb2 << endl;
    if (edgeP(v, nb2) == 1)return;  //no hop possible - its a triangle
    if (v == nb2)return; //trying to add a loop
    //cout << "Deleting " << v << " " << nb2 << endl;
    del(v, nb2); //delete the new edge
}

void graph::toggle(int a, int b) {//toggle an edge

    int u, v;

    if (M == 0)return;
    if ((a < 0) || (a >= V))a = ((a % V) + V) % V;  //failsafe vertex identity a
    if ((b < 0) || (b >= V))b = ((b % V) + V) % V;  //failsafe vertex identity b
    if (a == b)return;  //enforce simplicity
    if (nbr[a].memb(b) == 1) {//edge exists, turn off
        //cout << "Toggle " << a << " " << b << " off." << endl;
        u = nbr[a].remo(b);
        v = nbr[b].remo(a);
        //if(u!=v)cout << u << " " << v << endl;
        E--;
    } else {//edge does not exist, turn on
        //cout << "Toggle " << a << " " << b << " on." << endl;
        u = nbr[a].add(b);
        v = nbr[b].add(a);
        //if(u!=v)cout << u << " " << v << endl;
        E++;
    }

}

void graph::loggle(int v, int n1, int n2) {//hop an edge

    int nb1, nb2;  //identity of new and old neighbors
    int d1, d2;    //degrees of neighbors

    if (M == 0)return;  //empty graph, nothing to do
    //cout << "Not empty" << endl;
    if ((v < 0) || (v >= V))v = (v % V + V) % V; //force correct vertex number
    //cout << "Vertex is " << v << endl;
    d1 = degree(v);   //retrieve the degree
    //cout << "Its degree is " << d1 << endl;
    if (d1 < 1)return; //cannot hop with no neighbors
    if ((n1 < 0) || (n1 >= d1))
        n1 = (n1 % d1 + d1) % d1;  //force possible neighbor
    nb1 = nbrmod(v, n1);  //retrieve nieghbor
    //cout << "Neighbor is " << nb1 << endl;
    d2 = degree(nb1);    //get degree of neighbor
    if (d2 < 2)return;    //no place to hop
    //cout << "Its degree is " << d2 << endl;
    if ((n2 < 0) || (n2 >= d2))
        n2 = (n2 % d2 + d2) % d2;  //force possible neighbor, again
    nb2 = nbrmod(nb1, n2);  //retrieve nieghbor-squared
    //cout << v << " " << nb1 << " " << nb2 << endl;
    if (edgeP(v, nb2) == 1)return;  //no hop possible - its a triangle
    if (v == nb2)return; //trying to add a loop
    //cout << "Inserting " << v << " " << nb2 << endl;
    toggle(v, nb2); //toggle the new edge

}

void graph::simplexify(int a) {//simplexify at a

    if ((a < 0) || (a >= V))return; //don't attempt the impossible

    int qq;  //size of new simplex
    int tt;  //lower bound on new clique vertices
    int ss;  //saved neighbor of a

    qq = nbr[a].size(); //get new simplex size
    if (qq <= 1)return;  //nothing happens in this case

    int i, j; //loop index variables
    int *temp;

    if (M - V < qq)Enlarge(M + qq);  //if more space is needed, make it

    tt = V;       //record lower bound on new vertices
    V += (qq - 1);  //indicate that the new vertices exist
    temp = new int[qq];
    for (i = 0; i < qq; i++)temp[i] = nbr[a].memz(i); //save neighbots of i

    //create the cliqe
    for (i = tt; i < tt + qq - 1; i++)nbr[i].add(a);
    for (i = tt; i < tt + qq - 1; i++)
        for (j = tt; j < tt + qq - 1; j++)
            if (i != j)nbr[i].add(j);
    for (i = 1; i < qq; i++) {//move edges from neighbors of a to clique
        nbr[i + tt - 1].add(temp[i]);
        nbr[temp[i]].remo(a);
        nbr[temp[i]].add(tt + i - 1);
        nbr[a].remo(temp[i]);
    }
    //add the new vertices to a's list
    for (i = tt; i < tt + qq - 1; i++)nbr[a].add(i);
    delete[] temp;
    /*DONE*/
}

void graph::hop(int v, int n1, int n2) {//hop an edge

    int nb1, nb2;  //identity of new and old neighbors
    int d1, d2;    //degrees of neighbors

    if (M == 0)return;  //empty graph, nothing to do
    //cout << "Not empty" << endl;
    if ((v < 0) || (v >= V))v = (v % V + V) % V; //force correct vertex number
    //cout << "Vertex is " << v << endl;
    d1 = degree(v);   //retrieve the degree
    //cout << "Its degree is " << d1 << endl;
    if (d1 < 1)return; //cannot hop with no neighbors
    if ((n1 < 0) || (n1 >= d1))
        n1 = (n1 % d1 + d1) % d1;  //force possible neighbor
    nb1 = nbrmod(v, n1);  //retrieve nieghbor
    //cout << "Neighbor is " << nb1 << endl;
    d2 = degree(nb1);    //get degree of neighbor
    if (d2 < 2)return;    //no place to hop
    //cout << "Its degree is " << d2 << endl;
    if ((n2 < 0) || (n2 >= d2))
        n2 = (n2 % d2 + d2) % d2;  //force possible neighbor, again
    nb2 = nbrmod(nb1, n2);  //retrieve nieghbor-squared
    //cout << v << " " << nb1 << " " << nb2 << endl;
    if (edgeP(v, nb2) == 1)return;  //no hop possible - its a triangle
    if (v == nb2)return; //trying to add a loop
    //cout << "Deleting " << v << " " << nb1 << endl;
    del(v, nb1); //delete the old edge
    //cout << "Inserting " << v << " " << nb2 << endl;
    add(v, nb2); //add the new edge

}

//reject swaps if either vertex is too low degee or if there are more
//than two edges in the quartet
void graph::edgeswap(int a, int b, int k) {//decode and perform an edge swap
    //with degree bound k

    if (V < 4)return;  //we need to have four vertices

    int v1, v2;  //decoded vertices
    int n1, n2;  //decoded neighbors

    v1 = a % V;
    if (nbr[v1].size() <= k)return; //check degree bound, first vertex
    v2 = b % V;
    if (nbr[v2].size() <= k)return; //check degree bound, second vertex
    if (nbr[v1].memb(v2))return; //added edge in quartet
    n1 = (a / V) % nbr[v1].size(); //find first neighbor's index
    n1 = nbr[v1].memz(n1);     //acquire first neighbor
    n2 = (b / V) % nbr[v2].size(); //find second neighbor's index
    n2 = nbr[v2].memz(n2);     //acquire second neighbor
    if (nbr[v1].memb(n2))return; //added edge in quartet
    if (nbr[v2].memb(n1))return; //added edge in quartet
    if (nbr[n1].memb(n2))return; //added edge in quartet
    /****************ACTUALLY POSSIBLE TO SWAP*****************************/
    //cout << v1 << " " << n1 << " " << v2 << " " << n2 << endl;
    toggle(v1, n1);
    toggle(v2, n2);
    toggle(v1, n2);   //Perform the edge swap
    toggle(v2, n1);

}

//The absorb method adds a copy of a graph (other) to the current graph
void graph::Absorb(graph &other) {//add a copy of other to yourself

    int ofs;  //The unionizing offset
    int i, m;    //loop index variable

    if (other.V == 0)return; //no graph, nothing to add
    if (V + other.V > M) {//make sure the graph is large enough
        Enlarge(M + other.V);
        //cout << "Enlarge to " << M+other.V << endl;
    }

    ofs = V;         //save current number of vertices as offset
    V += other.V;    //compute new number of vertices
    E += other.E;    //compute new number of edges

    //cout << V << " " << E << "-ck" << endl;

    //write(cout);

    for (i = 0; i < other.V; i++) {//copy the other vertex
        m = i + ofs; //compute new vertex index
        //cout << "Upcopy " << m << " From " << i << endl;
        nbr[m].copyO(other.nbr[i], ofs); //shift copy the other graph
    }

    //write(cout);

}

void graph::Prism() {//create the prism of a graph

    int i;

    if (V == 0)return; //nothing to prismate
    if (2 * V > M)Enlarge(2 * V); //double the size
    for (i = 0; i < V; i++) {//loop over current vertices
        nbr[V + i].copyO(nbr[i], V); //create second copy of graph
    }
    for (i = 0; i < V; i++) {//now add in the prisim spokes
        nbr[i].add(i + V);  //up
        nbr[i + V].add(i);  //down
    }
    E = 2 * E + V; //update number of vertices
    V = 2 * V;   //update number of vertices
}

//information
int graph::size() {//number of vertices

    return (V);

}

int graph::edges() {//number of edges

    return (E);

}

int graph::edgeP(int a, int b) {//

    if ((a < 0) || (b < 0) || (a >= V) || (b >= V))
        return (0);  //non-vert are non-edge
    if (nbr[a].memb(b) == 1)return (1); else return (0); //check for exdge

}

void graph::dfrom(int z, int *ds) {//distances from z

    int i, j;   //loop index variables
    int fl;    //flag
    int cd;    //current distance
    int qq, vv; //scratch variables

    if ((z < 0) || (z >= V)) {//safety first
        cout << "Distance from invalid vertex requested!";
        return;
    }
    for (i = 0; i < V; i++)
        ds[i] = -1;  //negative one is the surrogate for infinity
    ds[z] = 0;
    cd = 0;
    do {
        fl = 0;
        for (i = 0; i < V; i++)
            if (ds[i] == cd) {//if we are at current distance
                qq =
                    nbr[i].size(); //get the size to prevent repeated function calls
                for (j = 0; j < qq; j++) {//loop over neighbors
                    vv = nbr[i].memz(j);  //get the neighbor
                    if (ds[vv] == -1) {//if we have not been here yet...
                        ds[vv] = cd + 1;  //assign the distance
                        fl = 1;         //and set the flag
                    }
                }
            }
        cd++;
    } while (fl == 1);  //until done
}

void graph::dfrom0(int z, int *ds) {//distances from z in color zero

    int i, j;   //loop index variables
    int fl;    //flag
    int cd;    //current distance
    int qq, vv; //scratch variables

    if ((z < 0) || (z >= V)) {//safety first
        cout << "Distance from invalid vertex requested!";
        return;
    }

    for (i = 0; i < V; i++)
        ds[i] = -1;  //negative one is the surrogate for infinity
    ds[z] = 0;
    cd = 0;
    do {
        fl = 0;
        for (i = 0; i < V; i++)
            if (ds[i] == cd) {//if we are at current distance
                qq =
                    nbr[i].size(); //get the size to prevent repeated function calls
                for (j = 0; j < qq; j++) {//loop over neighbors
                    vv = nbr[i].memz(j);  //get the neighbor
                    if (clr[vv] == 0) {//if the neighbor is color zero
                        if (ds[vv] == -1) {//if we have not been here yet...
                            ds[vv] = cd + 1;  //assign the distance
                            fl = 1;         //and set the flag
                        }
                    }
                }
            }
        cd++;
    } while (fl == 1);  //until done
}

int graph::ecc(int z) {//eccentricity of a vertex

    int *q;      //distance table
    int i, rv;    //loop index, return value
    int b;       //origin loop index

    if (M == 0)return (0);  //empty graphs yields no eccentricity
    if (V == 0)return (0);  //like it says
    q = new int[V];       //create distance array
    dfrom(z, q);         //compute distances from 0
    rv = 0;               //initialize the eccentricity
    for (i = 0; i < V; i++) {//loop over the vertices
        if (q[i] == -1) {
            delete[] q;
            return (-1);
        } //the graph is disconnected
        if (q[i] > rv)rv = q[i];  //find the maximum
    }
    delete[] q;        //delete distance array
    return (rv);         //return the eccentricity

}

int graph::diameter() {//eccentricity of a vertex

    int rv;  //return value
    int i;   //index variable
    int cv;  //comparison value

    rv = ecc(0); //initialize return value
    if (rv == -1)return (-1);  //disconnected graph
    for (i = 1; i < V; i++) {//Find the maximum eccentricity
        cv = ecc(i);         //compute
        if (cv > rv)rv = cv;    //update
    }
    return (rv);  //return diameter
}

int graph::radius() {//eccentricity of a vertex

    int rv;  //return value
    int i;   //index variable
    int cv;  //comparison value

    rv = ecc(0); //initialize return value
    if (rv == -1)return (-1);  //graph is disconnected
    for (i = 1; i < V; i++) {//Find the maximum eccentricity
        cv = ecc(i);          //compute
        if (cv < rv)rv = cv;     //update
    }
    return (rv);  //return diameter

}

int graph::connectedP() {//is the graph connected?

    int *q;
    int i, rv;

    if (M == 0)return (1);  //empty graphs are vacuously connected
    if (V == 0)return (1);  //like it says
    q = new int[V];       //create distance array
    dfrom(0, q);         //compute distances from 0
    rv = 1;               //initial hypothesis - connected
    for (i = 0; (i < V) && (rv == 1); i++)
        if (q[i] == -1)
            rv = 0;  //evidence of disconnection
    delete[] q;        //delete distance array
    return (rv);

}

void graph::eccSeq(int *ecs) {//compute the eccentricity sequence

    if ((M == 0) || (V == 0)) {//failsafe empty graphs
        ecs[0] = 0;
        return;
    }

    int i, j; //loop index variables
    int *ds;

    if (connectedP() == 1) {
        ds = new int[V];      //create distance buffer
        for (i = 0; i < V; i++) {   //loop over vertices
            ecs[i] = 0;         //zero eccentricity
            dfrom(i, ds);      //find distances
            for (j = 0; j < V; j++)
                if (ds[j] > ecs[i])
                    ecs[i] = ds[j]; //find eccentricity
        }
        delete[] ds;
    } else for (i = 0; i < V; i++)ecs[i] = -1; //all distances infinite

}

int graph::nbrmod(int v, int n) {//compute the n%degreeth neighbor of v

    if ((v >= V) || (v < 0))return (0);  //return zero for stupid request

    return (nbr[v].memz(n % nbr[v].size()));

}

int graph::degree(int v) {//report the degree of v

    if ((v >= V) || (v < 0))return (0);  //return zero for stupid request

    return (nbr[v].size());

}

double graph::meandegree() {//report the mean degree of the graph

    double accu;

    accu = 0.0;
    for (int i = 0; i < V; i++)accu += nbr[i].size();
    return (accu / ((double) V));

}

int graph::Nbrs(int v, int *nb) {//report the neighboors of v

    int q, i;   //degree and loop index

    if ((v >= V) || (v < 0))return (0);  //return zero for stupid request

    q = degree(v);
    for (i = 0; i < q; i++)nb[i] = nbr[v].memz(i);
    return (q);

}

int graph::MaxCol() {//report the maximal color

    if ((M == 0) || (clr == 0))return (0); //no colors, return zero value

    int i, b;   //loop index and maxcolor

    b = clr[0];  //initialize max color
    for (i = 1; i < V; i++)if (clr[i] > b)b = clr[i]; //scan for max color

    return (b);

}

void
graph::DiffChar(int v, double omega, double *dc) {//diffusion character at v

    double *sm;    //summation buffer
    int i, j, t;     //index variables
    double curd;   //current degree, real
    int k;         //current degree, integer

    if (V == 0)return;          //don't be silly
    if ((v < 0) || (v >= V))return; //also do not be silly
    sm = new double[V];  //allocate the sum buffer
    for (i = 0; i < V; i++)dc[i] = 0.0;  //zero the initial diffusion character
    for (t = 0; t < 6 * V; t++) {//loop over time steps
        for (i = 0; i < V; i++)sm[i] = 0.0;  //zero the summation buffer
        dc[v] += 1.0;  //add gas to the center vertex
        for (i = 0; i < V; i++) {//loop over vertives
            k = degree(i);           //get the degree
            curd = ((double) k) + 1.0;  //make the divisor
            for (j = 0; j < k; j++)
                sm[nbrmod(i, j)] += dc[i] / curd;  //move the gas
            sm[i] += dc[i]
                / curd;                            //and a little gas stays
        }
        for (i = 0; i < V; i++)dc[i] = omega * sm[i]; //apply the decay term
        //for(i=0;i<V;i++)cout << dc[i] << " ";cout << endl;
    }
    delete[] sm;  //give back summation butter

}

//Quality and weight handling

void graph::RecordQ(int num, int dex, int val) {//Q[num][dex]=val

    if ((num >= 0) && (num < MAXQ)
        && (quality[num] != 0)) {//if the quality is there
        if ((dex >= 0) && (dex < V)) {//if vertex is there
            quality[num][dex] = val;                  //record the value
        }
    }
}

void graph::RecordW(int num, int dex, double val) {//W[num][dex]=val

    if ((num >= 0) && (num < MAXQ)
        && (weights[num] != 0)) {//if the quality is there
        if ((dex >= 0) && (dex < V)) {//if vertex is there
            weights[num][dex] = val;                  //record the value
        }
    }
}

int graph::RetrieveQ(int num, int dex) {//return Q[num][dex]

    if ((num >= 0) && (num < MAXQ)
        && (quality[num] != 0)) {//if the quality is there
        if ((dex >= 0) && (dex < V)) {//if vertex is there
            return (quality[num][dex]);
        }
    }
    return -1;
}

double graph::RetrieveW(int num, int dex) {//return Q[num][dex]

    if ((num >= 0) && (num < MAXQ)
        && (weights[num] != 0)) {//if the quality is there
        if ((dex >= 0) && (dex < V)) {//if vertex is there
            return (weights[num][dex]);
        }
    }
    return -1.0;
}

//Simulation methods
/*Run and SIR epidemic with a given patient zero; return maximum number
 *of people infected, length of epidemic, total number infected
 *
 *  Note that the method uses the color buffer
 *  0=S 1=I 2=R 3=Exposed
 *
 *
 */
void graph::SIR(int p0, int &max, int &len, int &ttl, double alpha,
                vector<int> &prof) {
    //Sir Method

    int NI;     //number of infected individuals
    int i, j, k;  //loop index variables
    int *nin;   //number of infected neighbors

    max = len = ttl = 0; //zero the reporting statistics
    if ((V == 0) || (p0 < 0) || (p0 >= V))return; //no one was infected...

    nin = new int[V]; //create infected neioghbor counters
    setC2(0);    //set the population to susceptible
    clr[p0] = 1;   //infect patient zero
    NI = 1;        //initialize to one person currently infected
    len = 0;       //initialize length variable
    while (NI > 0) {//while there is still an epidemic
        //cout << "LEN=" << len << " NI=" << NI << endl;
        for (i = 0; i < V; i++)
            nin[i] = 0; //zero the number of infected neighbors buffer
        for (i = 0; i < V; i++)
            if (clr[i] == 1) {//found infected individual
                for (j = 0; j < nbr[i].size(); j++)
                    nin[nbr[i].memz(j)]++; //record exposure
            }
        //check for transmission
        for (i = 0; i < V; i++)
            if ((clr[i] == 0) && (nin[i] > 0))
                clr[i] = 3 * infected(nin[i], alpha);
        if (NI > max)max = NI; //check for updated maximum
        ttl += NI; //add the infected to the total
        NI = 0;  //zero the number infected counter
        for (i = 0; i < V; i++)
            switch (clr[i]) {//status update
                case 0: //susceptible, do nothing
                    break;
                case 1: //infected, move to removed
                    clr[i] = 2;
                    break;
                case 2: //removed, do nothing
                    break;
                case 3: //newly infected
                    clr[i] = 1;
                    NI++;
                    break;
            }
        len++; //record the time step
        prof.push_back(NI);
    }
    delete[] nin;  //return storage for nin buffer
}

void graph::create_new_variant(bitset<dna_len> &orig, bitset<dna_len> &variant, vector<int> rv, int low_bnd, int up_bnd) {
    // Randomness
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    shuffle(rv.begin(), rv.end(), default_random_engine(seed));
    int num_edits = (int) lrand48() % (up_bnd - low_bnd) + low_bnd;

    // Copy
    variant = orig;

    // Random edits
    for (int i = 0; i < num_edits; i++) {
        variant.flip(rv.at(i));
    }
}

struct sort_pred {
  bool operator()(const std::pair<double, pair<int, bitset<dna_len>>> &left,
                  const std::pair<double, pair<int, bitset<dna_len>>> &right) {
      return left.first < right.first;
  }
};

void graph::varInfected(bitset<dna_len> &immunity, vector<pair<int, bitset<dna_len>>> &strains, double alpha, int &str_id) {
    int one_cnt;
    int potent_ones;
    int max_ones = 0;
    int worst_str = 0;
    immunity.flip();
    vector<pair<double, pair<int, bitset<dna_len>>>> serv_strains;
    // Find worst strain
    for (auto &strain: strains) {
        one_cnt = (int) strain.second.count();
        bitset<dna_len> tmp = strain.second & immunity;
        potent_ones = (int) tmp.count();
        serv_strains.emplace_back((double) potent_ones / one_cnt, (strain));
    }
    immunity.flip();

    sort(serv_strains.rbegin(), serv_strains.rend(), sort_pred());
    for (auto &strain: serv_strains) {
        if (infected(1, strain.first * alpha) == 1) {
            str_id = strain.second.first;
            return (void) "DONE";
        }
    }
}

//Sir w Variants
void graph::varSIR(int p0, int &vcnt, vector<int> vprofs[], bitset<dna_len> variants[], int vorigs[], pair<int, int> vtimes[], double alpha, int lB, int uB, double var_prob) {
    if ((V == 0) || (p0 < 0) || (p0 >= V))return;

    int NI;     //number of infected individuals
    int new_inf[max_vars];
    bitset<dna_len> immunity[V];
    int inf_with[V];
    vector<int> randVec;
    int day = 0;
    vector<vector<int>> inf_neigh; // Nodes that are adjacent and can infect
    // Strins adjacent to susceptible node
    vector<pair<int, bitset<dna_len>>> strains;
    int strain_id;

    // Setup immunity and mde_var
    for (int i = 0; i < V; ++i) {
        inf_with[i] = -1;
        immunity[i].reset();
        variants[i].reset();
    }

    // Setup rand
    randVec.reserve(dna_len);
    for (int i = 0; i < dna_len; i++) {
        randVec.push_back(i);
    }
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    shuffle(randVec.begin(), randVec.end(),
            default_random_engine(seed));

    // Make initial variant
    for (int i = 0; i < init_bits; i++) {
        variants[0].set(randVec.at(i));
    }

    for (int i = 0; i < max_vars; i++) {
        vorigs[i] = -2;
        vprofs[i] = (vector<int>) {};
        vprofs[i].reserve(max_len);
        vtimes[i].first = -1;
        vtimes[i].second = -1;
        new_inf[i] = 0;
    }

    vcnt = 0;
    setC2(0); // All susceptible
    clr[p0] = 1; // Patient zero infected
    NI = 1;
    inf_with[0] = 0;
    day = 0;
    vtimes[0].first = 0;
    vprofs[0].push_back(1);
    vorigs[0] = -1;
    inf_neigh.reserve(V);
    immunity[0] = variants[0];

    while (NI > 0) { // While infected remain
        // Clear infected neighbours
        inf_neigh.clear();
        day++;
        if (day == max_len) {
            cout << "ERROR: Epidemic Length exceeds Maximum Epidemic Length" << endl;
        }
        for (int i = 0; i < V; i++) {
            inf_neigh.push_back((vector<int>) {});
        }

        // Find infected neighbours
        for (int inf = 0; inf < V; inf++) {
            if (clr[inf] == 1) { // Found infected individual
                for (int nei = 0; nei < nbr[inf].size(); nei++) {
                    inf_neigh.at(nbr[inf].memz(nei)).push_back(inf);
                }
            }
        }

        // Determine infections
        for (int node = 0; node < V; node++) {
            // If infected neighbours
            if (clr[node] == 0 && !inf_neigh.at(node).empty()) {
                strains.clear();
                for (int n: inf_neigh.at(node)) {
                    strains.emplace_back(inf_with[n], variants[inf_with[n]]);
                }

                strain_id = -1;
                varInfected(immunity[node], strains, alpha, strain_id);

                // Infected and update immunity
                if (strain_id != -1) { // Infected with inf_strain
                    clr[node] = 2; // Marked infected
                    new_inf[strain_id]++;
                    // New variant?
                    if (var_prob > 0 && drand48() < var_prob) { // Yes
                        vcnt++;
                        if (vcnt == max_vars) {
                            cout << "ERROR: Variant Count exceeds Maximum Variants" << endl;
                        }
                        create_new_variant(variants[strain_id],
                                           variants[vcnt],
                                           randVec, lB, uB);
                        new_inf[vcnt]++;
                        inf_with[node] = vcnt;
                        immunity[node] = immunity[node] | variants[vcnt];
                        vorigs[vcnt] = strain_id;
                        vtimes[vcnt].first = day;
                    } else { // No
                        inf_with[node] = strain_id;
                    }
                    immunity[node] = immunity[node] | variants[strain_id];
                }
            }
        }

        NI = 0;
        for (int node = 0; node < V; node++) {
            switch (clr[node]) {
                case 0:
                    break;
                case 1:
                    clr[node] = 0;
                    inf_with[node] = -1;
                    break;
                case 2:
                    clr[node] = 1;
                    NI++;
                    break;
            }
        }
        for (int s = 0; s <= vcnt; s++) {
            if (vtimes[s].second == -1) {
                if (new_inf[s] == 0) { // End of variant
                    vtimes[s].second = day;
                }
            }
            vprofs[s].push_back(new_inf[s]);
            new_inf[s] = 0;
        }
    }
}

//This is a modification of the SIR routine that adds profile - which
//should be a double array with as many positions as the number of
//vertices.  It returns the number of people infected in each time step.
void graph::SIRProfile(int p0, int &max, int &len, int &ttl, double alpha,
                       double *prof) {//Sir Method, with profile

    int NI;     //number of infected individuals
    int i, j, k;  //loop index variables
    int *nin;   //number of infected neighbors


    max = len = ttl = 0; //zero the reporting statistics
    if ((V == 0) || (p0 < 0) || (p0 >= V))
        return; //no one was infected...

    for (i = 0; i < V; i++)prof[i] = 0;  //zero the profile array

    nin = new int[V]; //create infected neioghbor counters
    setC2(0);      //set the population to infected
    clr[p0] = 1;      //infect patient zero
    NI = 1;           //initialize to one person currently infected
    len = 0;          //initialize length variable
    prof[len] = 1.0;  //record patient zero

    while (NI > 0) {//while there is still an epidemic
        //cout << "LEN=" << len << " NI=" << NI << endl;
        for (i = 0; i < V; i++)
            nin[i] = 0; //zero the number of infected neighbors buffer
        for (i = 0; i < V; i++)
            if (clr[i] == 1) {//found infected individual
                for (j = 0; j < nbr[i].size(); j++)
                    nin[nbr[i].memz(j)]++; //record exposure
            }
        //check for transmission
        for (i = 0; i < V; i++)
            if ((clr[i] == 0) && (nin[i] > 0))
                clr[i] = 3 * infected(nin[i], alpha);
        if (NI > max)max = NI; //check for updated maximum
        ttl += NI; //add the infected to the total
        NI = 0;  //zero the number infected counter
        for (i = 0; i < V; i++)
            switch (clr[i]) {//status update
                case 0: //susceptible, do nothing
                    break;
                case 1: //infected, move to removed
                    clr[i] = 2;
                    break;
                case 2: //removed, do nothing
                    break;
                case 3: //newly infected
                    clr[i] = 1;
                    NI++;
                    prof[len + 1] += 1.0;  //record the infection
                    break;
            }
        len++; //record the time step
    }
    delete[] nin;  //return storage for nin buffer
}

/* This routine duplicates SIR except that it assigns patient zero at random */
void
graph::SIRr(int &max, int &len, int &ttl, double alpha) {//Sir Method

    int NI;     //number of infected individuals
    int i, j, k;  //loop index variables
    int *nin;   //number of infected neighbors

    max = len = ttl = 0; //zero the reporting statistics

    nin = new int[V]; //create infected neioghbor counters
    setC2(0);    //set the population to infected
    clr[lrand48() % M] = 1;   //choose and infect patient zero
    NI = 1;        //initialize to one person currently infected
    len = 0;       //initialize length variable
    while (NI > 0) {//while there is still an epidemic
        //cout << "LEN=" << len << " NI=" << NI << endl;
        for (i = 0; i < V; i++)
            nin[i] = 0; //zero the number of infected neighbors buffer
        for (i = 0; i < V; i++)
            if (clr[i] == 1) {//found infected individual
                for (j = 0; j < nbr[i].size(); j++)
                    nin[nbr[i].memz(j)]++; //record exposure
            }
        //check for transmission
        for (i = 0; i < V; i++)
            if ((clr[i] == 0) && (nin[i] > 0))
                clr[i] = 3 * infected(nin[i], alpha);
        if (NI > max)max = NI; //check for updated maximum
        ttl += NI; //add the infected to the total
        NI = 0;  //zero the number infected counter
        for (i = 0; i < V; i++)
            switch (clr[i]) {//status update
                case 0: //susceptible, do nothing
                    break;
                case 1: //infected, move to removed
                    clr[i] = 2;
                    break;
                case 2: //removed, do nothing
                    break;
                case 3: //newly infected
                    clr[i] = 1;
                    NI++;
                    break;
            }
        len++; //record the time step
    }
    delete[] nin;  //return storage for nin buffer
}
/***********************************************************************/
/* The followinr SIR method superceeds the others, they are included
     for upward compatibility.
  */


//Alpha is the chance of transmission, advance is the I->R probability
//if q is -1 then it is assigned the first time the routine is called
//This routine returns the number of sick people
//The status of patients are stored in the quality
int
graph::SIRupdate(int &q, double alpha, double advance) {//one time step

    int i, j;        //loop indices
    int nsp;        //number of sick people
    int cp;         //current patient quality
    int nsn;        //number of size neighbors
    int deg;        //degree of a node
    int nq;         //neighbor's quality
    int scratch[V]; //new qualities


    if (q == -1) {//Allocate a quality
        q =
            RQquality(0);  //request a quality and allocate it to susceptable
        RecordQ(q, 0, 1);  //make the zeroth vertix patient zero
    }

    nsp = 0;  //zero the sick person counter
    for (i = 0; i < V; i++) {//loop over vertices
        cp = RetrieveQ(q, i);
        switch (cp) {//SIR
            case 0: //susceptible
                deg = degree(i);  //get the number of neighbors
                nsn = 0;          //zero the sick neighbor counter
                for (j = 0; j < deg; j++) {//loop over neighbors
                    nq = RetrieveQ(q,
                                   nbrmod(i,
                                          j));   //get the neighbor's status
                    if (nq == 1)
                        nsn++;                //count the sick neighbor
                }//now we know the number of sick neighbors
                if ((nsn > 0) && (infected(nsn, alpha))) {//new infected
                    scratch[i] = 1;  //note the infection
                    nsp++;         //record the infection
                } else scratch[i] = 0;  //stayed susceptible
                break;
            case 1: //infected
                if (drand48() < advance)scratch[i] = 2; //remove!
                else {
                    scratch[i] = 1;
                    nsp++;
                }         //preserve and record infection
                break;
            case 2: //removed
                scratch[i] = 2;     //stay removed
                break;
        }
    }
    for (i = 0; i < V; i++)
        RecordQ(q,
                i,
                scratch[i]);  //update the quality
    return (nsp);                //return the number of sick people


}

/* SEIR   S=0   E=1   I=2   R=3   */

//Alpha is the chance of transmission, creating an E, infec is the
//chance of E advancing to I, asym is the chance of E going direct to
//I -- that is an asymptomtic person, advance is the I->R probability
//if q is -1 then it is assigned the first time the routine is called
//This routine returns the number of sick people The status of
//patients are stored in the quality
int graph::SEIRupdate(int &q, double alpha,   //pairwise infection
                      double infec,   //E->I transition
                      double asym,    //E->R transition
                      double advance  //I->R transition
) {//one time step

    int i, j;        //loop indices
    int nsp;        //number of sick people
    int cp;         //current patient quality
    int nsn;        //number of size neighbors
    int deg;        //degree of a node
    int nq;         //neighbor's quality
    int scratch[V]; //new qualities
    double dart;    //used to resolve multiple probabilities


    if (q == -1) {//Allocate a quality
        q =
            RQquality(0);  //request a quality and allocate it to susceptable
        RecordQ(q, 0, 1);  //make the zeroth vertix patient zero
    }

    nsp = 0;  //zero the sick person counter
    for (i = 0; i < V; i++) {//loop over vertices
        cp = RetrieveQ(q, i);  //get the agent's status
        switch (cp) {//SIR
            case 0: //susceptible
                deg = degree(i);  //get the number of neighbors
                nsn = 0;          //zero the sick neighbor counter
                for (j = 0; j < deg; j++) {//loop over neighbors
                    nq = RetrieveQ(q,
                                   nbrmod(i,
                                          j));   //get the neighbor's status
                    if ((nq == 1) || (nq == 2))
                        nsn++;     //count the sick neighbors
                }//now we know the number of sick neighbors
                if ((nsn > 0) && (infected(nsn, alpha))) {//new infected
                    scratch[i] = 1;  //note the infection
                    nsp++;         //record the infection
                } else scratch[i] = 0;  //stayed susceptible
                break;
            case 1: //Exposed
                dart = drand48();
                if (dart < infec) {//exposed transitioned to infected
                    scratch[i] = 2;   //record infected status
                    nsp++;          //increment number of sick people
                } else if (dart
                    < infec + asym) {//transition directly to removed
                    scratch[i] = 3;   //record removed status
                } else {//stay exposed
                    scratch[i] = 1;   //record exposed status
                    nsp++;          //increment number of sick people
                }
                break;
            case 2: //infected
                if (drand48() < advance)scratch[i] = 3; //remove!
                else {
                    scratch[i] = 2;
                    nsp++;
                }         //preserve and record infection
                break;
            case 3: //removed
                scratch[i] = 3;     //stay removed
                break;
        }
    }
    for (i = 0; i < V; i++)
        RecordQ(q,
                i,
                scratch[i]);  //update the quality
    return (nsp);                //return the number of sick people

}

/* SEEpIR   S=0   E=1   Ep=2   I=3   R=4   */

//Alpha is the chance of transmission, creating an E, infec is the
//chance of E advancing to Ep, visible is the chance of the Ep-I
//transition, asym is the chance of Ep going to R directly, skipping
//I, and finally advance is the chance of the I->R transition
//---
//If q is -1 then it is assigned the first time the routine is
//called.  It is the quality that stores agent staus.  This routine
//returns the number of sick people The status of patients are
//stored in the quality
int graph::SEEpIRupdate(int &q, double alpha,   //pairwise infection
                        double infec,   //E->Ep transition
                        double visible, //Ep->I transition
                        double asym,    //Ep->R transition
                        double advance  //I->R transition
) {//one time step


    int i, j;        //loop indices
    int nsp;        //number of sick people
    int cp;         //current patient quality
    int nsn;        //number of size neighbors
    int deg;        //degree of a node
    int nq;         //neighbor's quality
    int scratch[V]; //new qualities
    double dart;    //used to resolve multiple probabilities


    if (q == -1) {//Allocate a quality
        q =
            RQquality(0);  //request a quality and allocate it to susceptable
        RecordQ(q, 0, 1);  //make the zeroth vertix patient zero
    }

    nsp = 0;  //zero the sick person counter
    for (i = 0; i < V; i++) {//loop over vertices
        cp = RetrieveQ(q, i); //get the agent's status
        switch (cp) {//SIR
            case 0: //susceptible
                deg = degree(i);  //get the number of neighbors
                nsn = 0;          //zero the sick neighbor counter
                for (j = 0; j < deg; j++) {//loop over neighbors
                    nq = RetrieveQ(q,
                                   nbrmod(i,
                                          j));   //get the neighbor's status
                    if ((nq == 2) || (nq == 3))
                        nsn++;     //count the sick neighbors
                }//now we know the number of sick neighbors
                if ((nsn > 0) && (infected(nsn, alpha))) {//new infected
                    scratch[i] = 1;  //note the infection
                    nsp++;         //record the infection
                } else scratch[i] = 0;  //stayed susceptible
                break;
            case 1: //Exposed
                dart = drand48();
                if (dart < infec) {//exposed transitioned to infected
                    scratch[i] =
                        2;   //record infected status, moves to Ep
                    nsp++;          //increment number of sick people
                } else {//stay exposed
                    scratch[i] = 1;   //record exposed status
                    nsp++;          //increment number of sick people
                }
                break;
            case 2: //Exposed prime
                dart = drand48();
                if (dart
                    < visible) {//exposed prime transitioned to infected
                    scratch[i] =
                        3;   //record infected status, moves to Ep
                    nsp++;          //increment number of sick people
                } else {//stay exposed prime
                    scratch[i] = 2;   //record exposed status
                    nsp++;          //increment number of sick people
                }
                break;
            case 3: //Infected
                if (drand48() < advance)scratch[i] = 4; //remove!
                else {
                    scratch[i] = 3;
                    nsp++;
                }         //preserve and record infection
                break;
            case 4: //removed
                scratch[i] = 4;     //stay removed
                break;
        }
    }
    for (i = 0; i < V; i++)
        RecordQ(q,
                i,
                scratch[i]);  //update the quality
    return (nsp);                              //return the number of sick people

}

//This version adds the longer term asymptomatic class A
//S=0, E=1, Ep=2, I=3, A=4, R=5
int graph::SEEIARupdate(int &q, double alpha,   //pairwise infection
                        double infec,   //E->Ep transition
                        double infecP,  //Ep->I transition
                        double cksee,   //Ep->A transition
                        double gamma,   //A,I -> R transition
                        double epsilon  //R->S transition
) {//one time step


    int i, j;        //loop indices
    int nsp;        //number of sick people
    int cp;         //current patient quality
    int nsn;        //number of size neighbors
    int deg;        //degree of a node
    int nq;         //neighbor's quality
    int scratch[V]; //new qualities
    double dart;    //used to resolve multiple probabilities


    if (q == -1) {//Allocate a quality
        q =
            RQquality(0);  //request a quality and allocate it to susceptable
        RecordQ(q, 0, 1);  //make the zeroth vertix patient zero
    }

    nsp = 0;  //zero the sick person counter
    for (i = 0; i < V; i++) {//loop over vertices
        cp = RetrieveQ(q, i); //get the agent's status
        switch (cp) {//SIR
            case 0: //susceptible
                deg = degree(i);  //get the number of neighbors
                nsn = 0;          //zero the sick neighbor counter
                for (j = 0; j < deg; j++) {//loop over neighbors
                    nq = RetrieveQ(q,
                                   nbrmod(i,
                                          j));   //get the neighbor's status
                    if ((nq != 0) && (nq != 5))
                        nsn++;     //count the sick neighbors
                }//now we know the number of sick neighbors
                if ((nsn > 0) && (infected(nsn, alpha))) {//new infected
                    scratch[i] = 1;  //note the infection
                    nsp++;         //record the infection
                } else scratch[i] = 0;  //stayed susceptible
                break;
            case 1: //Exposed
                dart = drand48();
                if (dart < infec) {//exposed transitioned to infected
                    scratch[i] =
                        2;   //record infected status, moves to Ep
                    nsp++;          //increment number of sick people
                } else {//stay exposed
                    scratch[i] = 1;   //record exposed status
                    nsp++;          //increment number of sick people
                }
                break;
            case 2: //Exposed prime
                dart = drand48();
                if (dart
                    < infecP) {//exposed prime transitioned to infected
                    scratch[i] =
                        3;   //record infected status, moves to Ep
                    nsp++;          //increment number of sick people
                } else if (dart
                    < infecP
                        + cksee) {//Exposed prime transition to long term A
                    scratch[i] = 4;
                    nsp++;
                } else {//stay exposed prime
                    scratch[i] = 2;   //record exposed status
                    nsp++;          //increment number of sick people
                }
                break;
            case 3: //Infected
                dart = drand48();  //throw the dart
                if (dart < gamma)scratch[i] = 5;   //remove!
                else {
                    scratch[i] = 3;
                    nsp++;
                }    //preserve and record infection
                break;
            case 4: //Long term asymptomatic
                dart = drand48();  //throw the dart
                if (dart < gamma)scratch[i] = 5;   //remove!
                else {
                    scratch[i] = 4;
                    nsp++;
                }    //preserve and record infection
                break;
            case 5: //removed
                dart = drand48();  //throw the dart
                if (dart < epsilon)scratch[i] = 0; //become susceptible
                else scratch[i] = 5;            //stay removed
                break;

        }
    }
    for (i = 0; i < V; i++)
        RecordQ(q,
                i,
                scratch[i]);  //update the quality
    return (nsp);                              //return the number of sick people
}

int graph::attack(double *pr) {//probabalistic attack method

    if ((M == 0) || (V == 0))
        return (0); //treat empty graphs as totally vulnerable

    int cnt;     //counter for number of knocked out vertices
    int *q;      //distance buffer
    int rv;      //return value
    int i, j;     //index variables
    int kill;    //vertex to kill
    double ttl;  //total of probability vector
    int fl;      //stop flag

    if (clr == 0)
        clr = new int[M];  //if the color buffer doesn't exist, create
    for (i = 0; i < V; i++)
        clr[i] = 0;   //color zero is ``vertex alive''
    q = new int[V];  //create distance array
    ttl = 0;
    for (i = 0; i < V; i++)
        ttl += pr[i];  //create total probability mass
    rv = 0; //initialize the return value
    fl = 0; //reset disconnection flag
    do {//knock out vertices
        i = 0;
        while ((i < V) && (clr[i] != 0))i++; //find first living vertex
        if (i >= V)break; //D'oh  Ate the whole graph
        dfrom0(i, q); //get zero-color distances from live vertex
        for (j = 0; j < V; j++)
            if ((clr[j] == 0) && (q[j] == -1))
                fl = 1; //disconnected!
        if (fl == 1)break;  //break out of the do loop
        rv++; //survived a connectedness check
        do { kill = rselect(pr, ttl, V); }
        while (clr[kill] == 1); //find living vertex
        clr[kill] = 1; //knock out the vertex
    } while (1); //infinite look with breakouts

    delete[] q;  //give back the storage
    return (rv); //return the number of cycles

}

/************ANT METHODS GO HERE*************************************/
//Start an ant at each vertex and run it for reps steps and return
//the number of ands entering each vertex
void graph::Ant(int reps, int *cnt) {

    int i;      //loop index variable
    int *psn;   //ant positions


}

//color methods
void graph::setC2(int vl) {//set colors to vl

    if (clr == 0)
        clr = new int[M];  //if the color buffer doesn't exist, create
    for (int i = 0; i < V; i++)clr[i] = vl; //set the value

}

void
graph::GDC(int *ord) {//run the greedy coloring algorithm with order ord

    int i, j, k, m;
    static int F[50];

    setC2(-1);  //set current colors to -1
    for (i = 0; i < V; i++) {//loop over vertices
        for (j = 0; j < 50; j++)F[j] = 0; //zero use buffer
        k = nbr[ord[i]].size(); //get degree
        for (j = 0; j < k; j++) {//loop over neighbors
            m = clr[nbr[ord[i]].memz(j)];
            if (m >= 0)F[m] = 1;
        }
        for (j = 0; (j < 50) && (F[j] == 1); j++);
        clr[i] = j;
    }
}

int graph::AUC(double *gn,
               int tg) {//run Austrailian coloring with target tg

    int *used, *aval;  //used and available colors
    int i, j, k;        //loop idex variables
    int v;            //vertex buffer

    used = new int[tg];  //allocate used color array
    aval = new int[tg];  //allocate available color array
    if (clr == 0)clr = new int[M]; //no color? allocate it
    for (i = 0; i < V; i++)clr[i] = -1; //mark as uncolored
    for (i = 0; i < V; i++) {//loop over vertices
        for (j = 0; j < tg; j++)used[j] = 0; //mark all colors unused
        for (j = 0; j < nbr[i].size(); j++) {//loop over neighbors
            v = nbr[i].memz(j); //get the current neighbor
            if (clr[v] >= 0)used[clr[v]] = 1; //make the color used
        }
        k = 0; //zero available color pointer
        for (j = 0; j < tg; j++)
            if (used[j] == 0)
                aval[k++] = j; //compile available colors
        //cout << "K=" << k << endl;
        //for(j=0;j<V;j++)cout << clr[j] << " ";cout << endl;
        if (k == 0)return (i); //stuck!
        clr[i] = aval[((int) (k * gn[i]))];
        //cout << i << " " << clr[i] << endl;
    }
    delete[] aval;    //return the storage
    delete[] used;    //return the storage
    return (V);         //managed to color the graph
}


//genetics
int
graph::FPN(int v, double *ft) {//fitness perprotional neighbor selection

    if ((v >= V) || (v < 0))
        return (0);  //return zero for stupid request
    return (nbr[v].FPS(ft));      //call the set fitness proportional selection
}

//I-O
void graph::write(ostream &aus) {//write the graph

    int i;

    aus << M << "\t" << V << "\t" << E << endl;
    for (i = 0; i < V; i++) {
        nbr[i].writememb(aus);
    }
    //Need to add code to write qualities and weights
}

void graph::read(istream &inp) {//read the graph

    char buf[1000];
    int k;

    if (M != 0)destroy(); //clear the graph if it is not already clear
    inp.getline(buf, 999);
    M = atoi(buf);
    k = 0;
    while (buf[k] != ' ')k++;
    while (buf[k] == ' ')k++;
    V = atoi(buf + k);
    while (buf[k] != ' ')k++;
    while (buf[k] == ' ')k++;
    E = atoi(buf + k);
    nbr = new set[M];
    for (int i = 0; i < V; i++) {
        nbr[i].setempty();    //safety first
        nbr[i].readmemb(inp); //read members
    }
    //Need to add code to read qualities and weights
}
void
graph::readadj(istream &inp, int numV) {//read a python adjacency graph

    char buf[1000];
    int source, target;
    char *tokenPtr;

    if (M != 0)destroy(); //clear the graph if it is not already clear
    //initialize graph
    empty(numV);
    //discard header
    do inp.getline(buf, 999); while (buf[0] == '#');
    //read graph
    for (int i = 0; i < V; i++) {
        tokenPtr = strtok(buf, " ");
        if (tokenPtr != NULL) {
            source = atoi(tokenPtr);
            tokenPtr = strtok(NULL, " ");
            while (tokenPtr != NULL) {
                target = atoi(tokenPtr);
                add(source, target);
                tokenPtr = strtok(NULL, " ");
            }
            inp.getline(buf, 999);
        }
    }
}

void graph::writeC(ostream &aus) {

    if (clr == 0)return;

    int i;  //loop index

    aus << clr[0];  //first color
    for (i = 1; i < V; i++)aus << " " << clr[i]; //rest of colors
    aus << endl; //endline

}

void graph::dotout(ostream &aus) {//write the dot format

    int v, e;  //vertex and edge loop indices

    aus << "graph {" << endl;
    for (v = 0; v < V; v++) {//loop over the vertices
        for (e = v + 1; e < V; e++) {//loop over possible neighbors
            if (edgeP(v, e)) {//if this is an edge
                aus << v << " -- " << e << ";" << endl;
            }
        }
    }
    aus << "}" << endl;

}

void
graph::adjout(ostream &aus,
              char sep) {//write the adjacency list with seperator

    if (V == 0)return;
    for (int i = 0; i < V; i++) {//loop over vertices
        nbr[i].writememb(aus, sep);
    }
}

void
graph::writeadj(ostream &aus) {//write the graph as a python adjacency graph

    time_t ttime = time(0);

    //header
    aus << "# Output from setu.cpp" << endl;
    aus << "# " << ctime(&ttime);
    aus << "#" << endl;
    for (int i = 0; i < V; i++) {
        aus << i;
        for (int j = 0; j < degree(i); j++) {
            if (nbr[i].memz(j) > i) aus << " " << nbr[i].memz(j);
        }
        aus << endl;
    }

}


