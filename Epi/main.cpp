#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <string>
#include <cstring>
#include <filesystem>

using namespace std;

#include "setu.h"
#include "stat.h"
#include "../BitSprayer/bitspray.h"

/*************************algorithm controls******************************/
#define PL 16
#define NSE 50
#define alpha 0.5
#define mepl 3              //  Minimum epidemic length
#define rse 5               //  Re-try short epidemics
//#define ftl 50            //  Final test length
#define verbose true
#define runs 30
#define mevs 40000
#define RIs 100
#define RE ((long)mevs/RIs)
#define NmC (long)9
#define EDGB 2              //  Minimum degree for swap
#define popsize 100
#define verts 160
#define GL 256
#define RNS 91207819
#define MAXL (long)pow(verts, 3)  //  Size of integer
#define MNM 4
#define tsize 5

#define omega 0.5
#define edgeAdd 1
#define triProb 0.5

#define states 12
#define Qz 100000
int Q[Qz];

#define ent_thres 1.0
bool dead[popsize];
double max_fit = verts;

/**************************Variable dictionary************************/
int pop[popsize][GL];           //  Population of command strings
bitspray *bPop[popsize];
double fit[popsize];            //  Fitness array
int dx[popsize];                //  Sorting index
double PD[PL + 1];              //  Profile dictionary
double CmD[NmC];                //  Command densities
bool profiles;
bool ringG;
bool variants;
graph iG;

/****************************Procedures********************************/
void initalg(const char *pLoc);           //initialize the algorithm
int validloci();                          //generate an acceptable large integer
int validloci(int &psn);
void express(graph &G, const int *cmd);   //express a command string
void initpop();                           //initialize population
void matingevent();                       //run a mating event
void report(ostream &aus);                //make a statistical report
void reportbest(ostream &aus, ostream &difc);
void createReadMe(ostream &aus);
void cmdLineIntro(ostream &aus);
void cmdLineRun(int run, ostream &aus);
double fitness(int *cmd);
double fitness(int idx, bitspray &A);
void developQ(bitspray &A);
bool getnum(int &val, int bits, int &psn);
void getGene(int idx, double *probs);

/****************************Main Routine*******************************/
int main(int argc, char *argv[]) {
  /**
   * Output Root, Profile Location, Profile Number
   */

  profiles = false;   //  F -> Epidemic Length; T -> Profile Matching
  ringG = false;  //  Ring graph or power law cluster
  variants = false;   //  With or without variants

  fstream stat, best, dchar, readme, iGOut;   //statistics, best structures
  char fn[60];          //file name construction buffer
  char *outLoc = new char[45];
  char *outRoot = argv[1];
  char *pLoc = argv[2];
  int pNum = stoi(argv[3]);

  // Create directory for output
  if (profiles) {
    sprintf(outLoc, "%sOutput - Profile%d w %dS, %02dP, %dM/",
            outRoot, pNum, states, popsize, MNM);
  } else {
    sprintf(outLoc, "%sOutput - Length w %dS, %02dP, %dM/",
            outRoot, states, popsize, MNM);
  }
  std::filesystem::create_directory(outLoc);

  initalg(pLoc);

  sprintf(fn, "%sbest.lint", outLoc);
  best.open(fn, ios::out);
  sprintf(fn, "%sdifc.dat", outLoc);
  dchar.open(fn, ios::out);
  sprintf(fn, "%sreadme.dat", outLoc);
  readme.open(fn, ios::out);
    createReadMe(readme);
    readme.close();
    if (!ringG) {
        sprintf(fn, "%sinitGraph.dat", outLoc);
        iGOut.open(fn, ios::out);
        iG.write(iGOut);
        iGOut.close();
    }
    if (verbose) {
        cmdLineIntro(cout);
    }
    if (!verbose) {
        cout << "Started" << endl;
    }
    for (int run = 0; run < runs; run++) {
        sprintf(fn, "%srun%02d.dat", outLoc, run); // File name
        stat.open(fn, ios::out);
        if (verbose) cmdLineRun(run, cout);
        initpop();
        report(stat); //report the statistics of the initial population
        for (int mev = 0; mev < mevs; mev++) {//do mating events
            matingevent();  //run a mating event
            if ((mev + 1) % RE == 0) {//Time for a report?
                if (verbose) {
                    cout << left << setw(5) << run;
                    cout << left << setw(4) << (mev + 1) / RE;
                }
                report(stat); //report statistics
            }
        }
        stat.close();
        reportbest(best, dchar);
        cout << "Done run " << run << " of " << runs - 1 << endl;
    }
    best.close();
    dchar.close();
    return (0);  //keep the system happy
}

void createReadMe(ostream &aus) {
  aus << "This file contains the info about the files in this folder."
      << endl;
  aus << "Graph Evolution Tool." << endl;

  aus << "Initial Graph: ";
  if (ringG) {
    aus << "Ring with m = 2" << endl;
  } else {
    aus << "Power Law Cluster Graph" << endl;
  }

  aus << "Fitness Function: ";
  if (profiles) {
    aus << "Profile Matching" << endl;
    aus << "Profile: ";
    for (double i: PD) {
      aus << i << " ";
    }
    aus << endl;
  } else {
    aus << "Epidemic Length" << endl;
  }

  aus << "Epidemic Model Used: ";
  if (variants) {
    aus << "SIR with Variants" << endl;
  } else {
    aus << "SIR" << endl;
  }

  aus << "Representation: Bitsprayers on Adjacency Matrix" << endl;
  aus << endl;
  aus << "The parameter settings are as follows: " << endl;
  aus << "Number of sample epidemics: " << NSE << endl;
  aus << "Alpha: " << alpha << endl;
  aus << "Minimum epidemic length: " << mepl << endl;
  aus << "Re-tries for short epidemics: " << rse << endl;
  aus << "Runs: " << runs << endl;
  aus << "Mating events: " << mevs << endl;
  aus << "Minimum degree for swap: " << EDGB << endl;
  aus << "Population size: " << popsize << endl;
  aus << "Number of vertices: " << verts << endl;
  aus << "Maximum number of mutations: " << MNM << endl;
  aus << "Tournament size: " << tsize << endl;
  aus << "Number of States: " << states << endl;
//    aus << "Entropy Threshold for Necrotic Filter: " << ent_thres << endl;
  aus << "Decay strength for diffusion characters: " << omega << endl;
  if (!ringG) {
    aus << "Omega value: " << omega << endl;
    aus << "Edges added at each step: " << edgeAdd << endl;
    aus << "Triangle creation probability: " << triProb << endl;
  }
  aus << endl;
  aus << "The file descriptions are as follows: " << endl;
  aus << "best.lint -> the best fitness and it's associated data for each "
         "run" << endl;
  aus << "difc.dat -> the diffusion characters of the best graph for each "
         "run" << endl;
  aus << "run##.dat -> population statistics for each run" << endl;
  if (!ringG) {
    aus << "initGraph.dat -> the initial power law clustering graph"
        << endl;
  }
}

void cmdLineIntro(ostream &aus) {
  aus << "Graph Evolution Tool." << endl;
  aus << "Initial Graph: ";
  if (ringG) {
    aus << "Ring with m = 2" << endl;
  } else {
    aus << "Power Law Cluster Graph" << endl;
  }

  aus << "Fitness Function: ";
  if (profiles) {
    aus << "Profile Matching" << endl;
    aus << "Profile: ";
    for (double i: PD) {
      aus << i << " ";
    }
    aus << endl;
  } else {
    aus << "Epidemic Length" << endl;
  }

  aus << "Epidemic Model Used: ";
  if (variants) {
    aus << "SIR with Variants" << endl;
  } else {
    aus << "SIR" << endl;
  }
  aus << "Representation: Bitsprayers on Adjacency Matrix" << endl;
  aus << "Check readme.dat for more information about parameters/output.";
  aus << endl;
}

void cmdLineRun(int run, ostream &aus) {
    aus << endl << "Beginning Run " << run << " of " << runs - 1 << endl;
    aus << left << setw(5) << "Run";
    aus << left << setw(4) << "RI";
    aus << left << setw(10) << "Mean";
    aus << left << setw(12) << "95% CI";
    aus << left << setw(10) << "SD";
    aus << left << setw(8) << "Best";
    aus << endl;
    aus << left << setw(5) << run;
    aus << left << setw(4) << "0";
}

void initalg(const char *pLoc) {//initialize the algorithm
    fstream inp;    //input file
    char buf[20];   //input buffer

    srand48(RNS);                 //read the random number seed
  if (!ringG) {
    iG.create(verts);
    iG.PCG(verts, edgeAdd, triProb);
  }
  if (profiles) {
    inp.open(pLoc, ios::in);      //open input file
    for (int i = 0; i < PL; i++) {
      PD[i] = 0;  //pre-fill missing values
    }
    PD[0] = 1;  //put in patient zero
    for (int i = 0; i < PL; i++) {      //loop over input values
      inp.getline(buf, 19);  //read in the number
      PD[i + 1] = strtod(buf, nullptr);    //translate the number
    }
    inp.close();
  }
  for (auto &i: bPop) {
    i = new bitspray(states);
  }
}

//This routine generates valid loci for the expression routine
int validloci() {//generate an acceptable large integer
    int cmd;       //command type generated
    double dart;   //Random command

    dart = drand48() - CmD[0];//throw the dart
    cmd = 0;                  //initialize the position on the dartboard
    while ((dart > 0) && (cmd < NmC - 1)) {
        dart -= CmD[++cmd]; //walk the board
    }
    cmd += (int) (NmC * (lrand48() % MAXL)); //add in the large integer part

    return (cmd);  //return the generated command
}

int validloci(int &psn) {
    int cmd;
    int lint;
    int size = ceil(log2(NmC));

    do {
        if (!getnum(cmd, size, psn)) return -1;
    } while (cmd > NmC - 1);
    size = ceil(log2(MAXL));
    if (!getnum(lint, size, psn)) return -1;
    cmd += (int) (NmC * (lint % MAXL));

    return cmd;
}

//This is expression of a large integer represenation
void express(graph &G, const int *cmd) {//express a command string
    int a, b, c;  //decoded values
    int cdv;    //command value
    int block;  //integer carving block

    if (ringG) {
        G.RNGnm(verts, 2);  //  Initial graph
    } else {
        G.copy(iG);
    }

    for (int i = 0; i < GL; i++) {//loop over the commands (genetic loci)
        block = cmd[i];  //get integer
        cdv = (int) (block % NmC); //slice of the command
        block /= NmC;    //clear command information from block
        switch (cdv) {  //What command is it?
            case 0://Toggle
                a = block % verts;
                b = (block / verts) % verts; //get vertex numbers
                G.toggle(a, b);    //toggle edge {a,b}
                break;
            case 1://Hop
                //get vertex numbers
                a = block % verts;
                b = (block / verts) % verts;
                c = (block / verts / verts) % verts;
                G.hop(a, b, c);
                break;
            case 2://Add
                a = block % verts;
                b = (block / verts) % verts; //get vertex numbers
                G.add(a, b);       //add edge {a,b}
                break;
            case 3://Delete
                a = block % verts;
                b = (block / verts) % verts; //get vertex numbers
                G.del(a, b);       //delete edge a,b
                break;
            case 4://Swap
                a = block % (verts * 10);
                b = (block / (verts * 10)) % (verts * 10); //get vertex numbers
                G.edgeswap(a, b, EDGB);
                break;
            case 5://Local Toggle
                //get vertex numbers
                a = block % verts;
                b = (block / verts) % verts;
                c = (block / verts / verts) % verts;
                G.loggle(a, b, c);
                break;
            case 6: // Local Add
                //get vertex numbers
                a = block % verts;
                b = (block / verts) % verts;
                c = (block / verts / verts) % verts;
                G.ladd(a, b, c);
                break;
            case 7: // Local Del
                //get vertex numbers
                a = block % verts;
                b = (block / verts) % verts;
                c = (block / verts / verts) % verts;
                G.ldel(a, b, c);
                break;
            case 8: // Null - Do nothing
                break;
            default:
                cout << "ERROR: default of switch";
                break;
        }
    }
}

int bitsToInt(const int bs[], int cnt) {
    int val = bs[0];
    for (int e = 1; e < cnt; e++) {
        val += (int) pow(bs[e], e);
    }
    return val;
}

bool necroticFilter() {
    double En = 0.0;
    int binSize = 10;
    int numBins = (int) pow(2, binSize);
    double counts[numBins];
    int vals = Qz / binSize;
    for (int i = 0; i < numBins; i++) {
        counts[i] = 0.0;
    }
    for (int i = 0; i < Qz; i += binSize) {
        counts[bitsToInt(Q + i, binSize)]++;
    }
    for (int i = 0; i < numBins; i++) {
        counts[i] /= vals;
    }
    for (int i = 0; i < numBins; i++) {
        if (counts[i] > 0.0) {
            En += -counts[i] * log(counts[i]);
        }
    }
    En /= log(2);
    return En < ent_thres;
}

double fitness(int *cmd) {//compute the epidemic length fitness
    graph G(verts);      //scratch graph
    int max, len, ttl;   //maximum, length, and total removed
    int cnt;             //counter for tries
    double prof[verts];  //profile variable
    double trials[NSE];  //stores squared error for each trial
    int en;              //epidemic number
    double delta;        //difference between profile and trial
    double accu = 0.0;         //accumulator
    express(G, cmd);               //create the graph

    if (mode == 0) {
        for (en = 0; en < NSE; en++) {
            cnt = 0;
            do {
                G.SIR(0, max, len, ttl, alpha);
                cnt++;
            } while (len < mepl && cnt < rse);
            trials[en] = len;
        }
        for (double trial: trials) {//loop over trials
            accu += trial;
        }
        accu = accu / NSE;
    } else {
        for (en = 0; en < NSE; en++) {//loop over epidemics
            cnt = 0;
            do {
                G.SIRProfile(0, max, len, ttl, alpha, prof);
                cnt++;
            } while (len < mepl && cnt < rse);
            trials[en] = 0;  //zero the current squared error
            if (len < PL + 1) {
                len = PL + 1;  //find length of epi/prof (longer)
            }
            for (int i = 0; i < len; i++) {//loop over time periods
                delta = prof[i] - PD[i];
                trials[en] += delta * delta;
            }
            trials[en] = sqrt(trials[en] / len); //convert to RMS error
        }

        for (double trial: trials) {
            accu += trial;
        }
        accu = accu / NSE;
    }
    return accu;  //return the fitness value
}

void developQ(bitspray &A) {//unpack the queue
    int h, t;  //head and tail of queue

    for (t = 0; t < Qz; t++) {
        Q[t] = 0;        //clear the queue
    }
    A.reset(Q, h, t);     //reset the self driving automata
    while (t < Qz - 2) {
        A.next(Q, h, t, Qz);  //run the automata
    }
}

bool getnum(int &val, int bits, int &psn) {//get a number from the queue
    if (psn + bits >= Qz - 1) {
//    cout << "ERROR: Q EMPTY" << endl;
        return false; //safety first
    }
    val = 0;  //zero the value
    for (int i = 0; i < bits; i++)
        val = val * 2 + Q[psn++];  //assemble the integer
    return true;
}

double fitness(int idx, bitspray &A) {
  int psn = 0;

  developQ(A);
//    if (necroticFilter()) {
//        dead[idx] = true;
//        return max_fit;
//    }
//    for (int &i: pop[idx]) {
//        i = validloci(psn);
//        if (i == -1) {
//            dead[idx] = true;
//            return max_fit;
//        }
//    }
  double fi = fitness(pop[idx]);
  return fi;
}

void getGene(int idx, double *probs) {
    fitness(idx, *bPop[idx]);
    for (int i = 0; i < NmC; i++) {
        probs[i] = 0.0;
    }
    for (int &c: pop[idx]) {
        probs[c % NmC]++;
    }
    for (int i = 0; i < NmC; i++) {
        probs[i] = (double) probs[i] / GL;
    }
}

void initpop() {//initialize population
  for (int i = 0; i < popsize; i++) {
    dead[i] = false;
    bPop[i]->randomize();
    fit[i] = fitness(i, *bPop[i]);
    dx[i] = i;
  }
}

void matingevent() {//run a mating event
    int rp, sw;   //loop index, random position, swap variable
    int cp1, cp2;   //crossover points

    //perform tournament selection, highest fitness first
    if (mode == 0) {
        tselect(fit, dx, tsize, popsize);
    } else {
        Tselect(fit, dx, tsize, popsize);
    }

    if (mode < 2) {
        //selection and crossover
        cp1 = (int) lrand48() % GL;
        cp2 = (int) lrand48() % GL;
        if (cp1 > cp2) {
            sw = cp1;
            cp1 = cp2;
            cp2 = sw;
        }
        for (int i = 0; i < cp1; i++) {
            pop[dx[0]][i] = pop[dx[tsize - 2]][i];
            pop[dx[1]][i] = pop[dx[tsize - 1]][i];
        }
        for (int i = cp1; i < cp2; i++) {
            pop[dx[0]][i] = pop[dx[tsize - 1]][i];
            pop[dx[1]][i] = pop[dx[tsize - 2]][i];
        }
        for (int i = cp2; i < GL; i++) {
            pop[dx[0]][i] = pop[dx[tsize - 2]][i];
            pop[dx[1]][i] = pop[dx[tsize - 1]][i];
        }

        //mutation
        rp = (int) lrand48() % MNM + 1;
        for (int i = 0; i < rp; i++) {
            pop[dx[0]][lrand48() % GL] = validloci();
        }
        rp = (int) lrand48() % MNM + 1;
        for (int i = 0; i < rp; i++) {
            pop[dx[1]][lrand48() % GL] = validloci();
        }

        //update fitness
        fit[dx[0]] = fitness(pop[dx[0]]);
        fit[dx[1]] = fitness(pop[dx[1]]);

        // Skeptical tournament selection
        if (mode == 0) {
            fit[dx[tsize - 1]] = fitness(pop[dx[tsize - 1]]);
            fit[dx[tsize - 2]] = fitness(pop[dx[tsize - 2]]);
        }
    } else {
        bPop[dx[0]]->copy(*bPop[dx[tsize - 2]]);
        bPop[dx[1]]->copy(*bPop[dx[tsize - 1]]);
        bPop[dx[0]]->tpc(*bPop[dx[1]]);
        rp = (int) lrand48() % MNM + 1;
        bPop[dx[0]]->mutate(rp);
        rp = (int) lrand48() % MNM + 1;
        bPop[dx[1]]->mutate(rp);
        // reset dead SDAs
        dead[dx[0]] = false;
        dead[dx[1]] = false;
        fit[dx[0]] = fitness(dx[0], *bPop[dx[0]]);
        fit[dx[1]] = fitness(dx[1], *bPop[dx[1]]);
    }
}

void report(ostream &aus) {//make a statistical report
    dset D;
    int deaths = 0;
    if (mode == 2) {
        for (bool i: dead) {
            if (i) {
                deaths++;
            }
        }
        double good_fit[popsize - deaths];
        for (double &i: good_fit) {
            i = 0.0;
        }
        int cnt = 0;
        for (int i = 0; i < popsize; i++) {
            if (!dead[i]) {
                good_fit[cnt++] = fit[i];
            }
        }
        D.add(good_fit, popsize - deaths);
        if (D.Rmax() != max_fit) {
            for (int i = 0; i < popsize; i++) {
                if (dead[i]) { // update max
                    fit[i] = D.Rmax();
                }
            }
            max_fit = D.Rmax();
        }
    } else {
        D.add(fit, popsize);  //load fitness
    }

    //print report
    if (mode == 0) {
        aus << left << setw(10) << D.Rmu();
        aus << left << setw(12) << D.RCI95();
        aus << left << setw(10) << D.Rsg();
        aus << left << setw(8) << D.Rmax() << endl;
        if (verbose) {
            cout << left << setw(10) << D.Rmu();
            cout << left << setw(12) << D.RCI95();
            cout << left << setw(10) << D.Rsg();
            cout << left << setw(8) << D.Rmax() << endl;
        }
    } else {
        aus << left << setw(10) << D.Rmu();
        aus << left << setw(12) << D.RCI95();
        aus << left << setw(10) << D.Rsg();
        aus << left << setw(8) << D.Rmin() << "\t";
        aus << "Dead: " << deaths << endl;
        if (verbose) {
            cout << left << setw(10) << D.Rmu();
            cout << left << setw(12) << D.RCI95();
            cout << left << setw(10) << D.Rsg();
            cout << left << setw(12) << D.Rmin() << "\t";
            cout << "Dead: " << deaths << endl;
        }
    }
}

void reportbest(ostream &aus, ostream &difc) {//report the best graph
    int b;       //loop indices and best pointer
    graph G(verts);  //scratch graph
    double En;
    static double M[verts][verts];
    static double Ent[verts];

    b = 0;
    if (mode == 0) {
        for (int i = 1; i < popsize; i++) {
            if (fit[i] > fit[b]) {
                b = i; //find best fitness
            }
        }
    } else {
        for (int i = 1; i < popsize; i++) {
            if (fit[i] < fit[b]) {
                b = i; //find best fitness
            }
        }
    }

    //output the fitness and the associated data
    aus << fit[b] << " -fitness" << endl;
    if (mode == 2) {
        double probs[NmC];
        getGene(b, probs);
        // Write the SDA
        aus << "Self-Driving Automata" << endl;
        bPop[b]->print(aus);
        // Command densities
        aus << "Resultant Command Densities" << endl;
        aus << probs[0];
        for (int i = 1; i < NmC; i++) {
            aus << " " << probs[i];
        }
        aus << endl;
    }
    // Write the gene
    aus << "Gene" << endl;
    aus << pop[b][0];
    for (int i = 1; i < GL; i++) {
        aus << " " << pop[b][i];
    }
    aus << endl;

    G.empty(verts);
    express(G, pop[b]);
    aus << "Graph" << endl;
    G.write(aus);
    aus << endl;
    for (int i = 0; i < G.size(); i++) {
        G.DiffChar(i, omega, M[i]);
    }

    for (int i = 0; i < G.size(); i++) {//loop over vertices
        En = 0.0;  //prepare En for normalizing
        for (int j = 0; j < G.size(); j++) {
            En += M[i][j];  //total the amount of gas at vertex i
        }
        for (int j = 0; j < G.size(); j++) {
            M[i][j] /= En;  //normalize so sum(gas)=1
        }
        En = 0.0;  //zero the entropy accumulator
        for (int j = 0; j < G.size();
             j++) {//build up the individual entropy terms
            if (M[i][j] > 0) {
                En += -M[i][j] * log(M[i][j]);  //this is entropy base E
            }
        }
        Ent[i] = En / log(2);  //convert entropy to Log base 2
    }

    //Now sort the entropy vector
    bool more = false;  //no swaps
    do {//swap out-of-order entries
        more = false;
        for (int i = 0; i < G.size() - 1; i++) {
            if (Ent[i] < Ent[i + 1]) {
                En = Ent[i];
                Ent[i] = Ent[i + 1];
                Ent[i + 1] = En;  //swap
                more = true;  //set the flag that a swap happened
            }
        }
    } while (more); //until in order

    difc << Ent[0];  //output first entropy value
    for (int i = 1; i < G.size(); i++) {
        difc << " " << Ent[i];  //output remaining values
    }
    difc << endl;  //end line
}

