/*
 * Code file for implementing bit sprayers, a.k.a. a self-driving automata class
 *
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "bitspray.h"

bitspray::bitspray() {//creates an unallocated bitspray
    init = St = cs = 0;
    trans = resp = nullptr;
}

bitspray::bitspray(int S) {//create buffer size B, states S
    init = St = cs = 0;
    trans = resp = nullptr;
    create(S);  //allocate a bitspray with B buffer and S states
}

bitspray::bitspray(bitspray &other) {//copy constructor
    init = St = cs = 0;
    trans = resp = nullptr;
    copy(other);  //call the copy routine
}

bitspray::~bitspray() {//destructor
    destroy();    //call the deallocation routines.
}

//utility routines
void bitspray::create(int S) {  //allocate B-sized buffer and S states
    if (St != 0) {   //safety first
        destroy();
    }
    init = 0;              //clear initial action
//  Nb = B;                //assign size
    St = S;                //assign states
//  buf = new int[Nb];     //allocate output buffer
    cs = 0;                //zero the current state
//  hd = tl = 0;           //not reset tag, head and tail at zero
    trans = new int *[St]; //allocate the spine of the transitions
    for (int i = 0; i < St; i++) {
        trans[i] = new int[2];      //allocate binary transitions
    }
    resp = new int *[St];         //allocate spine of responses
    for (int i = 0; i < St; i++) { //allocate Nb1b2 registers
        resp[i] = new int[3];
    }
}

void bitspray::randomize() {  //fill in new random values
    init = (int) lrand48() % 2;  //establish the initial action
    for (int i = 0; i < St; i++) {     //loop over the states
        trans[i][0] = (int) lrand48() % St;    //create zero transition
        trans[i][1] = (int) lrand48() % St;    //create one transition
        resp[i][0] = (int) lrand48() % 2 + 1;  //response length
        for (int j = 1; j <= resp[i][0]; j++)
            resp[i][j] = (int) lrand48() % 2;    //generate response bits
    }
}

void bitspray::copy(bitspray &other) {//copy routine
    if (other.St == 0)destroy();  //safety first!
    if (St != other.St) {         //not the same size
        destroy();
        create(other.St);
    }
    init = other.init;      //copy the initial action
//  Nb = other.Nb;          //copy size
    St = other.St;          //copy number of states
//  buf = new int[Nb];      //allocate output buffer
    cs = other.cs;          //copy the current state
//  hd = other.hd;          //copy the head
//  tl = other.tl;          //copy the tail
//  trans = new int *[St];  //allocate the spine of the transitions
    for (int i = 0; i < St; i++) {//loop over states
//    trans[i] = new int[2];  //alocate binary transitions
        trans[i][0] = other.trans[i][0];  //copy zero transition
        trans[i][1] = other.trans[i][1];  //copy one transition
    }
//  resp = new int *[St];   //allocate spine of responses
    for (int i = 0; i < St; i++) {//loop over responses
//    resp[i] = new int[3]; //allocate Nb1b2 registers
        resp[i][0] = other.resp[i][0];  //copy
        resp[i][1] = other.resp[i][1];  //the
        resp[i][2] = other.resp[i][2];  //response
    }
}

void bitspray::destroy() {//deallocate
    if (St > 0) {//if there is something to deallocate
//    delete[] buf;         //kill the memory buffer
        for (int i = 0; i < St; i++) { //loop over states
            delete[] trans[i];  //delete the transition information
            delete[] resp[i];   //delete the response information
        }
        delete[] trans;       //delete the spine of the transitions
        delete[] resp;        //delete the spine of the responses
    }
    init = St = cs = 0;
    trans = resp = nullptr;
//  Nb = 0;                 //mark as unallocated
}

//use routines
void bitspray::reset(int *Q, int &H, int &T) {//reset to initial state
    if (St > 0) {//if the automata has been properly set up
        cs = 0;            //set to starting state
//    buf[0] = init;     //load initial action into the buffer
        Q[0] = init;
        H = 0;            //head of the buffer queue is zero
        T = 1;            //tail of the buffer queue is one
    } //all reset!
}

void bitspray::next(int *Q, int &H, int &T, int qz) {//update the output queue
    int td; //transition driver buffer

    td = Q[H++];
    if (H == qz)H = 0;  //move the head forward
    for (int i = 0; i < resp[cs][0]; i++) {//put response into the output queue
        int d = resp[cs][i + 1];
        Q[T++] = d;  //transfer character
        if (T == qz)T = 0; //move the tail forward
    }
    cs = trans[cs][td];  //update state
}

//bool bitspray::hasnext(int b) const {
//  return (hd + b < Nb);
//}

//int bitspray::bit() {//get the next bit
//
//  int ret;  //return value
//  int i;    //loop index
//
//  if (Nb > 0 && hasnext(1)) { //intercept insanity
//    ret = buf[hd];       //compute the return value
//    hd++;                //advance head node
//    for (i = 0; i < resp[cs][0]; i++) {//transfer the response to the input
//      buf[tl] = resp[cs][i + 1]; //output->input
//      tl++;
//      if (tl == Nb)tl = 0;  //increment tail mod buffer length
//    }
//    cs = trans[cs][ret];  //update current state
//    return (ret);         //report the bit in question
//  } else return (0);      //empty automata returns zeros
//}

//int bitspray::val(int b) {//return an integer value based on b bits
//
//  int accu;  //value accumulator
//  int i;     //loop index
//
//  if (b < 1 || hasnext(b)) return (0);  //intercept psycho bits
//  accu = bit();          //get the first bit
//  for (i = 0; i < b - 1; i++) {
//    accu *= 2;
//    accu += bit();
//  }  //accumulate the bits
//  return (accu);         //return the integer to the user
//}

//double bitspray::rel(int b) {//return a [0,1) value based on b bits
//
//  int accu;  //value accumulator
//  int d;     //divisor
//  int i;     //loop index
//
//  if (b < 1 || hasnext(b)) return (0);  //intercept psycho bits
//  accu = bit();          //get the first bit
//  d = 2;                 //initialize the divisor
//  for (i = 0; i < b - 1; i++) {
//    d *= 2;
//    accu *= 2;
//    accu += bit();
//  }  //accumulate the bits
//  return (((double) accu) / d);      //return the integer to the user
//}

//Genetics
void bitspray::tpc(bitspray &other) {//two point crossover
    int cp1, cp2;   //crossover points
    int sw;         //swap variable

    if (St != other.St)return; //do not crossover things of different sizes
    cp1 = (int) lrand48() % St;
    cp2 = (int) lrand48() % St;   //generate crossover points
    if (cp1 > cp2) {
        sw = cp1;
        cp1 = cp2;
        cp2 = sw;
    }  //order them correctly
    //cout << "CP1=" << cp1 << " CP2=" << cp2 << endl;
    if (cp1 == 0) {
        sw = init;
        init = other.init;
        other.init = sw;
    }//init follow state 0
    for (int i = cp1; i < cp2; i++) {//loop over crossover area
        sw = trans[i][0];
        trans[i][0] = other.trans[i][0];
        other.trans[i][0] = sw;//swap
        sw = trans[i][1];
        trans[i][1] = other.trans[i][1];
        other.trans[i][1] = sw;//swap
        sw = resp[i][0];
        resp[i][0] = other.resp[i][0];
        other.resp[i][0] = sw;//swap
        sw = resp[i][1];
        resp[i][1] = other.resp[i][1];
        other.resp[i][1] = sw;//swap
        sw = resp[i][2];
        resp[i][2] = other.resp[i][2];
        other.resp[i][2] = sw;//swap
    }//Done crossing over!
}

void bitspray::mutate(int nm) {//mutate a response or transition
    int m;   //mutation type index

    for (int i = 0; i < nm; i++) {
        m = (int) lrand48() % (2 * St + 1);           //select mutation type
        //cout << m << " " << (m-1)/2 << endl;
        if (m == 0) {//flip initial action
            init = 1 - init;
            return;
        }
        m = (m - 1) / 2;  //find state to mutate
        if (lrand48() % 2 == 0) {//mutate transition
            if (lrand48() % 2 == 0) {//mutate first transition
                trans[m][0] = (int) lrand48() % St;  //mutate the transition
            } else {//mutate second transition
                trans[m][1] = (int) lrand48() % St;  //mutate the transition
            }
        } else {//mutate response
            resp[m][0] = (int) lrand48() % 2 + 1;  //generate new response size
            resp[m][1] = (int) lrand48() % 2;      //generate first bit of response
            resp[m][2] = (int) lrand48() % 2;      //generate second bit of response
        }
    }
}

//Input/Output
void bitspray::print(ostream &aus) {//write in human readable form
    aus << "->" << init << endl;
    for (int i = 0; i < St; i++) {//loop over states
        aus << i << ") ";  //print state number
        for (int j = 0; j < resp[i][0]; j++) aus << resp[i][j + 1]; //response
        aus << " 0->" << trans[i][0];  //zero transition
        aus << " 1->" << trans[i][1];  //one transition
        aus << endl;
    }
}

void bitspray::write(ostream &aus) {//write to output stream
//  aus << Nb << endl;          //save buffer size
    aus << St << endl;          //save number of states
    aus << init << endl;        //save the initial action
    for (int i = 0; i < St; i++) {//loop over states
        aus << trans[i][0] << endl;
        aus << trans[i][1] << endl;
        aus << resp[i][0] << endl;    //Save sate details
        aus << resp[i][1] << endl;
        aus << resp[i][2] << endl;
    }
}

void bitspray::read(istream &inp) {//read from input stream
    int B, S;        //code for buffers and states
    char buf[20];   //input buffer
    int i, j;        //loop index variables

//  inp.getline(buf, 19);
//  B = atoi(buf);     //read buffer size
    inp.getline(buf, 19);
    S = atoi(buf);     //read number of states
    create(S);                  //allocate the bitsprayer
    inp.getline(buf, 19);
    init = atoi(buf);  //read the initial action
    for (i = 0; i < St; i++) {//loop over states
        inp.getline(buf, 19);
        trans[i][0] = atoi(buf);  //read state information
        inp.getline(buf, 19);
        trans[i][1] = atoi(buf);  //read state information
        inp.getline(buf, 19);
        resp[i][0] = atoi(buf);   //read state information
        inp.getline(buf, 19);
        resp[i][1] = atoi(buf);   //read state information
        inp.getline(buf, 19);
        resp[i][2] = atoi(buf);   //read state information
    }
}

//void bitspray::getBuf(ostream &out) {
//  out << "Current Buffer: ";
//  for (int i = hd; i < tl; ++i) {
//    out << buf[i];
//  }
//  out << "\tH: " << hd << "\tT: " << tl;
//  out << endl;
//}
