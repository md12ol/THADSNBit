//
// Created by Michael Dubé on 2022-02-07.
//
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <string>
#include <bitset>
#include "setu.h"
#include "stat.h"

int main(int argc, char *argv[]) {
    int verts = 256;
    cout << "STARTED" << endl;
    graph G(verts);
    G.RNGnm(verts, 2);
    int var_count = 0;
    pair<int, int> var_lens[max_vars];
    vector<int> var_profs[max_vars];
    bitset<dna_len> variants[verts];
    int var_parents[max_vars];
    int epis = 100;
    bool justSIR = false;
    double lengths[epis];
    vector<double> v_lengths;
    v_lengths.reserve(epis);
    int max = 0, len = 0, ttl = 0;
    double var_prob = 0.01;
    int edit_lowB = 4;
    int edit_upB = 12;

//    G.varSIR(0, var_count, var_profs, mde_var, var_parents, var_lens, 0.5);
    for (int j = 0; j < epis; j++) {
        if (!justSIR) {
            G.varSIR(0, var_count, var_profs, variants,
                     var_parents, var_lens, 0.5, edit_lowB, edit_upB, var_prob);
            v_lengths.clear();
            for (int i = 0; i <= var_count; i++) {
                v_lengths.push_back((double) (var_lens[i].second));
            }
            sort(v_lengths.begin(), v_lengths.end());
            lengths[j] = v_lengths.at(var_count);
//            cout<<lengths[j]<<"\t";
            for (int i = 0; i < max_vars; i++) {
                if (!var_profs[i].empty()) {
                    cout << "V " << i << " [";
                    cout << left << setw(2) << var_lens[i].first << " -> " << setw(2) << var_lens[i].second << "]: ";
                    for (int k = 0; k < var_lens[i].first; k++) {
                        cout << "  ";
                    }
                    for (int l = 0; l <= var_lens[i].second - var_lens[i].first; l++) {
                        cout << var_profs[i][l] << " ";
                    }
                    cout << endl;
                }
            }
            cout << endl;
        } else {
            vector<int> prof;
            prof.reserve(100);
            G.SIR(0, max, len, ttl, 0.5, prof);
            lengths[j] = (double) len;
//            cout<<lengths[j]<<"\t";
//            for (int i = 0; i < prof.size(); i++) {
//                cout << prof.at(i);
//            }
//            cout << endl;
        }




        //        cout << "EPI " << j << endl;
//
    }
//    dset D;
//    D.add(lengths, epis);
//    cout << left << setw(10) << D.Rmu();
//    cout << left << setw(12) << D.RCI95();
//    cout << left << setw(10) << D.Rsg();
//    cout << left << setw(8) << D.Rmax() << endl;
    cout << "ENDED" << endl;
}