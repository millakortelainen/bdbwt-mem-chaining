#include <stdlib.h>
#include <iostream>
#include "mem.hh"
#include "chaining.hh"
#include "minimizer.hh"
#include "edlib.h"
#include <stdio.h>
#include "io.hh"
#include <math.h>
#include <algorithm>
#include <chrono>
#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"
#include "driver.hh"
#include "util.hh"
using namespace std;

int main(int argc, char *argv[]){

  
  // vector<pair<int, pair<int,int>>> retVec;
  // int n = 10;
  // int *keys_a = new int[n]; //Have to add extra space for 0 index.
  // vector<int> keysa =  {0, 2, 3, 4, 5, 7, 10, 11, 20, 30, 40};
  // vector<int> keys_v = {0, 0,   0,  9, 2, 1 , 2,  7, 0, 5, 6};
  // for(int i = 0; i < n; i++) keys_a[i] = keysa[i];

  // RMaxQTree *tree_a = new RMaxQTree(keys_a,n);

  // for(int ta = 0; ta < n; ta++) tree_a->update(keys_a[ta+1],ta,keys_v[ta]);
  // int ra, rb;
  // tie(ra,rb) = tree_a->query(0,30);
  // cout << "Returned " << keys_a[ra] << ", " << rb << endl;
  // return 0;


  auto conf = readConfiguration(argv[1]);
  string text = conf.text1;
  string text2 = conf.text2;
  //auto edlibConf = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0);
  auto edlibConf = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);
  
  conf.edlibConf = edlibConf;
  if(conf.mode == 3){ //edlib only
    chrono::steady_clock::time_point edlib1_begin = chrono::steady_clock::now();
    EdlibAlignResult result = edlibAlign(text.c_str(), text.size()-1, text2.c_str(), text2.size()-1, edlibConf);
    conf.originalEditDistance = result.editDistance;
    chrono::steady_clock::time_point edlib1_end = chrono::steady_clock::now();
    printf("edit_distance of whole strings = %d, took %ld milliseconds\n", result.editDistance,chrono::duration_cast<chrono::milliseconds>(edlib1_end - edlib1_begin).count());
    edlibFreeAlignResult(result);
    return 0;
  }
  BD_BWT_index<> index((uint8_t*)text.c_str());
  BD_BWT_index<> index2((uint8_t*)text2.c_str());
  conf.index1 = index;
  conf.index2 = index2;
  vector<Interval_pair> Ipairs;
  auto LFmapping1 = mapLF(index, true).first;
  auto LFmapping2 = mapLF(index2, true).first;
  //for (auto p : LFmapping)
  //{ 
    //cout << p.first << ", " << p.second << endl;
  //}
  set<char> alphabet;
  for(auto c : conf.text1){
    alphabet.insert(c);
  }
  conf.alphabet = alphabet;

  Ipairs = computeMemIntervals(conf);
  //for(auto i : Ipairs){
      //cout << i.forward.left << "," << i.reverse.left << "," << i.forward.right-i.forward.left +1<< endl;
      //cout << LFmapping1[i.reverse.left] << endl;
      //cout << LFmapping2[i.forward.left] << endl;
      //cout << "" << endl;
    //}
  //auto chains = computeChains(conf, Ipairs);
  //auto chainintsP = computeChainIntervals(conf, chains, Ipairs);
  //auto chainints = chainintsP.first; 
  //auto absentEdits = computeEditDistancesForAbsentIntervals(conf, chainintsP, Ipairs, conf.verboseEditDistances).first;
  //auto combined = combine_MEM_and_absent_with_editDistances(conf, absentEdits, chainints, conf.printAbsentAndChains);
}


