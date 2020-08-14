#include <stdlib.h>
#include <iostream>
#include "mem.hh"
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
using namespace std;

int main(int argc, char *argv[]){
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

  set<char> alphabet;
  for(auto c : conf.text1){
    alphabet.insert(c);
  }
  conf.alphabet = alphabet;

  Ipairs = computeMemIntervals(conf);
  auto chains = computeChains(conf, Ipairs);
  auto chainintsP = computeChainIntervals(conf, chains, Ipairs);
  auto chainints = chainintsP.first; 
  auto absentEdits = computeEditDistancesForAbsentIntervals(conf, chainintsP, Ipairs, false);
  auto combined = combine_MEM_and_absent_with_editDistances(conf, absentEdits, chainints, false);
}


