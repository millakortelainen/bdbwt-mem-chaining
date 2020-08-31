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
#include <sstream>
#include <random>

string convertTime(chrono::milliseconds ms){
  int milli = ms.count();
  int hours = milli / 3600000;
  milli = milli-(hours*3600000);
  int minutes = milli / 60000;
  milli = milli-(minutes*60000);
  int seconds = milli / 1000;
  milli = milli-(seconds*1000);

  string shours = to_string(hours);
  string smilli = to_string(milli);
  string sminutes = to_string(minutes);
  string sseconds = to_string(seconds);

  if(shours.length() < 2){
    shours = "0"+shours;
  }
  while(smilli.length() < 3){
    smilli = "0"+smilli;
  }
  if(sminutes.length() < 2){
    sminutes = "0"+sminutes;
  }
  if(sseconds.length() < 2){
    sseconds = "0"+sseconds;
  }
  ostringstream oss;
  oss << shours << "h" << sminutes << "min" << sseconds << "s" << smilli << "ms";
  string var = oss.str();

  return var;
}

Configuration swapRandomSymbol(Configuration conf, int num){
  vector<char> s;
  for(auto a : conf.alphabet){
    s.push_back(a);
  }
  random_device rd;
  mt19937 gen(rd());
  mt19937 gen2(rd());
  uniform_int_distribution<> distrib(0,conf.text2.length()-1);
  uniform_int_distribution<> distrib2(0,s.size()-1);
  for(int k = 0; k < num; k++){
    int i = distrib(gen);
    int j = distrib2(gen2);
    while(s[j] == conf.text2.at(i)){
      j = distrib2(gen2);
    }
    conf.text2.at(i) = s[j];
  }
  return conf;
}

pair<pair<chrono::milliseconds, int>,pair<chrono::milliseconds, int>> compute(Configuration conf){
  auto edlibConf = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);
  conf.edlibConf = edlibConf;
  chrono::steady_clock::time_point edlib_begin = chrono::steady_clock::now();
  EdlibAlignResult result = edlibAlign(conf.text1.c_str(), conf.text1.size()-1, conf.text2.c_str(), conf.text2.size()-1, edlibConf);
  conf.originalEditDistance = result.editDistance;
  chrono::steady_clock::time_point edlib_end = chrono::steady_clock::now();
  auto edlibTime = chrono::duration_cast<chrono::milliseconds>(edlib_end - edlib_begin);
  printf("Edlib edit_distance of whole strings = %d, took %s\n", result.editDistance, convertTime(edlibTime).c_str());

  chrono::steady_clock::time_point chaining_begin = chrono::steady_clock::now();
  BD_BWT_index<> index((uint8_t*)conf.text1.c_str());
  BD_BWT_index<> index2((uint8_t*)conf.text2.c_str());
  conf.index1 = index;
  conf.index2 = index2;

  auto Ipairs = computeMemIntervals(conf);
  auto chains = computeChains(conf, Ipairs);
  auto chainintsP = computeChainIntervals(conf, chains, Ipairs);
  auto chainints = chainintsP.first;
  auto absentEdits = computeEditDistancesForAbsentIntervals(conf, chainintsP, Ipairs, conf.verboseEditDistances);
  //auto combined = combine_MEM_and_absent_with_editDistances(conf, absentEdits.first, chainints, conf.printAbsentAndChains);
  chrono::steady_clock::time_point chaining_end = chrono::steady_clock::now();
  auto chainingTime = chrono::duration_cast<chrono::milliseconds>(chaining_end - chaining_begin);
  printf("Chaining edit_distance of whole strings = %d, took %s\n", absentEdits.second, convertTime(chainingTime).c_str());

  auto el = make_pair(edlibTime, result.editDistance);
  auto ch = make_pair(chainingTime, absentEdits.second);

  edlibFreeAlignResult(result);
  return make_pair(el,ch);
}
int main(int argc, char *argv[]){
  int editsCount = 0;
  int editAmount = atoi(argv[2]);
  auto conf = readConfiguration(argv[1]);
  set<char> alphabet;
  for(auto c : conf.text2){
    alphabet.insert(c);
  }
  conf.alphabet = alphabet;
  vector<pair<chrono::milliseconds,chrono::milliseconds>> timeResults;
  vector<pair<int,int>> editResults;


  for(int i = 0; i < atoi(argv[3]); i++){
    cout << "total Edits done: " << editsCount << endl;
    auto c = compute(conf);
    auto diff = c.first.first - c.second.first;
    cout << "edlibTime - chainingTime\t " << diff.count() << "ms." << endl;
    conf = swapRandomSymbol(conf, editAmount);
    editsCount += editAmount;
    cout << endl;
    timeResults.push_back(make_pair(c.first.first,c.second.first));
    editResults.push_back(make_pair(c.first.second, c.second.second));
  }
  cout << "Edlib-Time \t Edlib-Edit \t Chaining-Time \t Chaining-Edit" << endl;
  for(int i = 0; i < timeResults.size(); i++){
    cout << convertTime(timeResults[i].first) << "\t";
    cout << editResults[i].first << "\t";
    cout << convertTime(timeResults[i].second) << "\t";
    cout << editResults[i].second << "\t";
    cout << endl;
  }
  return 0;
}
