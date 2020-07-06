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
//#include "util.hh"
using namespace std;

int main(int argc, char *argv[]){
  string text;
  string text2;
  switch(0){
  case 0: {
    if(argc > 4){
      auto fileinput = readInputFromFasta(argv[3]);
      auto fileinput2 = readInputFromFasta(argv[4]);
      text = fileinput.at(0);
      text2 = fileinput2.at(0);
    }else{
      auto fileinput = readInputFromFasta(argv[3]);
      text = fileinput.at(0);
      text2 = fileinput.at(1);
    }
    break;
  }
  case 1: {
    text =  "CAATTTAAGGCCCGGG";
    text2 = "CAAAGTAAGGCCCTCC";
    //text  = "GTGCGTGATCATCATTT";
    //text2 = "AGTGCAAAGTGATTACC";
    break;
  }
  case 2:{
    text  = "ASDKISSAIKALAS";
    text2 = "ASDKASSAIAKALA";
    break;
  }
  case 3: {
    text  = "GTGCGTGTTCATCATTT";
    text2 = "GTGCGTGATCATCATTT";
    break;
  }
  case 4: {
    text  = "GTGCGTGATCATCATTT";
    text2 = "AGTGCGCGTGACATCTT";
    break;
  }
  case 5: {
    text = "CAATTTAAGGCCCGGGGTGCGTGATCATCATTTGTGCGTGTTCATCATTTGTGCGTGATCATCATTT";
    text2= "CAAAGTAAGGCCCTCCAGTGCAAAGTGATTACCGTGCGTGATCATCATTTAGTGCGCGTGACATCTT";
    break;
  }
  case 500: {
    text = "gttcaccatttaaataatcttcaatatcaacacgcgaagctcgcttgcagggatgaactgaatagacctgtttactccggaaaagcaagactatcctggtgctgatgctacggtacattgttcttggcacgattacggactattcacactgaatccgggtggggagggccttatggacacgtaatatgcgcgtactggttggcgttgtagacgcgcaacttcatcgataatctgactgcctgacaagctaccagcaatacgttactccatcccgctatcctcggtactgcttgcggtgtcaccccgttaagtgacgtcctgttcgcggctaggctacgagttgcgttaatgcactctgaatcagaattccgcagcgttaagctggcttcaccagcgtcttcggtctgacttaaacctactcccgacatttctacagtgactactgtgtacgccccacgaagtcaaccccgagctacacctaaccggcctccagcactgcc";
    text2 = "aattcgaatagttagctgacgtacgacatgttaccttaataatataactggtgtccgcgactgagtgctctcctacctcccacgagcctcaggaaaaacgtctttaaatctctacccggagctgtttaaggggaagccaactcgaacctagcagggcattaaatttgtattgcaccaaaacgaccggcttaacattccgtgtctcactggacggaaaaccaacctaagcagtatttggcctcctggtaggcgaaccatctacggtggaccgtataatcggactaaccggcaggtttacacttcgcaatgctacgctgcccagggccgggcccccagtaggtttgcactgtagagggagggccggagtgtatcccccatcggtaactctacatatgcgcaagccgccctgggcaagatcccatcccactcgtgtggctctcgcgccgggtggattgtacgatcggaatcctctggggacgcgcgttcagtaacttcgctta";
    break;
  }
  }
  //auto edlibConf = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0);
  auto edlibConf = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0);
  BD_BWT_index<> index((uint8_t*)text.c_str());
  BD_BWT_index<> index2((uint8_t*)text2.c_str());
  vector<Interval_pair> Ipairs;
  int originalEditDistance = -1;

  if(false){
    chrono::steady_clock::time_point edlib1_begin = chrono::steady_clock::now();
    EdlibAlignResult result = edlibAlign(text.c_str(), text.size()-1, text2.c_str(), text2.size()-1, edlibConf);
    originalEditDistance = result.editDistance;
    chrono::steady_clock::time_point edlib1_end = chrono::steady_clock::now();
    printf("edit_distance of whole strings = %d, took %ld milliseconds\n", result.editDistance,chrono::duration_cast<chrono::milliseconds>(edlib1_end - edlib1_begin).count());
    edlibFreeAlignResult(result);
  }
  
  auto depthargmin = 3; // Default minimum depth
  auto loga = ceil(log10(index.size()-1) / log10(4)+1)+log10(text.size()-1);
  auto deptharg = loga;
  
  cout << deptharg << endl;
  if(argc > 1 && argv[1] != 0){
    minimumDepth = strtol(argv[1],NULL,10);
  }else{
    minimumDepth = (depthargmin > deptharg)? depthargmin : deptharg; 
  }
  int minimizerWindowSize = minimumDepth-1;
  int computationMode = strtol(argv[2],NULL,10);

  Configuration conf;
  conf.mode = computationMode;
  conf.text1 = text;
  conf.text2 = text2;
  conf.index1 = index;
  conf.index2 = index2;
  conf.minimumDepth = minimumDepth;
  conf.minimizerWindowSize = minimizerWindowSize;
  conf.PLCPSparsity_q = 1;
  conf.edlibConf = edlibConf;
  conf.originalEditDistance = originalEditDistance;

  Ipairs = computeMemIntervals(conf);
  auto chains = computeChains(conf, Ipairs);
  auto chainintsP = computeChainIntervals(conf, chains, Ipairs);
  auto chainints = chainintsP.first;
  
  auto absentEdits = computeEditDistancesForAbsentIntervals(conf, chainintsP, Ipairs);
  auto combined = combine_MEM_and_absent_with_editDistances(conf, absentEdits, chainints);
  // //  int averageED = totalEditDistance / absentEdits.size()-1;
  // int ra = combinedED.at(0).first.forward.right;
  // int rb = combinedED.at(0).first.forward.right;
  // int tempED = 0;
  // for(auto e : combinedED){ // print all in order
  //   auto a = e.first.forward.left;
  //   auto b = e.first.forward.right;
  //   auto c = e.first.reverse.left;
  //   auto d = e.first.reverse.right;
  //   auto ed = e.second;
  //    if(ed > averageED){
  //     Interval_pair ip = Interval_pair(b+1,ra+1,d-1,rb-1);
  //     ra = a;
  //     rb = c;
      
  //     int length = (ip.forward.size() > ip.reverse.size())? ip.forward.size() : ip.reverse.size();
  //     cout << ip.toString() << "ed: " << tempED << " total lenght: " << length << endl;
  //     tempED = 0;
  //   }else{
  //     tempED += ed;
  //   }
  // }
}
  

