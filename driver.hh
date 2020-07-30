#ifndef BDBWT_MEM_CHAIN_DRIVER_HH
#define BDBWT_MEM_CHAIN_DRIVER_HH
#include "mem.hh"
#include "minimizer.hh"
#include "edlib.h"
#include "util.hh"
using namespace std;

vector<Interval_pair> computeMemIntervals(Configuration conf){
  chrono::steady_clock::time_point mems_begin = chrono::steady_clock::now();
  vector<Interval_pair> Ipairs;
  switch(conf.mode){
  case 0: { //BWT Only
    auto bo = bwt_to_int_tuples(conf);
    Ipairs = returnMemTuplesToIntervals(bo, false);
    break;
  }
  case 1: { //Minimizer Only
    vector<pair<string,int>> mini1;
    vector<pair<string,int>> mini2;
    cout << conf.minimumDepth << ", " << conf.minimizerWindowSize << endl;
#pragma opm parallel sections
    {
#pragma opm section
      {
	mini1 = minimizers(conf.text1, conf.minimumDepth, conf.minimizerWindowSize);
      }
#pragma opm section
      {
	mini2 = minimizers(conf.text2, conf.minimumDepth, conf.minimizerWindowSize);
      }
    }
    cout << mini1.size() << "\t " << mini2.size() << endl;
    auto mimimems = minimizerTuples(mini1,mini2, false, conf.text1, conf.text2).first;
    cout << mimimems.size() << endl;
    auto minimems = memifyMinimizers(mimimems, conf.text1, conf.text2);
    Ipairs = returnMemTuplesToIntervals(minimems, false);
    break;
  }
  case 2: { //hybrid
    vector<pair<string,int>> mini1;
    vector<pair<string,int>> mini2;
    vector<pair<int,int>> SA1, SA2, SA3, SA4;
#pragma opm parallel sections
    {
#pragma opm section
      {
        mini1 = minimizers(conf.text1, conf.minimumDepth, conf.minimizerWindowSize);
	cout << "mini1 done" << endl;
      }
#pragma opm section
      {
        mini2 = minimizers(conf.text2, conf.minimumDepth, conf.minimizerWindowSize);
	cout << "mini2 done" << endl;
      }
#pragma opm section
      {
        SA1 = buildSAfromBWT(conf.index1, true);
	cout << "SA1 done" << endl;
      }
#pragma opm section
      {
        SA2 = buildSAfromBWT(conf.index1, false);
	cout << "SA2 done" << endl;
      }
#pragma opm section
      {
        SA3 = buildSAfromBWT(conf.index2, true);
	cout << "SA3 done" << endl;
      }
#pragma opm section
      {
        SA4 = buildSAfromBWT(conf.index2, false);
	cout << "SA4 done " << endl;
      }
    }
    cout << "getting muts...";
    auto muts = minimizerTuples(mini1, mini2, true, conf.text1, conf.text2).second;
    mini1 = muts.first;
    mini2 = muts.second;
    cout << "got muts\t " << mini1.size() << "," << mini2.size() << endl;
    int k = conf.minimumDepth;
    int q = conf.PLCPSparsity_q;
    vector<pair<Interval_pair,string>> set1;
    vector<pair<Interval_pair,string>> set2;

#pragma opm parallel sections
    {
#pragma opm section
      {
	bool type = true;
	if(type){
	string rev = string(conf.text1.rbegin(),conf.text1.rend());
	auto lcp1 = createPLCP(conf.index1, conf.PLCPSparsity_q, conf.text1, true, SA1, true);
	auto lcp2 = createPLCP(conf.index1, conf.PLCPSparsity_q,      rev  , true, SA2, false);
	auto b1 = partitioning(conf.minimumDepth, conf.index1, lcp1);
	auto b2 = partitioning(conf.minimumDepth, conf.index1, lcp2);
	set1 = minimizerToBWTInterval(b1,b2,mini1,SA1,SA2,conf.index1,lcp1,lcp2, conf.text1);
	}else{
	//set1 = minimizerToBWTIntervalV2(mini1,conf.minimumDepth, SA1, SA2, conf.index1.get_global_c_array(), conf.text1);
	}
      }
#pragma opm section
      {
	bool type = true;
	if(type){
	string rev = string(conf.text2.rbegin(),conf.text2.rend());
	auto lcp1 = createPLCP(conf.index2, conf.PLCPSparsity_q, conf.text2, true, SA3, true);
	auto lcp2 = createPLCP(conf.index2, conf.PLCPSparsity_q,      rev  , true, SA4, false);
	auto b1 = partitioning(conf.minimumDepth, conf.index2, lcp1);
	auto b2 = partitioning(conf.minimumDepth, conf.index2, lcp2);
	set2 = minimizerToBWTInterval(b1,b2,mini2,SA3,SA4,conf.index2,lcp1,lcp2, conf.text2);
	}else{
	//set2 = minimizerToBWTIntervalV2(mini2,conf.minimumDepth, SA3, SA4, conf.index2.get_global_c_array(), conf.text2);
      }
    }
    }
    cout << "got BWT Intervals" << endl;
    cout << "set1 size: " << set1.size() << ", set2 size: " << set2.size() << endl;
    set<tuple<Interval_pair,Interval_pair,int>> seeds;
    // sort(set1.begin(),set1.end(),intervalSort);
    // sort(set2.begin(),set2.end(),intervalSort);
    int j = 0;
    for(int i = 0; i < set1.size(); i++){
      if(i < set2.size()){
        if(set1[i].second.compare(set2[i].second) != 0) continue;
        //cout << set1[i].second << set1[i].first.toString() << endl;
        //cout << set2[i].second << set2[i].first.toString() << endl;
        //cout << endl;
        seeds.insert(make_tuple(set1[i].first, set2[i].first, set1[i].second.length()));
      }
    }
    cout << "seeds size: " << seeds.size() << endl;
    auto bo = bwt_to_int_tuples(conf, seeds);
    Ipairs = returnMemTuplesToIntervals(bo, false);
      
    break;
  }
  }
  chrono::steady_clock::time_point mems_end = chrono::steady_clock::now();
  printf("mems took %ld seconds\n", chrono::duration_cast<chrono::seconds>(mems_end - mems_begin).count());
  return Ipairs;
}


vector<pair<int,pair<int,int>>> computeChains(Configuration conf, vector<Interval_pair> Ipairs){
  // for(auto asd : Ipairs){
//     int i,j,k;
//     i = asd.forward.left;
//     j = asd.reverse.left;
//     k = asd.forward.right-asd.forward.left+1;

//     cout << "Ipair: " << i << "," << j << "," << k << endl; 
// }
  chrono::steady_clock::time_point chains_begin = chrono::steady_clock::now();
  auto chains = chaining(Ipairs, conf.maxSize);
  chrono::steady_clock::time_point chains_end = chrono::steady_clock::now();
  return chains;
}
pair<vector<Interval_pair>,vector<int>> computeChainIntervals(Configuration conf, vector<pair<int,pair<int,int>>> chains, vector<Interval_pair> Ipairs){
  return chainingOutput(chains, Ipairs, conf.text1, conf.text2);
}
vector<pair<Interval_pair, int>> computeEditDistancesForAbsentIntervals(Configuration conf, pair<vector<Interval_pair>,vector<int>> chainintspair, vector<Interval_pair> Ipairs){
  auto chainints = chainintspair.first;
  auto absent = absentIntervals(chainints, conf.index1, conf.index2); 
  vector<pair<Interval_pair, int>> absentEdits;
  int totalEditDistance = 0;
  for(int i = 0; i < absent.size(); i++){
    int ed;
    if(absent[i].forward.left == -1){
      ed = absent[i].reverse.right - absent[i].reverse.left+1;
    }
    else if(absent[i].reverse.left == -1){
      ed = absent[i].forward.right - absent[i].forward.left+1;
    }else{
      EdlibAlignResult result = edlibAlign(conf.text1.substr(absent[i].forward.left, absent[i].forward.size()).c_str(), absent[i].forward.size(),
                                           conf.text2.substr(absent[i].reverse.left, absent[i].reverse.size()).c_str(), absent[i].reverse.size(),
                                           conf.edlibConf);
      ed = result.editDistance;
      edlibFreeAlignResult(result);
    }
    int a = absent[i].forward.left;
    int b = absent[i].forward.right;
    int c = absent[i].reverse.left;
    int d = absent[i].reverse.right;
    int length = ((b-a) > (d - c))? (b-a) : (d-c);
    
    cout << "I: " << absent[i].toString() << "..." << setw(50-absent[i].toString().size()-1) << right;
    
    cout << setw(50-absent[i].toString().size()-1) << "ED: " << ed << "..."
         << setw(18-to_string(ed).size()) <<"(ED)/|I|: " << to_string((round((ed / ((double)length+1))*10000)/10000)) << "..."
         << setw(20-to_string(round((ed / ((double)length+1)))).size()) << "max len: " << (length)+1 << endl;
    
    totalEditDistance += ed;
    absentEdits.push_back(make_pair(absent[i], ed));
  }
  cout << "total edit distance became: " << totalEditDistance << "/ " << conf.originalEditDistance << endl;
  return absentEdits;
}

vector<pair<Interval_pair, int>> combine_MEM_and_absent_with_editDistances(Configuration conf, vector<pair<Interval_pair, int>> absentEdits, vector<Interval_pair> chainints){
  int ci = 0;
  int ai = 0;
  bool absentFirst = true;
  vector<pair<Interval_pair, int>> combinedED;
  if(chainints[0].forward.left > absentEdits[0].first.forward.left){
    absentFirst = false;
  }
  while(ci < chainints.size() || ai < absentEdits.size()){ // print all in order
    auto cint = chainints[ci];
    auto aint = absentEdits[ai].first;

    if(absentFirst){
      if(absentEdits.size()-1 >= ai){
	cout << "A: " << absentEdits[ai].first.toString() << endl;
	combinedED.push_back(absentEdits[ai]);
      }
      if(chainints.size()-1 >= ci){
	cout << "C: " << chainints[ci].toString() << endl;
	combinedED.push_back(make_pair(chainints[ci], 0));
      }
    }else{
      if(chainints.size()-1 >= ci){
	cout << "C: " << chainints[ci].toString() << endl;
	combinedED.push_back(make_pair(chainints[ci], 0));
      }
      if(absentEdits.size()-1 >= ai){
	cout << "A: " << absentEdits[ai].first.toString() << endl;
	combinedED.push_back(absentEdits[ai]);
      }
    }
    ci++;
    ai++;
  }
  return combinedED;
}
#endif
