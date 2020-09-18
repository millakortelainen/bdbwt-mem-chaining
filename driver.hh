#ifndef BDBWT_MEM_CHAIN_DRIVER_HH
#define BDBWT_MEM_CHAIN_DRIVER_HH
#include "mem.hh"
#include "chaining.hh"
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
    if(conf.verbosity > 2) cout << conf.minimumDepth << ", " << conf.minimizerWindowSize << endl;
#pragma omp parallel sections
    {
#pragma omp section
      {
	mini1 = minimizers(conf.text1, conf.minimumDepth, conf.minimizerWindowSize);
      }
#pragma omp section
      {
	mini2 = minimizers(conf.text2, conf.minimumDepth, conf.minimizerWindowSize);
      }
    }
    if(conf.verbosity > 2) cout << mini1.size() << "\t " << mini2.size() << endl;
    auto mimimems = minimizerTuples(mini1,mini2,conf, false).first;
    if(conf.verbosity > 2) cout << mimimems.size() << endl;
    auto minimems = memifyMinimizers(mimimems, conf);
    Ipairs = returnMemTuplesToIntervals(minimems, false);
    break; 
  }
  case 2: { //hybrid
    vector<pair<string,int>> mini1;
    vector<pair<string,int>> mini2;
    vector<pair<int,int>> SA1, SA2, SA3, SA4;
#pragma omp parallel sections
    {
#pragma omp section
      {
        mini1 = minimizers(conf.text1, conf.minimumDepth, conf.minimizerWindowSize);
        if(conf.verbosity > 2) cout << "mini1 done" << endl;
      }
#pragma omp section
      {
        mini2 = minimizers(conf.text2, conf.minimumDepth, conf.minimizerWindowSize);
        if(conf.verbosity > 2) cout << "mini2 done" << endl;
      }
#pragma omp section
      {
        SA1 = buildSAfromBWT(conf.index1, true);
        if(conf.verbosity > 2) cout << "SA1 done" << endl;
      }
#pragma omp section
      {
        SA2 = buildSAfromBWT(conf.index1, false);
        if(conf.verbosity > 2) cout << "SA2 done" << endl;
      }
#pragma omp section
      {
        SA3 = buildSAfromBWT(conf.index2, true);
        if(conf.verbosity > 2) cout << "SA3 done" << endl;
      }
#pragma omp section
      {
        SA4 = buildSAfromBWT(conf.index2, false);
        if(conf.verbosity > 2) cout << "SA4 done " << endl;
      }
    }
    if(conf.verbosity > 2) cout << "getting muts...";
    auto muts = minimizerTuples(mini1, mini2, conf, true).second;
    mini1 = muts.first;
    mini2 = muts.second;
    if(conf.verbosity > 2)cout << "got muts\t " << mini1.size() << "," << mini2.size() << endl;
    int k = conf.minimumDepth;
    int q = conf.PLCPSparsity_q;
    vector<pair<Interval_pair,string>> set1;
    vector<pair<Interval_pair,string>> set2;

#pragma omp parallel sections
    {
#pragma omp section
      {
        string rev = string(conf.text1.rbegin(),conf.text1.rend());
        auto lcp1 = createPLCP(conf.index1, conf.PLCPSparsity_q, conf.text1, true, SA1, true);
        auto lcp2 = createPLCP(conf.index1, conf.PLCPSparsity_q,      rev  , true, SA2, false);
        auto b1 = partitioning(conf.minimumDepth, conf.index1, lcp1);
        auto b2 = partitioning(conf.minimumDepth, conf.index1, lcp2);
        set1 = minimizerToBWTInterval(b1,b2,mini1,SA1,SA2,conf.index1,lcp1,lcp2, conf.text1);
      }
#pragma omp section
      {
        string rev = string(conf.text2.rbegin(),conf.text2.rend());
        auto lcp1 = createPLCP(conf.index2, conf.PLCPSparsity_q, conf.text2, true, SA3, true);
        auto lcp2 = createPLCP(conf.index2, conf.PLCPSparsity_q,      rev  , true, SA4, false);
        auto b1 = partitioning(conf.minimumDepth, conf.index2, lcp1);
        auto b2 = partitioning(conf.minimumDepth, conf.index2, lcp2);
        set2 = minimizerToBWTInterval(b1,b2,mini2,SA3,SA4,conf.index2,lcp1,lcp2, conf.text2);

      }
    }
    if(conf.verbosity > 2)cout << "got BWT Intervals" << endl;
    if(conf.verbosity > 2)cout << "set1 size: " << set1.size() << ", set2 size: " << set2.size() << endl;
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
    if(conf.verbosity > 1)cout << "seeds size: " << seeds.size() << endl;
    auto bo = bwt_to_int_tuples(conf, seeds);
    Ipairs = returnMemTuplesToIntervals(bo, false);
    break;
  }
  }
  chrono::steady_clock::time_point mems_end = chrono::steady_clock::now();
  if(conf.verbosity > 2) printf("mems took %ld seconds\n", chrono::duration_cast<chrono::seconds>(mems_end - mems_begin).count());
  return Ipairs;
}

vector<pair<int,pair<int,int>>> computeChains(Configuration conf, vector<Interval_pair> Ipairs){
  chrono::steady_clock::time_point chains_begin = chrono::steady_clock::now();
  auto chains = chaining(Ipairs, conf.maxSize);
  chrono::steady_clock::time_point chains_end = chrono::steady_clock::now();
  return chains;
}
pair<vector<Interval_pair>,vector<int>> computeChainIntervals(Configuration conf, vector<pair<int,pair<int,int>>> chains, vector<Interval_pair> Ipairs){
  return chainingOutput(chains, Ipairs, conf);
}
pair<vector<pair<Interval_pair, int>>,int> computeEditDistancesForAbsentIntervals(Configuration conf, pair<vector<Interval_pair>,vector<int>> chainintspair, vector<Interval_pair> Ipairs, bool verbose){
  auto chainints = chainintspair.first;
  auto absent = absentIntervals(chainints, conf.index1, conf.index2); 
  vector<vector<pair<Interval_pair, int>>> absentEdits(omp_get_max_threads());
  vector<pair<Interval_pair,int>> absentRet;
  int totalEditDistance = 0;
  int recombEditDistance = 0;
  if(conf.recombAbsents){
    string reconF = "";
    string reconB = "";
    string appA = "";
    string appB = "";

    for(auto abs : absent){
      cout << abs.toString() << endl;
      if(abs.forward.left >= 0 && abs.forward.right >= 0){
        appA = conf.text1.substr(abs.forward.left, abs.forward.size());
        reconF.append(appA);
        cout << "append A: "<< appA << endl;
      }
      if(abs.reverse.left >= 0 && abs.reverse.right >= 0){
        appB = conf.text2.substr(abs.reverse.left, abs.reverse.size());
        reconB.append(appB);
        cout << "append B: "<< appB << endl;
      }
      // cout << endl;
  }
    if(reconF.length() > 0 && reconB.length() > 0){
    EdlibAlignResult result = edlibAlign(reconF.c_str(), reconF.length()-1,
                                         reconB.c_str(), reconB.length()-1,
                                         conf.edlibConf);
    cout << "Recombined ED: " << result.editDistance << endl;
    recombEditDistance = result.editDistance;
    cout << "Recombined Length: " << reconF.length() << ", " << reconB.length() << endl;
    //cout << reconF;
    edlibFreeAlignResult(result);
    }
  }
#pragma omp parallel for
  for(int i = 0; i < absent.size(); i++){
    int ed;
    if(absent[i].forward.left == -1 && absent[i].forward.right == -1){
      ed = (absent[i].reverse.right - absent[i].reverse.left)+1;
      if(ed > 0){
        totalEditDistance += ed;
      }
    }
    else if(absent[i].reverse.left == -1 && absent[i].reverse.right == -1){
      ed = (absent[i].forward.right - absent[i].forward.left)+1;
      if(ed > 0){
        totalEditDistance += ed;
      }
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
    if(conf.verbosity > 3) cout << "I: " << absent[i].toString() << "..." << setw(50-absent[i].toString().size()-1) << right;
    if(conf.verbosity > 3) cout << setw(50-absent[i].toString().size()-1) << "ED: " << ed << "..."
                     << setw(18-to_string(ed).size()) <<"(ED)/|I|: " << to_string((round((ed / ((double)length+1))*10000)/10000)) << "..."
                     << setw(20-to_string(round((ed / ((double)length+1)))).size()) << "max len: " << (length)+1 << endl;
    absentEdits[omp_get_thread_num()].push_back(make_pair(absent[i], ed));
    if(a != 0 && c == 0){
      totalEditDistance += a;
    }
    if(a == 0 && c != 0){
      totalEditDistance += c;
    }
  }
  for(auto a : absentEdits){
    for(auto b : a){
      totalEditDistance += b.second;
      absentRet.push_back(b);
    }
  }
  sort(absentRet.rbegin(), absentRet.rend(), intervalIntPairSort);
  if(conf.verbosity > 0) cout << "total edit distance: " << totalEditDistance << "/ " << conf.originalEditDistance << endl;
  if(conf.recombAbsents){
    return make_pair(absentRet,recombEditDistance);
  }
  return make_pair(absentRet,totalEditDistance);
}

vector<pair<Interval_pair, int>> combine_MEM_and_absent_with_editDistances(Configuration conf, vector<pair<Interval_pair, int>> absentEdits, vector<Interval_pair> chainints, bool verbose){
  int ci = 0;
  int ai = 0;
  bool absentFirst = true;
  vector<pair<Interval_pair, int>> combinedED;
  if(chainints[0].forward.left > absentEdits[0].first.forward.left){
    absentFirst = false;
  }
  if(conf.printAbsentAndChains){
    verbose = true;
  }
  while(ci < chainints.size() || ai < absentEdits.size()){ // print all in order
    auto cint = chainints[ci];
    auto aint = absentEdits[ai].first;

    if(absentFirst){
      if(absentEdits.size()-1 >= ai){
        if(verbose) cout << "A: " << absentEdits[ai].first.toString() << endl;
        combinedED.push_back(absentEdits[ai]);
      }
      if(chainints.size()-1 >= ci){
        if(verbose) cout << "C: " << chainints[ci].toString() << endl;
        combinedED.push_back(make_pair(chainints[ci], 0));
      }
    }else{
      if(chainints.size()-1 >= ci){
        if(verbose) cout << "C: " << chainints[ci].toString() << endl;
        combinedED.push_back(make_pair(chainints[ci], 0));
      }
      if(absentEdits.size()-1 >= ai){
        if(verbose) cout << "A: " << absentEdits[ai].first.toString() << endl;
        combinedED.push_back(absentEdits[ai]);
      }
    }
    ci++;
    ai++;
  }
  return combinedED;
}
#endif
