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
//#include "util.hh"
using namespace std;
void naiveOutput(BD_BWT_index<> index, BD_BWT_index<> index2, vector<tuple<int,int,int>> memVector, string text1 = "", string text2 = "",bool verbose = false){
  auto retSA = buildSAfromBWT(index); //RetSA builds SA array for the given text from it's BWT transform without having to use the extra space from permutating whole original text.
  auto retSA2 = buildSAfromBWT(index2);

  //Naive returning of the intervals
  for(auto a : memVector){
    int i, j, depth;
    tie(i,j,depth) = a;
    int begin_i= retSA[i]+1;
    int begin_j= retSA2[j]+1;
    
    //Handling the specific special case when SA^S[i] = text.size(); We use index.size()-1 instead so wouldn't need to keep the text in memory at all.
    if(begin_i >= index.size()-1){
      begin_i = index.size()-retSA[i]-1; //Index.size() will always be text.size()+1 due to the added END marker. Furthermore, we need to minus one to get proper offset from general case of begin_i = retSA[i]+1;
    }
    //Handling the specific special case when SA^T[j] = text2.size(); We use index2.size()-1 instead so wouldn't need to keep the text in memory at all.
    if(begin_j >= index2.size()-1){
      begin_j = index2.size()-retSA2[j]-1; //Analogously to above. 
    }
    int end_i = begin_i+depth-1;
    int end_j = begin_j+depth-1;
    
    cout << "Given tuple (i,j,d): "<< "(" << i << ","<<j <<"," << depth << ")" << "\n";
    
    cout << "S: ["<< begin_i <<","<< end_i <<"]" << "-->\t";
    if(text1.size() > 1 && verbose){
      for(int b = begin_i; b <= end_i; b++){
	cout << text1[b]; //For verificiation of results
      }
      cout << "\n";
    }
    if(text2.size() > 1 && verbose){
      cout << "T: ["<< begin_j <<","<< end_j <<"]" << "-->\t";
      for(int b = begin_j; b <= end_j; b++){
	cout << text2[b]; //For verification of results
      }
      cout << "\n";
    }
  }
}
vector<tuple<int,int,int>> batchOutput(BD_BWT_index<> index, BD_BWT_index<> index2, vector<tuple<int,int,int>> memVector, bool verbose = false){
  std::vector<struct occStruct> Ipairs;
  std::vector<struct occStruct> Ipairs2;
  std::vector<bool> marked1(index.size(),false);
  std::vector<bool> marked2(index2.size(),false);
  vector<tuple<int,int,int>> retVector;
  int p = 0;
  
  for(auto m : memVector){
    struct occStruct newOcc;
    struct occStruct newOcc2;
    int i,j,d;
    tie(i,j,d) = m;
    
    newOcc.key  = i;
    newOcc2.key = j;
    newOcc.primary  = p;
    newOcc2.primary = p;
    p++;
    
    marked1[i] = true;
    marked2[j] = true;
    Ipairs.push_back(newOcc);
    Ipairs2.push_back(newOcc2);
  }
  
  auto bl1 = batchLocate(Ipairs,marked1,index);
  auto bl2 = batchLocate(Ipairs2,marked2,index2);
  int maxkey = -1;
  for(int k = 0; k < memVector.size(); k++){
    int i,j,d;
    tie(i,j,d) = memVector[bl1[k].primary];
    int ik = bl1[k].key+1;
    int jk = bl2[k].key+1;
    if(ik+1 >= index.size()-1){ //Special case when interval begins at the beginning, SA[i]=text.size()
      ik = index.size()-ik;
    }
    if(jk+1 >= index2.size()-1){ //Special case when interval begins at the beginning, SA[i]=text.size()
      jk = index2.size()-jk;
    }
    retVector.push_back(make_tuple(ik, jk, d));
  }

  if(verbose){
    for(auto k : retVector){
      int i,j,d;
      tie(i,j,d) = k;
      cout << "triple: (" << i <<","<< j <<","<< d <<")\n";
    }
  }
  return retVector;
}
pair<vector<Interval_pair>,vector<int>> chainingOutput(vector<pair<int,pair<int,int>>> chains, vector<Interval_pair> Ipairs){
  int maxIndex = 0;
  int maxVal = 0;
  for(int i = 0; i < chains.size(); i++){
    int tempVal = chains[i].first;
    if(tempVal > maxVal){
      maxVal = tempVal;
      maxIndex = i;
    }
  }
  if(false){ //printing raw chains
    for(int i = 0; i < chains.size(); i++){
      cout <<"Chain["<< i << "]: " << chains[i].first << "," << chains[i].second.first << ":"<<chains[i].second.second << "\n";
    }
  }
  vector<Interval_pair> chainIntervals;
  vector<int> symcov;
  chainIntervals.push_back(Ipairs.at(chains.at(maxIndex).second.second)); //pushing first chain where we begin the traceback.
  symcov.push_back(chains.at(maxIndex).first);
  int last = -1;
  int i = chains[maxIndex].second.first;
  for(int j = chains.size()-1; j >= 0; j--){
    if(i < 0){
      break;
    }
    if(chains.at(i).second.second >= 0){
      auto I = Ipairs.at(chains.at(i).second.second);
      if(chainIntervals.size() > 0 && chains.size() > 1 &&
	 chainIntervals.at(chainIntervals.size()-1).forward.left >= I.forward.left &&
	 chainIntervals.at(chainIntervals.size()-1).reverse.left >= I.reverse.left){ //Ensuring (weak) precedence

	symcov.push_back(chains.at(i).first);
	chainIntervals.push_back(Ipairs.at(chains.at(i).second.second));
	last = i;
	i = chains[i].second.first;
	if(last == 0){ //would print out same index again => chain is done.
	  break;
	}
      }
    }
  }
  int count = 0;
  for(auto c : chainIntervals){
    cout << "Chain["<< count <<"]: "<< c.toString() <<"\t\t symcov:" << symcov.at(count) << endl;
    count++;
  }
  return (make_pair(chainIntervals, symcov));
}
vector<tuple<int,int,int>> bwt_to_int_tuples(BD_BWT_index<> index, BD_BWT_index<> index2){
  vector<tuple<int,int,int>> mems;
  bool threadedBWT = true;
  if(threadedBWT){
    set<tuple<int,int,int>> memFilter;
    cout << "finding mems between indexes...";
    auto enumLeft = enumerateLeft(index,Interval_pair(0, index.size()-1, 0, index.size()-1));
    if(enumLeft.at(0) == BD_BWT_index<>::END){
      enumLeft.erase(enumLeft.begin());
    }
    vector<vector<tuple<int,int,int>>> memThreads(omp_get_max_threads());
    for(auto i : enumLeft){
      cout << i << endl;
    }
#pragma omp parallel for
    for(int i = 0; i < enumLeft.size(); i++){
      auto retMem = bwt_mem2(index, index2, enumLeft.at(i));
      memThreads[omp_get_thread_num()].insert(memThreads[omp_get_thread_num()].end(), retMem.begin(), retMem.end());
    }
      cout << "done finding mems in threads, collapsing";
      for(auto a : memThreads){
	for(auto b : a){
	  memFilter.insert(b);
	}
      }
      mems.insert(mems.end(), memFilter.begin(), memFilter.end());
  }
  else{
    mems = bwt_mem2(index,index2, BD_BWT_index<>::END);
  }

  cout << "found mems," << mems.size() << "...";
  sort(mems.begin(), mems.end(), memSort); //Proper sorting of the tuples with priority order of i --> d --> j

  if(mems.size() == 0){
    cout << "could not find significiantly large enough MEMS " << endl;
    return mems;
  }
  auto bo = batchOutput(index, index2, mems, false);
  cout << "batchOutput into SA indices done" << "...";
  sort(bo.begin(), bo.end(), memSort); //overall Speed increase
  return bo;
}

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
    text  = "GTGCGTGATCATCATTTA";
    text2 = "AGTGCGCGTGACATCTTT";
    break;
  }
  case 500: {
    text = "gttcaccatttaaataatcttcaatatcaacacgcgaagctcgcttgcagggatgaactgaatagacctgtttactccggaaaagcaagactatcctggtgctgatgctacggtacattgttcttggcacgattacggactattcacactgaatccgggtggggagggccttatggacacgtaatatgcgcgtactggttggcgttgtagacgcgcaacttcatcgataatctgactgcctgacaagctaccagcaatacgttactccatcccgctatcctcggtactgcttgcggtgtcaccccgttaagtgacgtcctgttcgcggctaggctacgagttgcgttaatgcactctgaatcagaattccgcagcgttaagctggcttcaccagcgtcttcggtctgacttaaacctactcccgacatttctacagtgactactgtgtacgccccacgaagtcaaccccgagctacacctaaccggcctccagcactgcc";
    text2 = "aattcgaatagttagctgacgtacgacatgttaccttaataatataactggtgtccgcgactgagtgctctcctacctcccacgagcctcaggaaaaacgtctttaaatctctacccggagctgtttaaggggaagccaactcgaacctagcagggcattaaatttgtattgcaccaaaacgaccggcttaacattccgtgtctcactggacggaaaaccaacctaagcagtatttggcctcctggtaggcgaaccatctacggtggaccgtataatcggactaaccggcaggtttacacttcgcaatgctacgctgcccagggccgggcccccagtaggtttgcactgtagagggagggccggagtgtatcccccatcggtaactctacatatgcgcaagccgccctgggcaagatcccatcccactcgtgtggctctcgcgccgggtggattgtacgatcggaatcctctggggacgcgcgttcagtaacttcgctta";
    break;
  }
  }
  //chrono::steady_clock::time_point _begin = chrono::steady_clock::now();
  //chrono::steady_clock::time_point _end = chrono::steady_clock::now();
  //printf(", took %ld seconds\n", result.editDistance,chrono::duration_cast<chrono::seconds>(_end - _begin).count());
  
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
  
  chrono::steady_clock::time_point mems_begin = chrono::steady_clock::now();
  switch(computationMode){
  case 0: { //BWT Only
    auto bo = bwt_to_int_tuples(index, index2);
    Ipairs = returnMemTuplesToIntervals(bo, false);
    break;}
  case 1: { //Minimizer Only
    vector<pair<string,int>> mini1;
    vector<pair<string,int>> mini2;
#pragma opm parallel sections
    {
#pragma opm section
      {
	mini1 = minimizers(text,minimumDepth,minimizerWindowSize);
      }
#pragma opm section
      {
	mini2 = minimizers(text2,minimumDepth,minimizerWindowSize);
      }
    }
    auto mimimems = minimizerTuples(mini1,mini2);
    cout << endl;
    auto minimems = memifyMinimizers(mimimems, text, text2);
    Ipairs = returnMemTuplesToIntervals(minimems, false);
    break;}
  }
  chrono::steady_clock::time_point mems_end = chrono::steady_clock::now();
  printf("mems took %ld seconds\n", chrono::duration_cast<chrono::seconds>(mems_end - mems_begin).count());  
  
  cout << endl;
  // for(int i = 0; i < Ipairs.size(); i++){
  //   cout <<"Ipairs["<< i << "]: " << Ipairs[i].toString() << "\n";
  // }

  chrono::steady_clock::time_point chains_begin = chrono::steady_clock::now();
  auto chains = chaining(Ipairs, text2.size());
  chrono::steady_clock::time_point chains_end = chrono::steady_clock::now();
  printf("chains took %ld seconds\n", chrono::duration_cast<chrono::seconds>(chains_end - chains_begin).count());  
  auto chainints = chainingOutput(chains, Ipairs).first;
  auto absent = absentIntervals(chainints, index, index2);
  
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
      EdlibAlignResult result = edlibAlign(text.substr(absent[i].forward.left, absent[i].forward.size()).c_str(), absent[i].forward.size(),
					   text2.substr(absent[i].reverse.left, absent[i].reverse.size()).c_str(), absent[i].reverse.size(),
					   edlibConf);
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

  cout << "total edit distance became: " << totalEditDistance << "/ " << originalEditDistance << endl;
  return 0;
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
  int averageED = totalEditDistance / absentEdits.size()-1;
  int ra = combinedED.at(0).first.forward.right;
  int rb = combinedED.at(0).first.forward.right;
  int tempED = 0;
  for(auto e : combinedED){ // print all in order
    auto a = e.first.forward.left;
    auto b = e.first.forward.right;
    auto c = e.first.reverse.left;
    auto d = e.first.reverse.right;
    auto ed = e.second;
    if(ed > averageED){
      Interval_pair ip = Interval_pair(b+1,ra+1,d-1,rb-1);
      ra = a;
      rb = c;
      
      int length = (ip.forward.size() > ip.reverse.size())? ip.forward.size() : ip.reverse.size();
      cout << ip.toString() << "ed: " << tempED << " total lenght: " << length << endl;
      tempED = 0;
    }else{
      tempED += ed;
    }
  }
}
  

