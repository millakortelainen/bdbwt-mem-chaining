#include <stdlib.h>
#include <iostream>
#include "mem.hh"
#include "edlib.h"
#include <stdio.h>
#include "io.hh"
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
    int end_i  = begin_i+depth-1;
    int end_j  = begin_j+depth-1;
    
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
pair<vector<Interval_pair> ,vector<int>> chainingOutput(vector<pair<int,pair<int,int>>> chains, vector<Interval_pair> Ipairs){
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
    //cout << "chains[i].first= " << chains.at(i).first << "...\t\t";
    //cout << "chains[i].second= " << chains.at(i).second.first <<":"<<chains.at(i).second.second << "...\t";
    if(chains.at(i).second.second >= 0){
      auto I = Ipairs.at(chains.at(i).second.second);
      if(chainIntervals.size() > 0 &&
	 chainIntervals.at(chainIntervals.size()-1).forward.left >= I.forward.left &&
	 chainIntervals.at(chainIntervals.size()-1).reverse.left >= I.reverse.left){ //Ensuring (weak) precedence

	symcov.push_back(chains.at(i).first);
	chainIntervals.push_back(Ipairs.at(chains.at(i).second.second));
	last = i;
	i = chains[i].second.first;
	if(last == 0 || last == i){ //would print out same index again => chain is done.
	  break;
	}
      }
    }
  }
  int count = 0;
  for(auto c : chainIntervals){
    cout << "Chain["<<count<<"]: "<<c.toString() <<"\t\t symcov:" << symcov.at(count) << endl;
    count++;
  }
  return (make_pair(chainIntervals, symcov));
}
  
int main(int argc, char *argv[]){
  string text;
  string text2;
  auto fileinput = readInputFromFasta(argv[1]);
  switch(0){
  case 0: {
    text = fileinput.at(0);
    text2 = fileinput.at(1);
    break;
  }
  case 1: {
    text  = "GTGCGTGATCATCATTT";
    text2 = "AGTGCAAAGTGATTACC";
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

  EdlibAlignResult result = edlibAlign(text.c_str(), text.size()-1, text2.c_str(), text2.size()-1, edlibDefaultAlignConfig());
  printf("edit_distance('hello', 'world!') = %d\n", result.editDistance);
  edlibFreeAlignResult(result);

  auto deptharg = argv[2];
  minimumDepth = strtol(argv[2],NULL,10);  
  BD_BWT_index<> index((uint8_t*)text.c_str());
  BD_BWT_index<> index2((uint8_t*)text2.c_str());
  
  auto mems  = bwt_mem2(index, index2);//Find MEMS between two BDBWT indexes.
  sort(mems.begin(), mems.end(), memSort); //Proper sorting of the tuples with priority order of i --> d --> j
  auto filtered = filterMems(mems);
  
  //pretty_print_all(index,text);
  //pretty_print_all(index2,text2);

  //  naiveOutput(index,index2,filtered,text,text2, true);
  //cout << endl;
  auto bo = batchOutput(index, index2, filtered, false);
  cout << "batchOutput done" << endl;
  sort(bo.begin(), bo.end(), memSort); //overall Speed increase
  
  vector<Interval_pair> Ipairs = returnMemTuplesToIntervals(bo, false);
  Ipairs = filterIntervals(Ipairs);
 
  // for(int i = 0; i < Ipairs2.size(); i++){
  //   cout <<"Ipairs["<< i << "]: " << Ipairs2[i].toString() << "\n";
  // }
  
  auto chains = chaining(Ipairs, text2.size());
  cout << "Chaining done" << endl;

  auto chainints = chainingOutput(chains, Ipairs).first;

  auto absent = absentIntervals(chainints, index, index2);
  for(int i = 0; i < absent.size(); i++){
    cout <<"absent["<< i << "]: " << absent[i].toString() << "\n";
  }
}
  

