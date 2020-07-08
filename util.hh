#ifndef BDBWT_MEM_CHAIN_UTIL_HH
#define BDBWT_MEM_CHAIN_UTIL_HH
#include <iostream>
#include <string>
#include "include/BD_BWT_index.hh"
#include <stdio.h>
#include <cmath>
#include "rsa1d.hh"
#include "edlib.h"
using namespace std;

int suffcmp(struct suffix first, struct suffix last);
int memSort(tuple<int,int,int> set1, tuple<int,int,int> set2);
bool pairSort(const pair<int,int> &first, const pair<int,int> &second);
bool pairSortSecond(const pair<int,int> &first, const pair<int,int> &second);
bool chainPairSort(const pair<int,pair<int,int>> &first, const pair<int,pair<int,int>> &second);
bool minimizerLexSort(const pair<string,int> &first, const pair<string,int> &second);
bool intervalSort(const Interval_pair &a, const Interval_pair &b);
bool intervalSortDiff(const Interval_pair &a, const Interval_pair &b);
vector<struct occStruct> radixSort(vector<struct occStruct> list, int r);

vector<Interval_pair> returnMemTuplesToIntervals(vector<tuple<int,int,int>> tup, bool sortIntervals);
vector<Interval_pair> filterIntervals(vector<Interval_pair> Ipairs, int size);
vector<Interval_pair> absentIntervals(vector<Interval_pair> chains, BD_BWT_index<> idx1, BD_BWT_index<> idx2);
pair<map<int,int>,int> mapLF(BD_BWT_index<>& index, bool forward);
/* Naive function to build suffix array from input text 
 */
int* build_suffix_array(string text);
/* Building suffix array in space efficient manner directly from BWT.
 */
vector<pair<int,int>> buildSAfromBWT(BD_BWT_index<> idxS, bool revIndexSA);
/* Obsolete recursion function for building SA from BWT.
 */
pair< vector<int>,vector<int> > int_ret_recurse(BD_BWT_index<> idxS, map<int,int> LFI, int currIndex, int k, vector<int> retSA, vector<int> retISA, int i);

int chainingMax(int &a, int &b, int &c, int &d);
vector<tuple<int,int,int>> filterMems(vector<tuple<int,int,int>> mems);

/* Pretty print for certain informations. Mainly intended as debugging tool.
 */
void pretty_print_all(BD_BWT_index<>& index, string text1);
  



set<char> alphabet;
/* BEGIN_STRUCTS */
struct suffix{
  int index;
  string suffix;
};
struct occStruct{
  int first;
  int second;
};
struct Configuration {
  int mode;
  std::string text1;
  std::string text2;
  BD_BWT_index<> index1;
  BD_BWT_index<> index2;

  int minimumDepth;
  int minimizerWindowSize;
  int PLCPSparsity_q = 1;
  int maxSize = (index1.size() > index2.size())? index1.size() : index2.size();
  EdlibAlignConfig edlibConf;
  int originalEditDistance = -1;
  
};
/* END_STRUCTS */

/* BEGIN_SORT */
int suffcmp(struct suffix first, struct suffix last){
  return (first.suffix < last.suffix)? 1:0;
}
int memSort(tuple<int,int,int> set1, tuple<int,int,int> set2){
  int i,j,d;
  int x,y,z;
  tie(i,j,d) = set1;
  tie(x,y,z) = set2;
  
  if(i != x){
    return (i < x);
  }
  else{
    if(d != z){
      return (d > z);
    }
    else{
      if(j != y){
	return (j < y);
      }
    }
  }
  return 0;
}

bool pairSort(const pair<int,int> &first, const pair<int,int> &second){
  return (first.first < second.first);
}
bool pairSortSecond(const pair<int,int> &first, const pair<int,int> &second){
  return (first.second < second.second);
}
bool chainPairSort(const pair<int,pair<int,int>> &first, const pair<int,pair<int,int>> &second){
  return (first.first < second.first);
}
bool minimizerLexSort(const pair<string,int> &first, const pair<string,int> &second){
  return (first.first < second.first);
}
bool intervalSort(const Interval_pair &a, const Interval_pair &b){
  int aac = (a.reverse.left - a.forward.left);
  int bac = (b.reverse.left - a.forward.left);
  int maxa = (a.forward.right > a.reverse.right)? a.forward.right : a.reverse.right;
  int maxb = (b.forward.right > b.reverse.right)? b.forward.right : b.reverse.right;
  return (a < b);
}
bool intervalSortDiff(const Interval_pair &a, const Interval_pair &b){
  int aac = (a.reverse.left - a.forward.left);
  int bac = (b.reverse.left - b.forward.left);
  int maxa = (a.forward.right > a.reverse.right)? a.forward.right : a.reverse.right;
  int maxb = (b.forward.right > b.reverse.right)? b.forward.right : b.reverse.right;
  return (abs(aac) < abs(bac));
}

vector<struct occStruct> radixSort(vector<struct occStruct> list, int r){
  int rad = 10; 
  std::list<struct occStruct> radix[rad];
  for(int i = 0; i < r; i++){
    int r1 = pow(rad,i+1);
    int r2 = pow(rad,i);
    for(int j = 0; j < list.size(); j++){
      auto idx = ((list[j].second) % r1) / r2;
      radix[idx].push_back(list[j]);
    }
	
    int k = 0;
    for(int j = 0; j < rad; j++){
      while(!radix[j].empty()){
	list[k] = *(radix[j].begin());
	radix[j].erase(radix[j].begin());
	k++;
      }
    } 
  }
  return list;
}
/* END_SORT */

/*BEGIN_INTERVALS*/
vector<Interval_pair> returnMemTuplesToIntervals(vector<tuple<int,int,int>> tup, bool sortIntervals){
  vector<Interval_pair> Ipairs;
  for(auto b : tup){
    int i,j,d;
    tie(i,j,d) = b;
    Interval_pair temp(i,i+d-1, j,j+d-1);
    Ipairs.push_back(temp);	 
  }
  if(sortIntervals){
    sort(Ipairs.begin(), Ipairs.end(), intervalSort);
  }
  return Ipairs;
}

vector<Interval_pair> filterIntervals(vector<Interval_pair> Ipairs, int size){
  //  const int s = size;
  vector<int> bsl(size,-1);
  vector<int> bsr(size,-1);
  vector<Interval_pair> Ipairs2;
  sort(Ipairs.begin(), Ipairs.end(), intervalSortDiff);
  Ipairs2.push_back(Ipairs[0]);

  for(int i = 1; i < Ipairs.size(); i++){
    int count = 0;
    auto ip = Ipairs.at(i);
    auto io = Ipairs2.at(Ipairs2.size()-1);
    for(int j = ip.forward.left; j <= ip.forward.right; j++){
      if(bsl.at(j) > -1){
	count++;
      }
    }
    for(int j = ip.reverse.left; j <= ip.reverse.right; j++){
      if(bsr.at(j) > -1){
	count++;
      }
    }
    if(count > ((abs(ip.reverse.left - ip.forward.left)))){
      //continue;
    }
    if(ip.forward.left >= io.forward.left && ip.forward.right <= io.forward.right){ //nested case
      if(abs(ip.reverse.left-ip.forward.left) < abs(io.reverse.left-io.forward.left)){
	Ipairs2.pop_back();	
      }
      Ipairs2.push_back(ip);
      for(int j = ip.forward.left; j <= ip.forward.right; j++){
	bsl.at(j) = i;
      }
      for(int j = ip.reverse.left; j <= ip.reverse.right; j++){
	bsr.at(j) = i;
      }
    }else{
      Ipairs2.push_back(ip);
      for(int j = ip.forward.left; j <= ip.forward.right; j++){
	bsl.at(j) = i;
      }
      for(int j = ip.reverse.left; j <= ip.reverse.right; j++){
	bsr.at(j) = i;
      }
    }
  }
  return Ipairs2;
}
vector<Interval_pair> absentIntervals(vector<Interval_pair> chains, BD_BWT_index<> idx1, BD_BWT_index<> idx2){

  vector<Interval_pair> absents;
  int fwdLastEnd = idx1.size()-1;
  int bwdLastEnd = idx1.size()-1;
  
  if(chains[0].forward.right != fwdLastEnd || chains[0].reverse.right != bwdLastEnd){ // Either side of first(last) chain does not match the end of the input text.
    Interval_pair temp = Interval_pair(-1, -1, -1, -1);
    if(chains[0].forward.right != fwdLastEnd){
      temp.forward.left = chains[0].forward.right+1;
      temp.forward.right = fwdLastEnd;
      fwdLastEnd = chains[0].forward.left-1;
    }
    if(!chains[0].reverse.right != bwdLastEnd){
      temp.reverse.left = chains[0].reverse.right+1;
      temp.reverse.right = bwdLastEnd;
      bwdLastEnd = chains[0].reverse.left-1;
    }
    absents.push_back(temp);
  }
    
  for(int i = 1; i < chains.size(); i++){
    Interval_pair temp = Interval_pair(-1,-1,-1,-1);
    auto ch = chains[i];
    auto a = ch.forward.left;
    auto b = ch.forward.right;
    auto c = ch.reverse.left;
    auto d = ch.reverse.right;

    if(a == 0){
      temp.forward.left = fwdLastEnd;
      temp.forward.right = b-1;
    }
    if(b != fwdLastEnd){
      temp.forward.right = fwdLastEnd;
      temp.forward.left = b+1;
    }
    if(d != bwdLastEnd){
      temp.reverse.right = bwdLastEnd;
      temp.reverse.left = d+1;
    }
    if(c == 0){
      temp.reverse.left = bwdLastEnd;
      temp.reverse.right = d-1;
    }
    if(temp.forward.right < temp.forward.left){
      int t = temp.forward.right;
      temp.forward.right = temp.forward.left;
      temp.forward.left = t;
    }
    if(temp.reverse.right < temp.reverse.left){
      int t = temp.reverse.right;
      temp.reverse.right = temp.reverse.left;
      temp.reverse.left = t;
    }
    fwdLastEnd = a-1;
    bwdLastEnd = c-1;
    absents.push_back(temp);
  }
  if(fwdLastEnd > 0 || bwdLastEnd > 0){ //Either side of the first(last) chain does not match the beginning of the input text.
    Interval_pair temp = Interval_pair(-1, -1, -1, -1);
    if(fwdLastEnd != -1){
      temp.forward.left = 0;
      temp.forward.right = fwdLastEnd;
    }
    if(bwdLastEnd != -1){
      temp.reverse.left = 0;
      temp.reverse.right = bwdLastEnd;
    }
    absents.push_back(temp);
  }

  return(absents);
}
  
/*END_INTERVALS*/

/*BEGIN_AUX*/
/** Map LF values for the given index.
    return pair<map<int,int>,int> such that .first is the actual mapping, and .second is the zeroth index so that we don't need to search for it again.
*/
pair<map<int,int>,int> mapLF(BD_BWT_index<>& index, bool forward){
  map<char,int> m;
  map<int,int> lfmapping;
  auto C = index.get_global_c_array();
  int zeroth_index = -1;
  for(auto i : alphabet){ 
    m[i] = 0;
  }
  
  for(int i = 0; i < index.size(); i++){
    char l;
    if(forward){
      l = char(index.forward_bwt_at(i));
    }
    else{
      l = char(index.backward_bwt_at(i));
    }
    lfmapping[i] = C[l]+m[l];
    if(lfmapping[i] == 0){
      zeroth_index = i;
    }
    m[l] += 1;      
  }
  return make_pair(lfmapping,zeroth_index);
}

int* build_suffix_array(string text){
  text += BD_BWT_index<>::END;
  struct suffix suffix_array[text.size()+1];
  int *suffix_index_array = new int[text.size()+1];

  for(int i = 0; i < text.size(); i++){
    suffix_array[i].index = i;
    suffix_array[i].suffix = text.substr(i,text.size()+1);
  }
  sort(suffix_array,suffix_array+text.size(),suffcmp);

  for(int i = 0; i < text.size(); i++){
    cout << suffix_array[i].index << "\t" << suffix_array[i].suffix  << "\n";
    suffix_index_array[i] = suffix_array[i].index;
  }
  return suffix_index_array;
}

/** Creating SA with use of recursive LF mapping.
    The second in the pair is obsolete and should be removed. (from pretty print as well)
*/
pair< vector<int>,vector<int> > int_ret_recurse(BD_BWT_index<> idxS, map<int,int> LFI, int currIndex, int k, vector<int> retSA, vector<int> retISA, int i = 0){
  if(k < idxS.size()-1){
    cout << currIndex << "\t" << i << " ";
    cout << LFI.at(currIndex) << endl;
    //retSA[currIndex] = make_pair(idxS.size()-k-1, currIndex);
    retSA.at(currIndex) = idxS.size()-k-1;
    retISA.at(idxS.size()-k-1) = currIndex;
    return int_ret_recurse(idxS, LFI, LFI.at(currIndex), k+1, retSA, retISA, i+1);
  }else{
    retSA[currIndex] = idxS.size()-k-1;
    retISA[idxS.size()-k-1] = currIndex;
    return make_pair(retSA, retISA);
  }  
}
/** Creating SA with use of recursive LF mapping.
    Obsolete outside of naive output and pretty print
*/
vector<pair<int,int>> buildSAfromBWT(BD_BWT_index<> idxS, bool revIndexSA){
  //  cout << idxS.size()-1 << " <- size()" << endl;
  vector<int> retSA (idxS.size(),-1);
  vector<int> retISA(idxS.size(),-1);
  
  auto lfS = mapLF(idxS, revIndexSA);
  //vector<int64_t> C = idxS.get_global_c_array();
  vector<int> SA(idxS.size(),-1);
  vector<int> ISA(idxS.size(),-1);
  SA[lfS.second] = idxS.size()-1;
  ISA[idxS.size()-1] = lfS.second;
  int currIndex = lfS.first[lfS.second];
  int k = 0;
  //cout << "while" << endl;
  while(k < idxS.size()){
    retSA.at(currIndex) = idxS.size()-k-1;
    retISA.at(idxS.size()-k-1) = currIndex;
    k++;
    currIndex = lfS.first.at(currIndex);
  }
  //auto recursion = int_ret_recurse(idxS, lfS.first, lfS.first[lfS.second],  0, SA, ISA);
  //  retSA = recursion.first;
  //retISA = recursion.second;

  vector<pair<int,int>> ret;
  for(int i = 0; i < retSA.size(); i++){
    ret.push_back(make_pair(retSA[i], retISA[i]));
  }
  
  return ret;
}

int chainingMax(int &a, int &b, int &c, int &d){
  auto e = (a > b)? a:b;
  auto f = (c > d)? c:d;
  return (e > f)? e:f;
}

vector<tuple<int,int,int>> filterMems(vector<tuple<int,int,int>> mems){
  vector<tuple<int,int,int>> filtered;
  for(auto m : mems){
    int i,j,d;
    tie(i,j,d) = m;
    if(filtered.size() == 0){
      filtered.push_back(m);
    }else{
      int x,y,z;
      tie(x,y,z) = filtered[filtered.size()-1];
      if(i == x && y == j){
	
      }else{
	filtered.push_back(m);
      }
    }
  }
  return filtered;
}
/* END_AUX */



/** Print out BWT, LF and SA indices.
param index BD_BWT_index<>& Burrows Wheeler index.
param text1 string original text that the BWT is based on.
*/
void pretty_print_all(BD_BWT_index<>& index, string text1){
  string separator = "+------------------------------------------------+\n";
  cout << separator;
  auto SA = buildSAfromBWT(index, true);
  cout << separator;
  cout << "Text: " << text1 << "\n";
  cout << separator;
  cout << "| bwd\t" << "ln\t" << "fwd\t" << "LF\t" << "SA[i]\t" << "SA[i].s" << "\t |\n";
  auto mapping = mapLF(index, true).first;
  
  for(int i = 0; i < index.size(); i++){
    char t = (char)index.forward_bwt_at(i);
    char tr= (char)index.backward_bwt_at(i);
    cout << "| "<< tr << "\t("<<i<<")\t" << t << "\t" << mapping[i] << "\t" << SA[i].first << "\t" << SA[i].second << "\t |\n";
  }
  cout << separator;
}
#endif
