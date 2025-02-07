#ifndef BDBWT_MEM_CHAIN_UTIL_HH
#define BDBWT_MEM_CHAIN_UTIL_HH
#include <iostream>
#include <string>
#include "include/BD_BWT_index.hh"
#include <stdio.h>
#include <cmath>
#include "edlib.h"
using namespace std;

int memSort(tuple<int,int,int> set1, tuple<int,int,int> set2);
bool pairSort(const pair<int,int> &first, const pair<int,int> &second);
bool minimizerLexSort(const pair<string,int> &first, const pair<string,int> &second);
bool intervalSort(const Interval_pair &a, const Interval_pair &b);
vector<struct occStruct> radixSort(vector<struct occStruct> list, int r);

vector<Interval_pair> returnMemTuplesToIntervals(vector<tuple<int,int,int>> tup, bool sortIntervals);
vector<Interval_pair> absentIntervals(vector<Interval_pair> chains, BD_BWT_index<> idx1, BD_BWT_index<> idx2);
pair<map<int,int>,int> mapLF(BD_BWT_index<>& index, bool forward);
/* Naive function to build suffix array from input text 
 */
vector<pair<int,int>> buildSAfromBWT(BD_BWT_index<> idxS, bool revIndexSA);
/* Obsolete recursion function for building SA from BWT.
 */
int chainingMax(int &a, int &b, int &c, int &d);

set<char> alphabet;
/* BEGIN_STRUCTS */
struct occStruct{
  int first;
  int second;
};
struct Configuration {
  int mode;
  int verbosity = 1;
  std::string text1;
  std::string text2;
  BD_BWT_index<> index1;
  BD_BWT_index<> index2;
  set<char> alphabet;
  
  int minimumDepth;
  int minimizerWindowSize;
  int miniMergerCount = 0;
  bool threadedBWT = true;

  int PLCPSparsity_q = 1;
  int maxSize = (index1.size() > index2.size())? index1.size() : index2.size();
  EdlibAlignConfig edlibConf;
  int originalEditDistance = -1;

  bool printAbsentAndChains = false;
  bool verboseEditDistances = false;
  bool rawChains = false;
  bool chainStringSegments = false;
  bool recombAbsents = false;
  bool useLinearRangeMax = false;

};

struct cell1d{
  pair<int,int> primary;
  int value;
};
static int cellcmp1d(struct cell1d first, struct cell1d last){
  return (first.primary.first < last.primary.first)? 1:0;
}

struct memTupleSortStruct {
  bool operator() (const tuple<int,int,int>& set1, const tuple<int,int,int>& set2) const{
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
        if(j != y ){
          return (j < y);
        }
      }
    }
    return 0;
}
};

/* END_STRUCTS */

/* BEGIN_SORT */
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
bool newChainingSort(const Interval_pair &a, const Interval_pair &b){
  return (a.forward.left < b.forward.left);
}
bool intervalIntPairSort(const pair<Interval_pair,int> &c, const pair<Interval_pair,int> &d){
  auto a = c.first;
  auto b = d.first;
  int aac = (a.reverse.left - a.forward.left);
  int bac = (b.reverse.left - a.forward.left);
  int maxa = (a.forward.right > a.reverse.right)? a.forward.right : a.reverse.right;
  int maxb = (b.forward.right > b.reverse.right)? b.forward.right : b.reverse.right;
  return (a < b);
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
    Interval_pair temp(i,i+d, j,j+d);
    Ipairs.push_back(temp);
  }
  if(sortIntervals){
    sort(Ipairs.begin(), Ipairs.end(), intervalSort);
  }
  return Ipairs;
}

vector<Interval_pair> absentIntervalsV2(vector<Interval_pair> chains, BD_BWT_index<>idx1, BD_BWT_index<> idx2){
  vector<Interval> abF;
  vector<Interval> abR;
  vector<Interval_pair> abRet;
  int i = 0;
  int j = 0;
  int fwdEnd = idx1.size()-1;
  int bwdEnd = idx2.size()-1;
  Interval temp = Interval(-1, -1);

  while(i < chains.size()){
    temp = Interval(-1,-1);
    auto a = chains[i].forward.left;
    auto b = chains[i].forward.right;
    if(b > fwdEnd){
      abF.push_back(temp);
      fwdEnd = a-1;
    }else{
    temp.right = fwdEnd;
    temp.left  = b+1;
    fwdEnd = a-1;
    abF.push_back(temp);
    }

    temp = Interval(-1,-1);
    auto c = chains[i].reverse.left;
    auto d = chains[i].reverse.right;
    if(d > bwdEnd){
      abR.push_back(temp);
      bwdEnd = c-1;
    }else{
    temp.right = bwdEnd;
    temp.left  = d+1;
    bwdEnd = c-1;
    abR.push_back(temp);
  }
    i++;
  }
  cout << "abF size: " << abF.size() << "abR size: " << abR.size() << endl;
  cout << "abF end left: " << abF.back().left << "abR end left: " << abR.back().left << endl;
  auto temp2 = Interval_pair(-1,-1,-1,-1);
  if(abF.back().left != 0){
    temp2.forward.right = fwdEnd;
    temp2.forward.left = 0;
    abF.push_back(temp2.forward);
    abR.push_back(temp2.reverse);
  }else if(abR.back().left != 0){
    temp2.reverse.right = bwdEnd;
    temp2.reverse.left = 0;
    abF.push_back(temp2.forward);
    abR.push_back(temp2.reverse);
  }
  for(int i = 0; i < abF.size(); i++){
    //if(abF[i].right > idx1.size()-1 || abR[i].right > idx2.size()-1) continue;
    abRet.push_back(Interval_pair(abF[i], abR[i]));
    cout << "Chain["<<i<<"]: " << chains[i].toString() << endl;
    cout << "Absent["<<i<<"]: " << abRet.back().toString() << endl;
  }
  return abRet;
}
vector<Interval_pair> absentIntervals(vector<Interval_pair> chains, BD_BWT_index<> idx1, BD_BWT_index<> idx2){
  if(false){
    return absentIntervalsV2(chains,idx1,idx2);
  }
  vector<Interval_pair> absents;
  //vector<pair<Interval_pair,string>> both;
  int fwdLastEnd = idx1.size()-1;
  int bwdLastEnd = idx2.size()-1;
  Interval_pair temp = Interval_pair(-1, -1, -1, -1);

  if(chains[0].forward.right != fwdLastEnd || chains[0].reverse.right != bwdLastEnd){ // Either side of first(last) chain does not match the end of the input text.
    if(chains[0].forward.right != fwdLastEnd){
      temp.forward.left = chains[0].forward.right+1;
      temp.forward.right = fwdLastEnd;
      fwdLastEnd = chains[0].forward.left-1;
    }
    if(chains[0].reverse.right != bwdLastEnd){
      temp.reverse.left = chains[0].reverse.right+1;
      temp.reverse.right = bwdLastEnd;
      bwdLastEnd = chains[0].reverse.left-1;
    }
    absents.push_back(temp);
    //both.push_back(make_pair(temp, "Absent"));
    //both.push_back(make_pair(chains[0],"Chain"));
  }
  for(int i = 1; i < chains.size(); i++){
    //  int j = 1;
    if(i >= 1){
      temp = Interval_pair(-1,-1,-1,-1);
    }
    auto ch = chains[i];
    auto a = ch.forward.left;
    auto b = ch.forward.right;
    auto c = ch.reverse.left;
    auto d = ch.reverse.right;
  //  cout << i << "FWDLastEnd: "<< fwdLastEnd << " BWDLastEnd: " << bwdLastEnd << "Current Chain: " << ch.toString() << endl;
    if(a == 0){
//	    cout << "case one" << endl;
      temp.forward.left = fwdLastEnd;
      temp.forward.right = b-1;
    }
    if(b != fwdLastEnd){
//	    cout << "case two" << endl;
      temp.forward.right = fwdLastEnd;
      temp.forward.left = b+1;
    }
    if(c == 0){
//	    cout << "case three" << endl;
      temp.reverse.left = bwdLastEnd;
      temp.reverse.right = d-1;
    }
    if(d != bwdLastEnd){
//	    cout << "case four" << endl;
      temp.reverse.right = bwdLastEnd;
      temp.reverse.left = d+1;
    }
    if(temp.forward.right < temp.forward.left){
	    continue;
//	    cout << "case five" << endl;
      int t = temp.forward.right;
      temp.forward.right = temp.forward.left;
      temp.forward.left = t;
    }
    if(temp.reverse.right < temp.reverse.left){
	    continue;
//	    cout << "case six" << endl;
      int t = temp.reverse.right;
      temp.reverse.right = temp.reverse.left;
      temp.reverse.left = t;
    }
    fwdLastEnd = a-1;
    bwdLastEnd = c-1;
    absents.push_back(temp);
   // cout << endl;
   // both.push_back(make_pair(temp,"Absent"));
   // both.push_back(make_pair(chains[i], "Chain"));
	
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
    //both.push_back(make_pair(temp,"Absent"));
  }
  int lastEndF = idx1.size();
  /*for(int a = 0; a < both.size(); a++){
	  auto b = both[a];
	  cout << b.second <<"\t"<< b.first.toString() << "\t" << lastEndF-b.first.forward.right-1 << endl;
	  if(b.first.forward.left != -1){
	  lastEndF = b.first.forward.left;
	  }
	  a++;
	  if(a < both.size()){
	  b = both[a];
	  cout << b.second <<"\t"<< b.first.toString() << "\t" << lastEndF-b.first.forward.right-1 << endl << endl;
	  if(b.first.forward.left != -1){
	  lastEndF = b.first.forward.left;
	  }
	  }
  }
  cout << "return" << endl;*/
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
/** Creating SA with through BWT.
*/
vector<pair<int,int>> buildSAfromBWT(BD_BWT_index<> idxS, bool revIndexSA){
  //  cout << idxS.size()-1 << " <- size()" << endl;
  vector<int> retSA (idxS.size(),-1);
  vector<int> retISA(idxS.size(),-1);

  auto lfS = mapLF(idxS, revIndexSA);
  //vector<int64_t> C = idxS.get_global_c_array();
  // vector<int> SA(idxS.size(),-1);
  // vector<int> ISA(idxS.size(),-1);
  // SA[lfS.second] = idxS.size()-1;
  // ISA[idxS.size()-1] = lfS.second;
  int currIndex = 0;
  int k = 0;
  //cout << "currIndex: " << currIndex << endl;
  while(k < idxS.size()){
    //cout << "currIndex= " << currIndex << endl;
    //cout<<idxS.forward_bwt_at(currIndex)<<endl;
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

/* END_AUX */
#endif
 
