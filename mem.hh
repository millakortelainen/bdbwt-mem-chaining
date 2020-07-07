#ifndef BDBWT_MEM_CHAIN_MEM_HH
#define BDBWT_MEM_CHAIN_MEM_HH
#include <iostream>
#include "include/BD_BWT_index.hh"
#include "include/Iterators.hh"
#include <string>
#include <set>
#include <stdio.h>
#include <inttypes.h>
#include <tuple>
#include <stack>
#include <bits/stdc++.h>
#include <future>
#include "util.hh"
#include "rsa1d.hh"
#include "rsa2d.hh"
#include <omp.h>

using namespace std;
  
bool verboseEnumeration = false;
bool verboseSigma = false;
bool verboseI = false;
bool verboseMaximum = false;
bool verboseSubroutine = false;
bool verboseElement = false;
bool verboseChaining = false;
bool naiveCrossProduct = false;

int minimumDepth = 1;

/** Setting up for, and processing MEM's via BDBWT computation. 
Takes the two BWT indices as an input, and optionally a vector containing interval tuple of seeds.
The seeds are primarily meant for the hybrid computation mode, but other use for it could also be found.
*/
vector<tuple<int,int,int>> bwt_to_int_tuples(BD_BWT_index<> index, BD_BWT_index<> index2, set<tuple<Interval_pair,Interval_pair,int>> seeds);
/** Main driver for finding MEM-matches with bi-directional BWT index.
*/
vector<tuple<int,int,int>> bwt_mem2(BD_BWT_index<> idxS, BD_BWT_index<> idxT, uint8_t startLabel, set<tuple<Interval_pair,Interval_pair,int>> seed);
/** Subroutine called from bwt_mem2, extends the given interval to the left, and right for all valid combinations gained from cross_product.
Returns location tuples for positions where a MEM match can be confirmed to be in.
*/
vector<tuple<int,int,int>> bwt_mem2_subroutine(BD_BWT_index<> idxS, BD_BWT_index<> idxT, pair<Interval_pair,Interval_pair> pr, int depth);

/** Translating bwt-mem tuples into tuples corresponding to the locations in the original texts.
 */
vector<tuple<int,int,int>> batchOutput(BD_BWT_index<> index, BD_BWT_index<> index2, vector<tuple<int,int,int>> memVector, bool verbose);
/** Subroutine of batchOuput.
 */
vector<struct occStruct> batchLocate(vector<struct occStruct>  pairs, vector<bool> marked, BD_BWT_index<> bwt);

/** Chaining function implements the chaining algorithm as defined in the paper "Chaining with Overlaps Revisited".
The algorithm chains together all MEM matches given as input, and finds the best chaining alignment such that the chains do not themselves, while allowing overlaps on the text
*/
vector<pair<int,pair<int,int>>> chaining(vector<Interval_pair> A, int size);

/** Output function for the chains.
 */
pair<vector<Interval_pair>,vector<int>> chainingOutput(vector<pair<int,pair<int,int>>> chains, vector<Interval_pair> Ipairs, string text, string text2);


/** Enumerates the unique characters on the forward index of the BWT.
    param idx BD_BWT_index<> Burrows-wheeler transform
    param ip Interval_pair range to find the unique characters in.
    return std::vector<uint8_t> array containing each unique character.
*/
vector<uint8_t> enumerateLeft(BD_BWT_index<> idx, Interval_pair ip);
/** Analogously to enumerateLeft 
*/
vector<uint8_t> enumerateRight(BD_BWT_index<> idx, Interval_pair ip);

/** More efficient implementation for cross-product 
*/
vector<tuple<uint8_t,uint8_t,uint8_t,uint8_t>> cross_product(vector<pair<uint8_t,uint8_t>> A, vector<pair<uint8_t,uint8_t>> B);


void naiveOutput(BD_BWT_index<> index, BD_BWT_index<> index2, vector<tuple<int,int,int>> memVector, string text1, string text2,bool verbose);









/** Enumerates the unique characters on the forward index of the BWT.
    param idx BD_BWT_index<> Burrows-wheeler transform
    param ip Interval_pair range to find the unique characters in.
    return std::vector<uint8_t> array containing each unique character.
*/
vector<uint8_t> enumerateLeft(BD_BWT_index<> idx, Interval_pair ip){
  vector<uint8_t> ret;
  if(ip.forward.left > idx.size() || ip.forward.right > idx.size()){
    return ret;
  }
  set<uint8_t> s;
  for(int i = ip.forward.left; i <= ip.forward.right; i++){
    auto c = idx.forward_bwt_at(i);
    s.insert(c); //Set is always sorted
  }
  if(verboseEnumeration) cout << "Enumerated left on forward interval of " << ip.toString() << "\n";
  ret.insert(ret.end(), s.begin(), s.end());
  return ret;
}
/** Enumerates the unique characters on the backward index of the BWT.
    param idx BD_BWT_index<> Burrows-wheeler transform
    param ip Interval_pair range to find the unique characters in.
    return std::vector<uint8_t> array containing each unique character.
*/
vector<uint8_t> enumerateRight(BD_BWT_index<> idx, Interval_pair ip){
  vector<uint8_t> ret;
  if(ip.reverse.left > idx.size() || ip.reverse.right > idx.size()){
    return ret;
  }
  set<uint8_t> s;
  for(int i = ip.reverse.left; i <= ip.reverse.right; i++){
    auto c = idx.backward_bwt_at(i);
    s.insert(c); //Set is always sorted
  }
  if(verboseEnumeration) cout << "Enumerated right on backward interval of " << ip.toString() << "\n";
  ret.insert(ret.end(), s.begin(), s.end());
  return ret;
}
vector<tuple<uint8_t,uint8_t,uint8_t,uint8_t>> cross_product(vector<pair<uint8_t,uint8_t>> A, vector<pair<uint8_t,uint8_t>> B){
  set<tuple<uint8_t,uint8_t,uint8_t,uint8_t>> bl;
  vector<tuple<uint8_t,uint8_t,uint8_t,uint8_t>> retbl;
  vector<pair<uint8_t,uint8_t>> aa;
  vector<pair<uint8_t,uint8_t>> ab;
  auto end = BD_BWT_index<>::END;
  while(A.size() > 0){
    auto a1 = A[0];
    auto a2 = A[0]; int a2i = 0;
    bool allIncompatible = true;
    for(int i = 1; i < A.size(); i++){
      if(a1.first != A[i].first && a1.second != A[i].second){
	a2 = A[i];
	a2i = i;
	break;
      }
    }
    
    for(auto bi : B){
      if((bi.first != a1.first && bi.second != a1.second)){
	allIncompatible = false;
	bl.insert(make_tuple(a1.first, a1.second, bi.first, bi.second));
      }
      if((bi.first != a2.first && bi.second != a2.second)){
	allIncompatible = false;
	bl.insert(make_tuple(a2.first, a2.second, bi.first, bi.second));
      }
      if((a1.first == end && bi.first == end) || (a1.second == end && bi.second == end)){ //Handle END marker
      	allIncompatible = false;
      	bl.insert(make_tuple(a1.first, a1.second, bi.first, bi.second));
      }
      if((a2.first == end && bi.first == end) || (a2.second == end && bi.second == end)){ //Handle END marker
      	allIncompatible = false;
      	bl.insert(make_tuple(a2.first, a2.second, bi.first, bi.second));
      }
    }
    if(allIncompatible){
      break;
    }
    A.erase(A.begin());
    if(a2i != 0){
      A.erase(A.begin()+a2i);
    }
  }
  
  for(int i = 0; i < A.size(); i++){
    aa.push_back(make_pair(A[i].first, '1'));
    ab.push_back(make_pair('2', A[i].second));

    for(auto bi : B){
      for(auto al : aa){
	if(bi.first != A[i].first && bi.second != A[i].second && al.second != bi.second){
	  bl.insert(make_tuple(al.first, al.second, bi.first, bi.second));
	}
      }
      for(auto al : ab){
	if(bi.first != A[i].first && bi.second != A[i].second && al.first != bi.first){
	  bl.insert(make_tuple(al.first, al.second, bi.first, bi.second));
	}
      }
    }
  }
  for(auto r : bl){
    retbl.push_back(r);
  }
  return retbl;
}
/** Subroutine for printing out indexes containing MEM's as triples (i,j,d). Where i, and j refer to starting indexes on the two texts, and d the depth (lenght) of the string.
param idxS BD_BWT_index<> BWT index built on the first text (S).
param idxT BD_BWT_index<> BWT index built on the second text (T).
param pr pair<Interval_pair, Interval_pair> tuple of interval_pairs corresponding to both BWT indexes. pr.first refers to forward and backward indexes of the first index (IdxS) and pr.second to the latter analogously.
param depth int

return vector<tuple<int,int,int>> vector containing tuples (i,j,d) indicating MEM matches.
*/
vector<tuple<int,int,int>> bwt_mem2_subroutine(BD_BWT_index<> idxS, BD_BWT_index<> idxT, pair<Interval_pair,Interval_pair> pr, int depth){
  vector<tuple<int,int,int>> ret;
  vector<pair<uint8_t,uint8_t>> A;
  vector<pair<uint8_t,uint8_t>> B;

  if(verboseSubroutine) cout << "Subroutine interval: " <<pr.first.toString() << pr.second.toString() << "\n";
#pragma omp parallel sections
  {
#pragma omp section
    {
      if(pr.first.forward.right > idxT.size()-1 || pr.first.reverse.right > idxT.size()-1){
      }else{
	for(auto a : enumerateLeft(idxS,pr.first)){
	  auto sea = idxS.left_extend(pr.first,a);
	  if(sea.forward.right > idxS.size()-1 || sea.reverse.right > idxS.size()-1){
	    continue;
	  }
	  for(auto b : enumerateRight(idxS,sea)){
	    A.push_back(make_pair(a,b));
	    if(verboseSubroutine) cout << "A: " << a << ", " << b << "\n";
	  }
	}
      }
    }
#pragma omp section
    {
      if(pr.second.forward.right > idxT.size()-1 || pr.second.reverse.right > idxT.size()-1){
      }else{
	for(auto c : enumerateLeft(idxT,pr.second)){
	  auto sea = idxT.left_extend(pr.second,c);
	  if(sea.forward.right > idxT.size()-1 || sea.reverse.right > idxT.size()-1){
	    continue;
	  }
	  for(auto d : enumerateRight(idxT,sea)){
	    B.push_back(make_pair(c,d));
	    if(verboseSubroutine) cout << "B: " << c << ", " << d << "\n";
	  }
	}
      }
    }
    
  }
  std::vector<tuple<uint8_t,uint8_t,uint8_t,uint8_t>> cross;
  uint8_t a,b,c,d;
  if(naiveCrossProduct){ //Set naiveCrossProduct to true to enable naive cross-product calculation (mainly for debugging).
    for(auto i : A){
      a = i.first;
      b = i.second;
      for(auto j: B){
	c = j.first;
	d = j.second;
	if((a != c && b != d) ||
	   (a == BD_BWT_index<>::END && c == BD_BWT_index<>::END && (b != d)) ||
	   (b == BD_BWT_index<>::END && d == BD_BWT_index<>::END && (a != b))) { 
	  cross.push_back(make_tuple(a,b,c,d));
	}
      }
      if(verboseSubroutine) cout << "---" << "\n";
    }
  }else{
    cross = cross_product(A,B);
  }
  for(int k = 0; k < cross.size(); k++){
    tie(a,b,c,d) = cross[k];
    Interval_pair i_1 = idxS.left_extend(pr.first,a);
    if(verboseSubroutine) cout << "Extended interval " << pr.first.toString() << " to the left with " << a << " and got: " << i_1.toString() << "\n";
    
    Interval_pair i_2 = idxS.right_extend(i_1,b);
    if(verboseSubroutine) cout << "Extended interval " << i_1.toString() << " to the right with " << b << " and got: " << i_2.toString() << "\n";
    
    Interval_pair i_3 = idxT.left_extend(pr.second,c);
    if(verboseSubroutine) cout << "Extended interval " << pr.second.toString() << " to the left with " << c << " and got: " << i_3.toString() << "\n";
    
    Interval_pair i_4 = idxT.right_extend(i_3,d);
    if(verboseSubroutine) cout << "Extended interval " << i_3.toString() << " to the right with " << d << " and got: " << i_4.toString() << "\n";

    for(int i = i_2.forward.left; i <= i_2.forward.right; i++){
      for(int j = i_4.forward.left; j <= i_4.forward.right; j++){
	ret.push_back(make_tuple(i,j,depth));
	if(verboseSubroutine) cout <<"[" <<omp_get_thread_num()<<"]" << "pushed result of: (" << i << ", " << j << ", " << depth << ")" << "from interval" << i_2.toString() << "\n";
      }
    }
  }
  return ret;
}

/** Find and return all MEM's between two BWT indexes.
    param idxS BD_BWT_index<> first BWT index
    param idxT BD_BWT_index<> second BWT index.

    return vector<tuple<int,int,int>> returns tuples (i,j,d) where i,j correspond to starting points of MEM's in the two indexes.
*/
vector<tuple<int,int,int>> bwt_mem2(BD_BWT_index<> idxS, BD_BWT_index<> idxT, uint8_t startLabel = BD_BWT_index<>::END, tuple<Interval_pair,Interval_pair,int> seed = make_tuple(Interval_pair(-1,-2,-1,-2),Interval_pair(-1,-2,-1,-2),-1)){
  vector<pair<pair<Interval_pair,Interval_pair>,int>> collectedSubroutineCalls;
  vector<tuple<int,int,int>> ret;
  set<tuple<Interval_pair,Interval_pair,int>> S;
  set<tuple<Interval_pair,Interval_pair,int>> processed;
  Interval_pair ip0, ip1; int depth = -1;
  int itrl1Size = idxS.size()-1;
  int itrl2Size = idxT.size()-1;
  int maxDepth = -1;
  
  if(get<0>(seed).forward.left <= 0){
    S.insert(make_tuple(Interval_pair(0,itrl1Size,0,itrl1Size),Interval_pair(0,itrl2Size,0,itrl2Size),0));
  }else{
    S.insert(seed);
  }
  while(!S.empty()){
    tie(ip0,ip1,depth) = *S.begin();
    
    if(depth > maxDepth && depth % 5 == 0){
      maxDepth = depth;
    }
    S.erase(S.begin());
    
    if((ip0.forward.right - ip0.forward.left+1) < 1 || (ip1.forward.right - ip1.forward.left+1) < 1){
      continue;
    }

    if(idxS.is_left_maximal(ip0) || idxT.is_left_maximal(ip1) ||
       (enumerateLeft(idxS, ip0) !=  enumerateLeft(idxT, ip1))  ||
       (enumerateLeft(idxS, ip0).size()==1 && enumerateLeft(idxS, ip0)[0] == BD_BWT_index<>::END)){ //Handle END symbols
      if(depth >= minimumDepth){
	collectedSubroutineCalls.push_back(make_pair(make_pair(ip0,ip1),depth));
      }
    }
    set<pair<Interval_pair,Interval_pair>> I;
    if(startLabel != BD_BWT_index<>::END && depth <= 0){
      Interval_pair i1 = idxS.left_extend(ip0,startLabel);
      Interval_pair i2 = idxT.left_extend(ip1,startLabel);
      if(i1.forward.left < 0 || i2.forward.left < 0){ //left_extend() returns [-1,-2] interval if extending with character c is not possible.
	continue; //no need to bother with invalid index
      }
      I.insert(make_pair(i1,i2));
    }else{
      auto Sigma = enumerateLeft(idxS,ip0);
      for(auto c : Sigma){
	if(c == BD_BWT_index<>::END){
	  continue;
	}
	Interval_pair i1 = idxS.left_extend(ip0,c);
	Interval_pair i2 = idxT.left_extend(ip1,c);
	if(i1.forward.left < 0 || i2.forward.left < 0){ //left_extend() returns [-1,-2] interval if extending with character c is not possible.
	  continue; //no need to bother with invalid index
	}
	I.insert(make_pair(i1,i2));
      }
    }
    pair<Interval_pair,Interval_pair> x;
    
    if(I.size() != 0){
      x = *I.begin(); //Initialize first item for comparison
    }
    else{
      continue;
    }
    int xForwardDelta 	= x.first.forward.right - x.first.forward.left;
    int xForwardDelta2 	= x.second.forward.right - x.second.forward.left;
    int maxDelta 	= xForwardDelta + xForwardDelta2;

    for(auto y : I){

      int yForwardDelta = y.first.forward.right - y.first.forward.left;
      int yForwardDelta2= y.second.forward.right - y.second.forward.left;
      int zDelta 	= yForwardDelta + yForwardDelta2;
      
      if(zDelta > maxDelta){
	x = y;
	maxDelta = zDelta;

      }
    }
    ip0 = x.first;
    ip1 = x.second;
    if(I.size() == 0){
      if(verboseI) cout << "\t" << "No values in I" << "\n";
    }
    else{
      I.erase(x);

      if(I.size() == 0){
      }
      else{
	for(auto y : I){
	  if(y.first.forward.right > idxS.size()-1 || y.second.forward.right > idxT.size()-1){
	    continue;
	  }
	  if(y.first.reverse.right > idxS.size()-1 || y.second.reverse.right > idxT.size()-1){
	    continue;
	  }
	  if(idxS.is_right_maximal(y.first) || idxT.is_right_maximal(y.second) ||
	     (enumerateRight(idxS, y.first) != enumerateRight(idxT, y.second)) ||
	     (enumerateRight(idxS, y.first).size()==1 && enumerateRight(idxS,y.first)[0] == BD_BWT_index<>::END)){  //Handle END symbols
	    
	    S.insert(make_tuple(y.first,y.second, depth+1));
	  }
	}
      }
    }
    //Have to take into the special case where we it is impossible to extend in one direction, but other direction might still have valid extensions left. If left side doesn't have any valid extensions, it will get filtered out on the next iteration before pushing anything into the stack. Only checking for enumerateLeft() could result in case where left side only has one possible extending character, and would not be reliable here.
    //std::vector<uint8_t> e1 = enumerateRight(idxS, x.first);
    // cout << "b" << endl;
    if(x.first.forward.right > idxS.size()-1 || x.second.forward.right > idxT.size()-1){
      continue;
    }
    if(x.first.reverse.right > idxS.size()-1 || x.second.reverse.right > idxT.size()-1){
      continue;
    }
    if(idxS.is_right_maximal(x.first) || idxT.is_right_maximal(x.second) ||
       (enumerateRight(idxS, x.first) != enumerateRight(idxT, x.second)) ||
       (enumerateRight(idxS, x.first).size()==1 && enumerateRight(idxS, x.first)[0] == BD_BWT_index<>::END)){  //Handle END symbols
      S.insert(make_tuple(x.first,x.second, depth+1));
    }
  }

  if(collectedSubroutineCalls.size() == 0){
    //cout << "Could not find any MEM's with significiant enough length" << endl;
    return ret;
  }
  //  cout << "collected subroutine (" <<collectedSubroutineCalls.size()-1<<") elements with min depth of: " << minimumDepth << endl;
  // int dmin = INT_MAX;
  // int dmax = INT_MIN;
  // int sum = 0;
  // for(auto s : collectedSubroutineCalls){
  //   if(s.second < dmin){
  //     dmin = s.second;
  //   }
  //   if(s.second > dmax){
  //     dmax = s.second;
  //   }
  //   sum += s.second;
  // }
  // int dfilter = ((sum/collectedSubroutineCalls.size())-1);
  // dfilter = (dfilter > minimumDepth)? dfilter: minimumDepth;
  // dfilter = minimumDepth;
  // cout << "minimum depth based on sampling " << dfilter << endl;
  // ret.reserve(collectedSubroutineCalls.size());
  vector<vector<tuple<int,int,int>>> rettempThreadContainer(omp_get_max_threads());
  
#pragma omp parallel for
  for(int i = 0; i < collectedSubroutineCalls.size(); i++){
    auto p = collectedSubroutineCalls[i];
    if(p.second < minimumDepth){
      continue;
    }
    if(verboseSubroutine) cout << "Enter subroutine with depth: " << p.second << "\n";
    auto rettemp = bwt_mem2_subroutine(idxS,idxT,p.first,p.second);
    rettempThreadContainer[omp_get_thread_num()].insert(rettempThreadContainer[omp_get_thread_num()].end(), rettemp.begin(), rettemp.end());
    for(auto rt : rettemp){
      rettempThreadContainer[omp_get_thread_num()].push_back(rt);
    }
  }
  for(auto i : rettempThreadContainer){
    ret.insert(ret.end(), i.begin(), i.end());
  }
  return ret;
}
	       
vector<struct occStruct> batchLocate(vector<struct occStruct>  pairs, vector<bool> marked, BD_BWT_index<> bwt){
  int i = 0;
  int n = bwt.size();
  vector<struct occStruct> translate;
  vector<struct occStruct> ret;
  auto LFindex = mapLF(bwt, true).first;
  
  for(int j = n; j > 0; j--){
    if(marked[i]){
      struct occStruct temp;      
      temp.first = i;
      temp.second = j-1;
      translate.push_back(temp);	
    }
    i = LFindex[i];
  }
  
  #pragma omp parallel sections
  {
    #pragma omp section
    {
      translate = radixSort(translate, 2);
    }
    #pragma omp section
    {
      pairs = radixSort(pairs,2);
    }
  }

  int x = 0;
  int y = 0;
  while(x < pairs.size()){
    auto a = pairs[x];
    if(a.first == translate[y].first){
      struct occStruct temp;
      temp.first = translate[y].second;
      temp.second = a.second;
      ret.push_back(temp);
      x++;
    }else{
      if(y == translate.size()-1){
	y = 0;
      }else{
	y++;
      }
    }
  }
  return ret;
}

vector<pair<int,pair<int,int>>> chaining(vector<Interval_pair> A, int size){
  cout << "beginning chaining for " << A.size() << " intervals...";
  rsa1d T_a = rsa1d(size); 
  rsa1d T_b = rsa1d(size);
  rsa2d T_c = rsa2d(size);
  rsa2d T_d = rsa2d(size);
  vector<pair<int, int>> C_a(A.size(), make_pair(0,0));
  vector<pair<int, int>> C_b(A.size(), make_pair(0,0));
  vector<pair<int, int>> C_c(A.size(), make_pair(0,0));
  vector<pair<int, int>> C_d(A.size(), make_pair(0,0));
  vector<pair<int, pair<int,int>>> C_p(A.size(), make_pair(0,make_pair(0,0)));
  vector<pair<int, int>> C(A.size(), make_pair(0,0));
  vector<pair<int, int>> E_1;
  vector<pair<int, int>> E_2;
  T_a.insertCell(make_pair(0,0),INT_MIN);
  T_b.insertCell(make_pair(0,0),INT_MIN);
  for(int j = 0; j < A.size(); j++){
    T_a.insertCell(make_pair(A[j].reverse.right,j), INT_MIN);
    T_b.insertCell(make_pair(A[j].reverse.right,j), INT_MIN);

    auto key = make_pair(A.at(j).reverse.left - A.at(j).forward.left,j);
    auto skey1 = A.at(j).forward.right;
    auto skey2 = A.at(j).reverse.right;
    
    if(verboseChaining) cout << "Initializer updating : " << key.first << ", " << key.second << ", with secondary keys: " << skey1 << " & " << skey2 << endl;
    T_c.insertCell(key, skey1, INT_MIN);
    T_d.insertCell(key, skey2, INT_MIN);
    if(verboseChaining) cout << endl;

    auto p1 = make_pair(A[j].forward.left, j);
    auto p2 = make_pair(A[j].forward.right, j);
    E_1.push_back(p1);
    E_1.push_back(p2);
    if(verboseChaining) cout << "Pushed E_1: " << p1.first << ", " << p1.second << " & " << p2.first << ", " << p2.second << endl; 
  }
  T_a.sortArray();
  T_b.sortArray();
  T_c.sortArray();
  T_d.sortArray();
  
  if(verboseChaining) cout << "E_1 size: " << E_1.size() << endl;
  T_a.upgrade(make_pair(0,0),0);
  sort(E_1.begin(), E_1.end(), pairSort);

  auto a = T_a.rangeMax(-1,-1); //Initializing variable types.
  auto b = T_b.rangeMax(-1,-1);
  auto c = T_c.rangeMax(-1,-1,-1,-1);
  auto d = T_d.rangeMax(-1,-1,-1,-1);
      
  cout << "chaining initializer done" << endl;
  for(int i = 0; i < E_1.size(); i++){
    auto e = E_1[i];
    if(verboseChaining) cout << E_1[i].first << endl;
    int j = e.second;
      
    Interval_pair I = A[j];
    auto rlfl = I.reverse.left - I.forward.left;
    
    if(verboseChaining) cout << A[j].toString() << endl;
    if(verboseChaining) cout << "i: " << i << ", j: " << j << endl;
    if(verboseChaining) cout << "I.forward.left == E_1[i].first == " << I.forward.left << " == " << E_1[i].first << endl;
    
    if(I.forward.left == E_1[i].first){
      if(verboseChaining) cout <<"Interval: " << I.toString() << endl; 
      
#pragma omp parallel sections
      {
#pragma omp section
	{
	  a = T_a.rangeMax(0,			I.reverse.left); //I.reverse.left-1 causes failure when (reverse.left = forward.left = 0) with zeroth-index-indexing
	}
#pragma omp section
	{
	  b = T_b.rangeMax(I.reverse.left,	I.reverse.right);
	}
#pragma omp section
	{
	  c = T_c.rangeMax(INT_MIN,		rlfl,	0, I.forward.right);
	}
#pragma omp section
	{
	  d = T_d.rangeMax(rlfl+1,		INT_MAX,0, I.reverse.right);
	}
      }
      
      C_a[j].first = a.first.primary.second;	C_a[j].second = a.second; 
      C_b[j].first = b.first.primary.second;	C_b[j].second = I.reverse.left + b.second; 
      C_c[j].first = c.first.primary.second;	C_c[j].second = I.forward.left + c.second;
      C_d[j].first = d.first.primary.second;	C_d[j].second = I.reverse.left + d.second;
      
      if(verboseChaining) cout << "maxCandinates: "
			       << C_a[j].second << "("<< C_a[j].first << "),"
			       << C_b[j].second << "("<< C_b[j].first << "),"
			       << C_c[j].second << "("<< C_c[j].first << "),"
			       << C_d[j].second << "("<< C_d[j].first << ")," << endl;
      
      auto max = chainingMax(C_a[j].second, C_b[j].second,C_c[j].second,C_d[j].second);
      if     (C_a[j].second == max) C[j] = C_a[j];
      else if(C_b[j].second == max) C[j] = C_b[j];
      else if(C_c[j].second == max) C[j] = C_c[j];
      else if(C_d[j].second == max) C[j] = C_d[j];

      auto cpsum = C[j].second+I.forward.right-I.forward.left+1;
      C_p[j] = make_pair(cpsum, make_pair(C[j].first,j));
      
      if(verboseChaining) cout << "max C[j] = " << C[j].second << endl;
      if(verboseChaining) cout << "C_p= " << cpsum << ", " << j << endl;
      if(verboseChaining) cout << "I.forward.right = " << I.forward.right << " I forward.left = " << I.forward.left << " C_p = " << cpsum << endl;

      auto upgsumL = (int)C[j].second-(int)I.forward.left;
      T_c.upgrade(make_pair(rlfl, j), I.forward.right, upgsumL);
      if(verboseChaining) cout << "upgraded T_c (" << I.reverse.left - I.forward.left << ", " << I.forward.right << "), j= " << j << A[j].toString() <<" into: " << upgsumL << endl;

      auto upgsumR = (int)C[j].second-(int)I.reverse.left;
      T_d.upgrade(make_pair(rlfl, j), I.forward.right, upgsumR);
      if(verboseChaining) cout << "upgraded T_d (" << I.reverse.left - I.forward.left << ", " << I.forward.right << "), j= " << j << A[j].toString() <<" into: " << upgsumR << endl;  
    }
    else{
      if(verboseChaining) cout << "Else statement, upgrade and update!" << endl;
      T_a.upgrade(make_pair(I.reverse.right, j), C_p[j].first);
      T_b.upgrade(make_pair(I.reverse.right, j), C[j].second- I.reverse.left);
      T_c.update(make_pair(rlfl, j), I.forward.right, INT_MIN);
      T_d.update(make_pair(rlfl, j), I.forward.right, INT_MIN);
    }
    if(verboseChaining) cout << "end \n";
  }
  if(verboseChaining) cout << "end2 \n";
  return C_p;
}

void naiveOutput(BD_BWT_index<> index, BD_BWT_index<> index2, vector<tuple<int,int,int>> memVector, string text1 = "", string text2 = "",bool verbose = false){
  auto retSA = buildSAfromBWT(index, true); //RetSA builds SA array for the given text from it's BWT transform without having to use the extra space from permutating whole original text.
  auto retSA2 = buildSAfromBWT(index2, true);

  //Naive returning of the intervals
  for(auto a : memVector){
    int i, j, depth;
    tie(i,j,depth) = a;
    int begin_i= retSA[i].first+1;
    int begin_j= retSA2[j].first+1;
    
    //Handling the specific special case when SA^S[i] = text.size(); We use index.size()-1 instead so wouldn't need to keep the text in memory at all.
    if(begin_i >= index.size()-1){
      begin_i = index.size()-retSA[i].first-1; //Index.size() will always be text.size()+1 due to the added END marker. Furthermore, we need to minus one to get proper offset from general case of begin_i = retSA[i]+1;
    }
    //Handling the specific special case when SA^T[j] = text2.size(); We use index2.size()-1 instead so wouldn't need to keep the text in memory at all.
    if(begin_j >= index2.size()-1){
      begin_j = index2.size()-retSA2[j].first-1; //Analogously to above. 
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
    
    newOcc.first  = i;
    newOcc2.first = j;
    newOcc.second  = p;
    newOcc2.second = p;
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
    tie(i,j,d) = memVector[bl1[k].second];
    int ik = bl1[k].first+1;
    int jk = bl2[k].first+1;
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
pair<vector<Interval_pair>,vector<int>> chainingOutput(vector<pair<int,pair<int,int>>> chains, vector<Interval_pair> Ipairs, string text, string text2){
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
    auto a = c.forward;
    auto b = c.reverse;
    int d = c.forward.right - c.forward.left+1;
    cout << "Chain["<< count <<"]: "<< c.toString() <<"\t\t symcov:" << symcov.at(count) << endl;
    //    cout << text.substr(a.left,d) <<",\t "<< text2.substr(b.left,d) << endl;
    count++;
  }
  return (make_pair(chainIntervals, symcov));
}
// Giving seeds an default value so we can easily leave it out if we don't use them.
vector<tuple<int,int,int>> bwt_to_int_tuples(BD_BWT_index<> index, BD_BWT_index<> index2, set<tuple<Interval_pair,Interval_pair,int>> seeds = {make_tuple(Interval_pair(-1,-2,-1,-2),Interval_pair(-1,-2,-1,-2),-1)}){
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
    if(seeds.size() > 1){
      vector<tuple<Interval_pair,Interval_pair,int>> seedsVector;
      copy(seeds.begin(), seeds.end(), back_inserter(seedsVector));
#pragma omp parallel for
      for(int i = 0; i < seedsVector.size(); i++){
	auto retMem = bwt_mem2(index, index2, BD_BWT_index<>::END, seedsVector[i]);
	memThreads[omp_get_thread_num()].insert(memThreads[omp_get_thread_num()].end(), retMem.begin(), retMem.end());
      }
    }
    else{
#pragma omp parallel for
    for(int i = 0; i < enumLeft.size(); i++){
      auto retMem = bwt_mem2(index, index2, enumLeft.at(i));
      memThreads[omp_get_thread_num()].insert(memThreads[omp_get_thread_num()].end(), retMem.begin(), retMem.end());
    }
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

#endif
