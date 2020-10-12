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
#include "util.hh"
#include <omp.h>
using namespace std;
bool verboseEnumeration = false;
bool verboseSigma = false;
bool verboseI = false;
bool verboseMaximum = false;
bool verboseSubroutine = false;
bool verboseElement = false;
bool naiveCrossProduct = false;

int minimumDepth = 1;

/** Setting up for, and processing MEM's via BDBWT computation. 
Takes the two BWT indices as an input, and optionally a vector containing interval tuple of seeds.
The seeds are primarily meant for the hybrid computation mode, but other use for it could also be found.
*/
vector<tuple<int,int,int>> bwt_to_int_tuples(Configuration conf, set<tuple<Interval_pair,Interval_pair,int>> seeds);
/** Main driver for finding MEM-matches with bi-directional BWT index.
*/
vector<tuple<int,int,int>> bwt_mem2(Configuration conf, uint8_t startLabel, set<tuple<Interval_pair,Interval_pair,int>> seed);
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


/** Checks for difference between the enumerations for the two indexes and intervals
*/
bool enum_diff(BD_BWT_index<> idxS, BD_BWT_index<> idxT, Interval_pair ip0, Interval_pair ip1, bool dir){ 
  if(dir){
    if(ip0.forward.left == ip0.forward.right == 1 || ip1.forward.left == ip1.forward.right == 1){
      if(idxS.forward_bwt_at(ip0.forward.left) == BD_BWT_index<>::END){
        return true;
      }
      if(idxT.forward_bwt_at(ip1.forward.left) == BD_BWT_index<>::END){
        return true;
      }
    }
    unordered_set<uint8_t> s;
    for(int i = ip0.forward.left; i <= ip0.forward.right; i++){
      auto c = idxS.forward_bwt_at(i);
      s.insert(c);
    }
    for(int j = ip1.forward.left; j <= ip1.forward.right; j++){
      auto c = idxT.forward_bwt_at(j);
      if(s.count(c) == 0){
        return true;
      }
    }
  }else{
    if(ip0.reverse.left == ip0.reverse.right == 1 || ip1.reverse.left == ip1.reverse.right == 1){
      if(idxS.backward_bwt_at(ip0.reverse.left) == BD_BWT_index<>::END){
        return true;
      }
      if(idxT.backward_bwt_at(ip1.reverse.left) == BD_BWT_index<>::END){
        return true;
      }
    }
    unordered_set<uint8_t> s;
    for(int i = ip0.reverse.left; i <= ip0.reverse.right; i++){
      auto c = idxS.backward_bwt_at(i);
      s.insert(c);
    }
    for(int j = ip1.reverse.left; j <= ip1.reverse.right; j++){
      auto c = idxT.backward_bwt_at(j);
      if(s.count(c) == 0){
        return true;
      }
    }
  }
  return false;
}


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

/** Computes Cartesian product between A,(a,b) and B,(c,d), such that a =! c and b != d.
 */
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
      else if((bi.first != a2.first && bi.second != a2.second)){
        allIncompatible = false;
        bl.insert(make_tuple(a2.first, a2.second, bi.first, bi.second));
      }
      else if((a1.first == end && bi.first == end) || (a1.second == end && bi.second == end)){ //Handle END marker
      	allIncompatible = false;
      	bl.insert(make_tuple(a1.first, a1.second, bi.first, bi.second));
      }
      else if((a2.first == end && bi.first == end) || (a2.second == end && bi.second == end)){ //Handle END marker
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
        ret.push_back(make_tuple(i,j,depth-1));
        if(verboseSubroutine) cout <<"[" <<omp_get_thread_num()<<"]" << "pushed result of: (" << i << ", " << j << ", " << depth << ")" << "from interval" << i_2.toString() << "\n";
      }
    }
  }
  //cout << "returning " << ret.size() << " mems from sub" << endl;
  return ret;
}

/** Find and return all MEM's between two BWT indexes.
    param idxS BD_BWT_index<> first BWT index
    param idxT BD_BWT_index<> second BWT index.

    return vector<tuple<int,int,int>> returns tuples (i,j,d) where i,j correspond to starting points of MEM's in the two indexes.
*/
vector<tuple<int,int,int>> bwt_mem2(Configuration conf, uint8_t startLabel = BD_BWT_index<>::END, tuple<Interval_pair,Interval_pair,int> seed = make_tuple(Interval_pair(-1,-2,-1,-2),Interval_pair(-1,-2,-1,-2),-1)){
  auto idxS = conf.index1;
  auto idxT = conf.index2;
  vector<pair<pair<Interval_pair,Interval_pair>,int>> collectedSubroutineCalls;
  vector<tuple<int,int,int>> ret;
  set<tuple<int,int,int>, memTupleSortStruct> retSet;
  vector<tuple<Interval_pair,Interval_pair,int>> S;
  set<tuple<Interval_pair,Interval_pair,int>> processed;
  Interval_pair ip0, ip1; int depth = -1;
  int itrl1Size = idxS.size()-1;
  int itrl2Size = idxT.size()-1;
  int maxDepth = -1;
  bool seeded = false;

  if(get<0>(seed).forward.left <= 0){
  S.push_back(make_tuple(Interval_pair(0,itrl1Size,0,itrl1Size),Interval_pair(0,itrl2Size,0,itrl2Size),0));
  }else{
    S.push_back(seed);
    seeded = true;
  }
  while(!S.empty()){
    tie(ip0,ip1,depth) = *S.begin();
    S.erase(S.begin());
    if((ip0.forward.right - ip0.forward.left+1) < 1 || (ip1.forward.right - ip1.forward.left+1) < 1){
       continue;
    }
    bool diff = enum_diff(idxS,idxT,ip0,ip1,true);
    if(diff && !seeded) { //Handle END symbols
      if(depth >= minimumDepth && !seeded){
        collectedSubroutineCalls.push_back(make_pair(make_pair(ip0,ip1),depth));
      }
    }
    if(seeded){
      //bool diff2 = enum_diff(idxS,idxT,ip0,ip1,false);
      if(diff){
        auto rettemp = (bwt_mem2_subroutine(idxS,idxT,make_pair(ip0,ip1),depth));
        for(auto t : rettemp){ 
          ret.push_back(t);
        }
      }
    }
    set<pair<Interval_pair,Interval_pair>> I;
    set<pair<Interval_pair,Interval_pair>> Ir;
    if(startLabel != BD_BWT_index<>::END && depth <= 0){
      Interval_pair i1 = idxS.left_extend(ip0,startLabel);
      Interval_pair i2 = idxT.left_extend(ip1,startLabel);
      if(i1.forward.left < 0 || i2.forward.left < 0){ //left_extend() returns [-1,-2] interval if extending with character c is not possible.
        continue; //no need to bother with invalid index
      }
      I.insert(make_pair(i1,i2));
    }else{
      auto Sigma = conf.alphabet;
      for(auto c : Sigma){
        if(c == BD_BWT_index<>::END){
          //cout << "could not extend " << ip0.toString() << "," << ip1.toString() << " to left with " << c << endl;
          continue;
        }
        Interval_pair i1 = idxS.left_extend(ip0,c);
        Interval_pair i2 = idxT.left_extend(ip1,c);
        if(i1.forward.left < 0 || i2.forward.left < 0){ //left_extend() returns [-1,-2] interval if extending with character c is not possible.
          //cout << "could not extend " << ip0.toString() << "," << ip1.toString() << " to left with " << c << endl;
          continue; //no need to bother with invalid index
        }
        I.insert(make_pair(i1,i2));
      }
    }
      // pair<Interval_pair,Interval_pair> x;
      
      //   if(I.size() != 0){
      //     x = *I.begin(); //Initialize first item for comparison
  //   }
  //   else{
  //     continue;
  //   }
  //   int xForwardDelta 	= x.first.forward.right - x.first.forward.left;
  //   int xForwardDelta2 	= x.second.forward.right - x.second.forward.left;
  //   int maxDelta 	= xForwardDelta + xForwardDelta2;

  //   for(auto y : I){
  //     int yForwardDelta = y.first.forward.right - y.first.forward.left;
  //     int yForwardDelta2= y.second.forward.right - y.second.forward.left;
  //     int zDelta 	= yForwardDelta + yForwardDelta2;
      
  //     if(zDelta > maxDelta){
	// x = y;
	// maxDelta = zDelta;

  //     }
  //   }
    //ip0 = x.first;
    // ip1 = x.second;
    if(I.size() == 0){
      //if(true) cout << "\t" << "No values in I, going right" << "\n";



    }
    else{
      //  I.erase(x);
        
      // if(I.size() == 0){
      // }
      // else{
      // cout << "coming form ->" << ip0.toString() << ip1.toString() << endl;
      for(auto y : I){
        //cout << "value in I: " << y.first.toString() << y.second.toString() << endl;
        // if(y.first.forward.right > idxS.size()-1 || y.second.forward.right > idxT.size()-1){
        //   continue;
        // }
        // if(y.first.reverse.right > idxS.size()-1 || y.second.reverse.right > idxT.size()-1){
        //   continue;
        // }
        if(idxS.is_right_maximal(y.first) || idxT.is_right_maximal(y.second) ||
           enum_diff(idxS,idxT,y.first,y.second,false)
           //(enum_diff(idxS,idxT,y.first,y.second,true) && seeded)
           ){  //Handle END symbols
          S.push_back(make_tuple(y.first,y.second, depth+1));
          // continue;
        }
      }
    }
    if(seeded){
      for(auto c : conf.alphabet){
        if(c == BD_BWT_index<>::END){
          continue;
        }
        Interval_pair i1 = idxS.right_extend(ip0,c);
        Interval_pair i2 = idxT.right_extend(ip1,c);
        if(i1.forward.left < 0 || i2.forward.left < 0){
          continue;
        }
        Ir.insert(make_pair(i1,i2));
        //cout << "ir insert" << endl;
      }
      for(auto y : Ir){
        if(enum_diff(idxS,idxT,y.first,y.second,true)
           //  || enumerateLeft(idxS,y.first).back() != BD_BWT_index<>::END
           //|| enumerateLeft(idxT,y.second).back() != BD_BWT_index<>::END
           || !idxS.is_left_maximal(y.first)
           || !idxT.is_left_maximal(y.second)
           // || true
           ){  //Handle END symbols
          S.push_back(make_tuple(y.first,y.second, depth+1));
          //cout << "ir s insert" << endl;
        }
      }
      }
    }
    //Have to take into the special case where we it is impossible to extend in one direction, but other direction might still have valid extensions left. If left side doesn't have any valid extensions, it will get filtered out on the next iteration before pushing anything into the stack. Only checking for enumerateLeft() could result in case where left side only has one possible extending character, and would not be reliable here.
    //std::vector<uint8_t> e1 = enumerateRight(idxS, x.first);
    // cout << "b" << endl;
    // if(x.first.forward.right > idxS.size()-1 || x.second.forward.right > idxT.size()-1){
    //   continue;
    // }
    // if(x.first.reverse.right > idxS.size()-1 || x.second.reverse.right > idxT.size()-1){
    //   continue;
    // }
    // if(idxS.is_right_maximal(x.first) || idxT.is_right_maximal(x.second) ||
    //    (enumerateRight(idxS, x.first) != enumerateRight(idxT, x.second)) ||
    //    (enumerateRight(idxS, x.first).size()==1 && enumerateRight(idxS, x.first)[0] == BD_BWT_index<>::END)){  //Handle END symbols
    //   S.insert(make_tuple(x.first,x.second, depth+1));
    // }
  //}

  if(collectedSubroutineCalls.size() == 0 || seeded){
    return ret;
  }
  if(conf.verbosity > 2) cout << "collected subroutine (" <<collectedSubroutineCalls.size()<<") elements with min depth of: " << minimumDepth << endl;
  vector<vector<tuple<int,int,int>>> rettempThreadContainer(omp_get_max_threads());
  int totalcount = 0;
#pragma omp parallel for
  for(int i = 0; i < collectedSubroutineCalls.size(); i++){
    auto p = collectedSubroutineCalls[i];
    if(p.second < minimumDepth){
      continue;
    }
    if(verboseSubroutine) cout << "Enter subroutine with depth: " << p.second << "\n";
    auto rettemp = bwt_mem2_subroutine(idxS,idxT,p.first,p.second);
    int count = 0;
    for(auto rt : rettemp){
      rettempThreadContainer[omp_get_thread_num()].push_back(rt);
      count++;
      totalcount++;
    }
  }
  if(conf.verbosity > 2) cout<< "total count " << totalcount << endl;
  for(auto i : rettempThreadContainer){
    for(auto j : i){
      retSet.insert(retSet.end(), j);
    }
  }
  ret.insert(ret.end(), retSet.begin(),retSet.end());
  return ret;
}

/** Subroutine to batchOutput: Convert BDBWT MEM tuples into text-interval tuples
 */
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
/** Conversion from BDBWT MEM tuples into text-interval tuples.
 */
vector<tuple<int,int,int>> batchOutput(BD_BWT_index<> index, BD_BWT_index<> index2, vector<tuple<int,int,int>> memVector, bool verbose = false){
  std::vector<struct occStruct> Ipairs;
  std::vector<struct occStruct> Ipairs2;
  std::vector<bool> marked1(index.size(),false);
  std::vector<bool> marked2(index2.size(),false);
  vector<struct occStruct> bl1;
  vector<struct occStruct> bl2;
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
#pragma omp parallel sections
  {
#pragma omp section
    {
      bl1 = batchLocate(Ipairs,marked1,index);
    }
#pragma omp section
    {
      bl2 = batchLocate(Ipairs2,marked2,index2);
    }
  }
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
/** Giving seeds an default value so we can easily leave it out if we don't use them.
 */
vector<tuple<int,int,int>> bwt_to_int_tuples(Configuration conf, set<tuple<Interval_pair,Interval_pair,int>> seeds = {make_tuple(Interval_pair(-1,-2,-1,-2),Interval_pair(-1,-2,-1,-2),-1)}){
  auto index = conf.index1;
  auto index2 = conf.index2;
  minimumDepth = conf.minimumDepth;
  vector<tuple<int,int,int>> mems;
  if(conf.threadedBWT){
    vector<tuple<int,int,int>> memFilter;
    if(conf.verbosity > 2) cout << "finding mems between indexes...";
    auto enumLeft = enumerateLeft(index,Interval_pair(0, index.size()-1, 0, index.size()-1));
    if(enumLeft.at(0) == BD_BWT_index<>::END){
      enumLeft.erase(enumLeft.begin());
    }
    vector<vector<tuple<int,int,int>>> memThreads(omp_get_max_threads());
    if(seeds.size() > 1){
      vector<tuple<Interval_pair,Interval_pair,int>> seedsVector;
      copy(seeds.begin(), seeds.end(), back_inserter(seedsVector));
#pragma omp parallel for
      for(int i = 0; i < seedsVector.size(); i++){
        auto retMem = bwt_mem2(conf, BD_BWT_index<>::END, seedsVector[i]);
        memThreads[omp_get_thread_num()].insert(memThreads[omp_get_thread_num()].end(), retMem.begin(), retMem.end());
      }
    }
    else{
#pragma omp parallel for
      for(int i = 0; i < enumLeft.size(); i++){
        auto retMem = bwt_mem2(conf, enumLeft.at(i));
        memThreads[omp_get_thread_num()].insert(memThreads[omp_get_thread_num()].end(), retMem.begin(), retMem.end());
      }
    }
    int count = 0;
#pragma omp single
    for(auto a : memThreads){
      for(auto b : a){
        memFilter.push_back(b);
      }
    }
    sort(memFilter.begin(),memFilter.end(), memSort);
    set<tuple<int,int,int>, memTupleSortStruct> uniqueTuples;
    for(auto a : memFilter){
      uniqueTuples.insert(uniqueTuples.end(), a);
    }
    mems.reserve(uniqueTuples.size()+1);
    for(auto a : uniqueTuples){
      mems.emplace_back(a);
    }
  }
  else{
    mems = bwt_mem2(conf, BD_BWT_index<>::END);
  }

  if(conf.verbosity > 2) cout << "found mems," << mems.size() << "...";
  sort(mems.begin(), mems.end(), memSort); //Proper sorting of the tuples with priority order of i --> d --> j

  if(mems.size() == 0){
    cout << "could not find significiantly large enough MEMS " << endl;
    return mems;
  }
  auto bo = batchOutput(index, index2, mems, false);
  if(conf.verbosity > 2) cout << "batchOutput into SA indices done" << "...";
  sort(bo.begin(), bo.end(), memSort); //overall Speed increase
  return bo;
}
#endif
