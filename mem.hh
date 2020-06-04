#include <iostream>
#include "include/BD_BWT_index.hh"
#include "include/Iterators.hh"
#include "util.hh"
#include "rsa1d.hh"
#include "rsa2d.hh"
#include <string>
#include <set>
#include <stdio.h>
#include <inttypes.h>
#include <tuple>
#include <stack>
#include <bits/stdc++.h>
#include <future>
int minimumDepth = 1;
using namespace std;

bool verboseEnumeration = false;
bool verboseSigma = false;
bool verboseI = false;
bool verboseMaximum = false;
bool verboseSubroutine = false;
bool verboseElement = false;
bool verboseChaining = false;
bool naiveCrossProduct = false;


/** Enumerates the unique characters on the forward index of the BWT.
    param idx BD_BWT_index<> Burrows-wheeler transform
    param ip Interval_pair range to find the unique characters in.
    return std::vector<uint8_t> array containing each unique character.
*/
std::vector<uint8_t> enumerateLeft(BD_BWT_index<> idx, Interval_pair ip){
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
  for(auto j : s){
    ret.push_back(j);
    if(verboseEnumeration){
      if(j == BD_BWT_index<>::END){
	cout << "$" << ", ";
      }else{
	cout << j <<", ";
      }
      cout << "\n";
    }
  }
  return ret;
}
/** Enumerates the unique characters on the backward index of the BWT.
    param idx BD_BWT_index<> Burrows-wheeler transform
    param ip Interval_pair range to find the unique characters in.
    return std::vector<uint8_t> array containing each unique character.
*/
std::vector<uint8_t> enumerateRight(BD_BWT_index<> idx, Interval_pair ip){
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
  for(auto j : s){
    ret.push_back(j);
    if(verboseEnumeration){
      if(j == BD_BWT_index<>::END){
	cout << "$" << ", ";
      }else{
	cout << j <<", ";
      }
      cout << "\n";
    }
  }
  return ret;
}
vector<tuple<uint8_t,uint8_t,uint8_t,uint8_t>> cross_product(vector<pair<uint8_t,uint8_t>> A, vector<pair<uint8_t,uint8_t>> B){
  vector<tuple<uint8_t,uint8_t,uint8_t,uint8_t>> bl;
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
	bl.push_back(make_tuple(a1.first, a1.second, bi.first, bi.second));
      }
      if((bi.first != a2.first && bi.second != a2.second)){
	allIncompatible = false;
	bl.push_back(make_tuple(a2.first, a2.second, bi.first, bi.second));
      }
      if((a1.first == end && bi.first == end) || (a1.second == end && bi.second == end)){ //Handle END marker
      	allIncompatible = false;
      	bl.push_back(make_tuple(a1.first, a1.second, bi.first, bi.second));
      }
      if((a2.first == end && bi.first == end) || (a2.second == end && bi.second == end)){ //Handle END marker
      	allIncompatible = false;
      	bl.push_back(make_tuple(a2.first, a2.second, bi.first, bi.second));
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
	  bl.push_back(make_tuple(al.first, al.second, bi.first, bi.second));
	}
      }
      for(auto al : ab){
	if(bi.first != A[i].first && bi.second != A[i].second && al.first != bi.first){
	  bl.push_back(make_tuple(al.first, al.second, bi.first, bi.second));
	}
      }
    }
  }
  return bl;
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
  for(auto a : enumerateLeft(idxS,pr.first)){
    auto sea = idxS.left_extend(pr.first,a);
    for(auto b : enumerateRight(idxS,sea)){
      A.push_back(make_pair(a,b));
      if(verboseSubroutine) cout << "A: " << a << ", " << b << "\n";
    }
  }
  for(auto c : enumerateLeft(idxT,pr.second)){
    auto sea = idxT.left_extend(pr.second,c);
    for(auto d : enumerateRight(idxT,sea)){
      B.push_back(make_pair(c,d));
      if(verboseSubroutine) cout << "B: " << c << ", " << d << "\n";
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
	if(verboseSubroutine) cout << "pushed result of: (" << i << ", " << j << ", " << depth << ")" << "\n";
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
vector<tuple<int,int,int>> bwt_mem2(BD_BWT_index<> idxS, BD_BWT_index<> idxT){
  vector<tuple<int,int,int>> ret;
  stack<tuple<Interval_pair,Interval_pair,int>> S;
  Interval_pair ip0, ip1; int depth = -1;
  int itrl1Size = idxS.size()-1;
  int itrl2Size = idxT.size()-1;
  
  S.push(make_tuple(Interval_pair(0,itrl1Size,0,itrl1Size), Interval_pair(0,itrl2Size,0,itrl2Size),0));
  while(!S.empty()){
    if(verboseElement) cout << "\n\n---Popping new element:" << "";
    tie(ip0,ip1,depth) = S.top(); S.pop();
    if(verboseElement) cout << ip0.toString() << ip1.toString() << ", current depth is: " << depth << " \nStack has: " << S.size() << " elements left."  << "\n\n";
    
    if((ip0.forward.right - ip0.forward.left+1) < 1 || (ip1.forward.right - ip1.forward.left+1) < 1){
      if(verboseElement) cout << "Invalid Element" << ip0.toString() << ip1.toString() << "\n";
      continue;
    }
    if(idxS.is_left_maximal(ip0) || idxT.is_left_maximal(ip1) ||
       (enumerateLeft(idxS,ip0)  != enumerateLeft(idxT,ip1))  ||
       (enumerateLeft(idxS,ip0).size()==1 && enumerateLeft(idxS,ip0)[0] == BD_BWT_index<>::END)){ //Handle END symbols
      if(depth >= minimumDepth){
	if(verboseSubroutine) cout << "Enter subroutine with depth: " << depth << "\n";
	auto rettemp = bwt_mem2_subroutine(idxS,idxT,make_pair(ip0,ip1),depth);
	for(auto i : rettemp){
	  ret.push_back(i);
	}
      }
    }
    
    if(verboseSigma) cout << "---\n";
    auto Sigma = enumerateLeft(idxS,ip0);
    set<pair<Interval_pair,Interval_pair>> I;
    for(auto c : Sigma){
      uint8_t cp = c;
      Interval_pair i1 = idxS.left_extend(ip0,c);
      Interval_pair i2 = idxT.left_extend(ip1,c);
      if(c == BD_BWT_index<>::END){
      	continue;
      }
      if(verboseSigma){
	if(cp == BD_BWT_index<>::END){
	  cp = '$';
	}
	cout << "Sigma: " << cp << "\n";
	cout << "\tExtending interval to left on IdxS" << ip0.toString() << " with character: " << cp << " ... -> " << i1.toString()<<"\n";
	cout << "\tExtending interval to left on IdxT" << ip1.toString() << " with character: " << cp << " ... -> " << i2.toString()<<"\n";
      }
      if(i1.forward.left < 0 || i2.forward.left < 0){ //left_extend() returns [-1,-2] interval if extending with character c is not possible.
	continue; //no need to bother with invalid index
      }
      I.insert(make_pair(i1,i2));
      if(verboseI) cout << "\t\tI size: " << I.size() << " Added: " << i1.toString() << i2.toString() << "\n";
    }
    if(verboseSigma) cout << "---\n";
    pair<Interval_pair,Interval_pair> x;
    
    if(I.size() != 0){
      x = *I.begin(); //Initialize first item for comparison
    }
    else{
      continue; // I is empty, nothing to do.
    }
    int xForwardDelta 	= x.first.forward.right - x.first.forward.left;
    int xForwardDelta2 	= x.second.forward.right - x.second.forward.left;
    int maxDelta 	= xForwardDelta + xForwardDelta2;
    if(verboseMaximum) cout << "\n";
    for(auto y : I){
      if(verboseMaximum) cout << "Comparing x:\t " << x.first.toString() << "\t" << x.second.toString() << "\nagainst y:\t " << y.first.toString() << "\t" << y.second.toString() << "\n";

      int yForwardDelta = y.first.forward.right - y.first.forward.left;
      int yForwardDelta2= y.second.forward.right - y.second.forward.left;
      int zDelta 	= yForwardDelta + yForwardDelta2;
      
      if(zDelta > maxDelta){
	x = y;
	maxDelta = zDelta;
	if(verboseMaximum) cout << "Found greater:\t " << y.first.toString() << "\t" << y.second.toString() << "\n";
      }
    }
    if(verboseMaximum) cout << "\nMaximum interval computed was: x = " << x.first.toString() << x.second.toString() << " Proceeding to remove from I..."<< "\n";
    ip0 = x.first;
    ip1 = x.second;
    
    if(I.size() == 0){
      if(verboseI) cout << "\t" << "No values in I" << "\n";
    }
    else{
      I.erase(x);
      if(verboseI) cout << "Removed x from I, I is now: " << "\n";
      if(I.size() == 0){
	if(verboseI) cout << "\t" << "No values in I" << "\n";
      }
      else{
	if(verboseI){
	  for(auto i : I){
	    cout << "\t" << i.first.toString() << i.second.toString() << "\n";	  
	  }
	}
	for(auto y : I){
	  if(verboseElement){
	    cout << "Is " << y.first.toString() << y.second.toString() << " right maximal on idxS? " << idxS.is_right_maximal(y.first) << ", on idxT? " << idxT.is_right_maximal(y.second) << "\n";
	    for(auto e : enumerateRight(idxS, y.first)){
	      cout << e << ",";
	    }
	    cout << "\n";
	    for(auto e : enumerateRight(idxT,y.second)){
	      cout << e << ",";
	    }
	    cout << "\n";															       
	  }	 
	  if(idxS.is_right_maximal(y.first) || idxT.is_right_maximal(y.second) ||
	     (enumerateRight(idxS, y.first) != enumerateRight(idxT, y.second)) ||
	     (enumerateRight(idxS, y.first).size()==1 && enumerateRight(idxS,y.first)[0] == BD_BWT_index<>::END)){  //Handle END symbols
	    
	    S.push(make_tuple(y.first,y.second, depth+1));
	    if(verboseElement) cout << "Pushed y:" << y.first.toString() << y.second.toString() << "\n";
	  }
	}
      }
    }
    if(verboseElement){
      cout << "Is x right maximal on idxS? " << idxS.is_right_maximal(x.first) << ", on idxT? " << idxT.is_right_maximal(x.second) << "\n";
      for(auto e : enumerateRight(idxS, x.first)){
	cout << e << ",";
      }
      cout << "\n";
      for(auto e : enumerateRight(idxT, x.second)){
	cout << e << ",";
      }
      cout << "\n";
    }
    //Have to take into the special case where we it is impossible to extend in one direction, but other direction might still have valid extensions left. If left side doesn't have any valid extensions, it will get filtered out on the next iteration before pushing anything into the stack. Only checking for enumerateLeft() could result in case where left side only has one possible extending character, and would not be reliable here.
    if(idxS.is_right_maximal(x.first) || idxT.is_right_maximal(x.second) ||
       (enumerateRight(idxS, x.first) != enumerateRight(idxT, x.second)) ||
       (enumerateRight(idxS, x.first).size()==1 && enumerateRight(idxS,x.first)[0] == BD_BWT_index<>::END)){  // Handle END symbols

      S.push(make_tuple(x.first,x.second, depth+1));
      if(verboseElement) cout << "Pushed x:" << x.first.toString() << x.second.toString() << "\n";
    }
  }
  return ret;
}
	       
vector<struct occStruct> batchLocate(vector<struct occStruct>  pairs, vector<bool> marked, BD_BWT_index<> bwt){
  int i = 0;
  int n = bwt.size();
  vector<struct occStruct> translate;
  vector<struct occStruct> ret;
  auto LFindex = mapLF(bwt).first;
  
  for(int j = n; j > 0; j--){
    if(marked[i]){
      struct occStruct temp;      
      temp.key = i;
      temp.primary = j-1;
      translate.push_back(temp);
    }
    i = LFindex[i];
  }  
  translate = radixSort(translate, 2);
  pairs = radixSort(pairs,2);

  int x = 0;
  int y = 0;
  while(x < pairs.size()){
    auto a = pairs[x];
    if(a.key == translate[y].key){
      struct occStruct temp;
      temp.key = translate[y].primary;
      temp.primary = a.primary;
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
  rsa1d T_a = rsa1d(size); 
  rsa1d T_b = rsa1d(size);
  rsa2d T_c = rsa2d(size);
  rsa2d T_d = rsa2d(size);
  vector<pair<int, int>> C_a(A.size(), make_pair(0,0));
  vector<pair<int, int>> C_b(A.size(), make_pair(0,0));
  vector<pair<int, int>> C_c(A.size(), make_pair(0,0));
  vector<pair<int, int>> C_d(A.size(), make_pair(0,0));
  vector<pair<int,pair<int,int>>> C_p(A.size(), make_pair(0,make_pair(0,0)));
  vector<pair<int, int>> C(A.size(), make_pair(0,0));
  vector<pair<int,int>> E_1;
  vector<pair<int,int>> E_2;
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
  T_c.sortArray();
  T_d.sortArray();
  T_a.sortArray();
  T_b.sortArray();
  if(verboseChaining) cout << "E_1 size: " << E_1.size() << endl;
  T_a.upgrade(make_pair(0,0),0);
  sort(E_1.begin(), E_1.end(), pairSort);

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
      
      auto a = T_a.rangeMax(0,			I.reverse.left); //I.reverse.left-1 causes failure when (reverse.left = forward.left = 0) with zeroth-index-indexing
      auto b = T_b.rangeMax(I.reverse.left,	I.reverse.right);
      auto c = T_c.rangeMax(INT_MIN,		rlfl,	0, I.forward.right);
      auto d = T_d.rangeMax(rlfl+1,		INT_MAX,0, I.reverse.right);
      
      C_a[j].first = a.first.primary.second;	C_a[j].second = a.second; 
      C_b[j].first = b.first.primary.second;	C_b[j].second = I.reverse.left + b.second; 
      C_c[j].first = c.first.primary.second;	C_c[j].second = I.forward.left + c.second;
      C_d[j].first = d.first.primary.second;	C_d[j].second = I.reverse.left + d.second;
      
      if(verboseChaining) cout << "maxCandinates: "   << C_a[j].second << "(" << C_a[j].first <<"),"
			       << C_b[j].second << "("<< C_b[j].first  << "),"
			       << C_c[j].second << "("<< C_c[j].first  << "),"
			       << C_d[j].second << "("<< C_d[j].first  << ")," << endl;
      
      auto max = chainingMax(C_a[j].second, C_b[j].second,C_c[j].second,C_d[j].second);
      if     (C_c[j].second == max) C[j] = C_c[j];
      else if(C_b[j].second == max) C[j] = C_b[j];
      else if(C_d[j].second == max) C[j] = C_d[j];
      else if(C_a[j].second == max) C[j] = C_a[j];

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