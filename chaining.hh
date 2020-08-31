#ifndef BDBWT_MEM_CHAIN_CHAIN_HH
#define BDBWT_MEM_CHAIN_CHAIN_HH
#include "rsa1d.hh"
#include "rsa2d.hh"
#include <omp.h>

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


/** Chaining function implements the chaining algorithm as defined in the paper "Chaining with Overlaps Revisited".
    The algorithm chains together all MEM matches given as input, and finds the best chaining alignment such that the chains do not themselves, while allowing overlaps on the text
*/
vector<pair<int,pair<int,int>>> chaining(vector<Interval_pair> A, int size);

/** Output function for the chains.
 */
pair<vector<Interval_pair>,vector<int>> chainingOutput(vector<pair<int,pair<int,int>>> chains, vector<Interval_pair> Ipairs, Configuration conf);


vector<pair<int,pair<int,int>>> chaining(vector<Interval_pair> A, int size){
  //cout << "beginning chaining for " << A.size() << " intervals...";
  rsa1d T_a = rsa1d(size); 
  rsa1d T_b = rsa1d(size);
  rsa1d T_c = rsa1d(size);
  rsa1d T_d = rsa1d(size);
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
    T_c.insertCell(key,INT_MIN);
    T_d.insertCell(key, INT_MIN);

    auto p1 = make_pair(A[j].forward.left, j);
    auto p2 = make_pair(A[j].forward.right, j);
    E_1.push_back(p1);
    E_1.push_back(p2);
  }
  T_a.sortArray();
  T_b.sortArray();
  T_c.sortArray();
  T_d.sortArray();
  T_a.upgrade(make_pair(0,0),0);
  sort(E_1.begin(), E_1.end(), pairSort);

  auto a = T_a.rangeMax(-1,-1); //Initializing variable types.
  auto b = T_b.rangeMax(-1,-1);
  auto c = T_c.rangeMax(-1,-1);
  auto d = T_d.rangeMax(-1,-1);
  //cout << "chaining initializer done" << endl;
  for(int i = 0; i < E_1.size(); i++){
    auto e = E_1[i];
    int j = e.second;
    Interval_pair I = A[j];
    auto rlfl = I.reverse.left - I.forward.left;
    if(I.forward.left == E_1[i].first){
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
          c = T_c.rangeMax(INT_MIN,		rlfl);
        }
#pragma omp section
        {
          d = T_d.rangeMax(rlfl+1,		INT_MAX);
        }
      }
      C_a[j].first = a.first.primary.second;	C_a[j].second = a.second; 
      C_b[j].first = b.first.primary.second;	C_b[j].second = I.reverse.left + b.second; 
      C_c[j].first = c.first.primary.second;	C_c[j].second = I.forward.left + c.second;
      C_d[j].first = d.first.primary.second;	C_d[j].second = I.reverse.left + d.second;

      auto max = chainingMax(C_a[j].second, C_b[j].second,C_c[j].second,C_d[j].second);
      if     (C_a[j].second == max) C[j] = C_a[j];
      else if(C_c[j].second == max) C[j] = C_c[j];
      else if(C_b[j].second == max) C[j] = C_b[j];
      else if(C_d[j].second == max) C[j] = C_d[j];

      auto cpsum = C[j].second+I.forward.right-I.forward.left+1;
      C_p[j] = make_pair(cpsum, make_pair(C[j].first,j));

      auto upgsumL = (int)C[j].second-(int)I.forward.left;
      T_c.upgrade(make_pair(rlfl, j), upgsumL);
      auto upgsumR = (int)C[j].second-(int)I.reverse.left;
      T_d.upgrade(make_pair(rlfl, j), upgsumR);
    }
    else{
      T_a.upgrade(make_pair(I.reverse.right, j), C_p[j].first);
      T_b.upgrade(make_pair(I.reverse.right, j), C[j].second- I.reverse.left);
      T_c.update(make_pair(rlfl, j), INT_MIN);
      T_d.update(make_pair(rlfl, j), INT_MIN);
    }
  }
  return C_p;
}
pair<vector<Interval_pair>,vector<int>> chainingOutput(vector<pair<int,pair<int,int>>> chains, vector<Interval_pair> Ipairs, Configuration conf){
  auto text = conf.text1;
  auto text2 = conf.text2;
  int maxIndex = 0;
  int maxVal = 0;
  for(int i = 0; i < chains.size(); i++){
    int tempVal = chains[i].first;
    if(tempVal > maxVal){
      maxVal = tempVal;
      maxIndex = i;
    }
  }
  if(conf.rawChains){ //printing raw chains
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
  for(int j = chains.size(); j >= 0; j--){
    if(i < 0){
      break;
    }
    if(chains.at(i).second.second >= 0){
      auto I = Ipairs.at(chains.at(i).second.second);
      if(chainIntervals.size() > 0 && chains.size() > 1 &&
         chainIntervals.back().forward.left >= I.forward.left &&
         chainIntervals.back().reverse.left >= I.reverse.left){ //Ensuring (weak) precedence

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
    int d = c.forward.right - c.forward.left;
    if(conf.verbosity > 0) cout << "Chain["<< count <<"]: "<< c.toString() <<"\t\t symcov:" << symcov.at(count) << endl;
    if(conf.chainStringSegments){
      cout << text.substr(a.left,d+1) <<",\t "<< text2.substr(b.left,d+1) << endl;
    }
    count++;
  }
  if(conf.verbosity > 0) cout << "total number of chains: " << count << ", total symmetric coverage: " << symcov.front() << endl;
  return (make_pair(chainIntervals, symcov));
}
#endif
