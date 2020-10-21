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
#include "RMaxQTree.h"


/** Chaining function implements the chaining algorithm as defined in the paper "Chaining with Overlaps Revisited".
    The algorithm chains together all MEM matches given as input, and finds the best chaining alignment such that the chains do not themselves, while allowing overlaps on the text
*/
vector<pair<int,pair<int,int>>> chaining(vector<Interval_pair> A, int size);

/** Output function for the chains.
 */
pair<vector<Interval_pair>,vector<int>> chainingOutput(vector<pair<int,pair<int,int>>> chains, vector<Interval_pair> Ipairs, Configuration conf);

vector<pair<int,pair<int,int>>> chainingNew(vector<Interval_pair> B){
  //cout << "chain init...B size = " << B.size();
  sort(B.begin(), B.end(), newChainingSort);

  vector<pair<int, pair<int,int>>> retVec;
  int n = B.size();
  int minkey_ab = INT_MAX;
  int maxkey_ab = INT_MIN;
  int minkey_cd = INT_MAX;
  int maxkey_cd = INT_MIN;
  int *keys_a = new int[n+1]; //Have to add extra space for 0 index.
  int *keys_b = new int[n+1];
  int *keys_c = new int[n];
  int *keys_d = new int[n];
  int Ap[n][4];
  int C[n], Cp[n][2];
  //retVec.push_back(make_pair(0,make_pair(0,0)));
  //cout << " retVec size " << retVec.size() << endl;
  for(int ra = 0; ra < n; ra++){
    auto I = B[ra];
    Ap[ra][0] = I.forward.left;
    Ap[ra][1] = I.forward.right;
    Ap[ra][2] = I.reverse.left;
    Ap[ra][3] = I.reverse.right;
    retVec.emplace_back(make_pair(0,make_pair(0,0)));
    //cout << "pushed " << ra << "elements retVec size = "<< retVec.size() << endl;
  }
  //cout << "retVec size 1 = " << retVec.size() << endl;
  keys_a[0] = 0;
  keys_b[0] = 0;
  minkey_ab = 0;
  maxkey_ab = 0;
  vector<int> keys1;
  vector<int> keys2;
  keys1.push_back(0);
  keys2.push_back(0);
  for(int ii = 0; ii < n; ii++){
    keys1.push_back(Ap[ii][3]);
    keys2.push_back(Ap[ii][2] - Ap[ii][0]);
  }
  sort(keys1.begin(), keys1.end());
  sort(keys2.begin(), keys2.end());
  for(int ri = 0; ri < n+1; ri++){
    keys_a[ri] = keys1[ri];
    keys_b[ri] = keys1[ri];
    if(keys1[ri] < minkey_ab) minkey_ab = keys1[ri]; //necessary since RMax returns -1 if query bounds are outside key-range.
    if(keys1[ri] > maxkey_ab) maxkey_ab = keys1[ri];
  }
  for(int ri = 0; ri < n; ri++){
    keys_c[ri] = keys2[ri];
    keys_d[ri] = keys2[ri];
    if(keys2[ri] < minkey_cd) minkey_cd = keys2[ri];
    if(keys2[ri] > maxkey_cd) maxkey_cd = keys2[ri];
     //cout << "ri loop " << ri << ", retVec size = "<< retVec.size() << endl;
  }
  //cout << "minmax keys: [" << minkey_ab << ", " << maxkey_ab << ",]" <<  ",[" << minkey_cd << ", " << maxkey_cd << ",]" << endl;
  // cout << "retVec size 2 = " << retVec.size() << endl;
  RMaxQTree *tree_a = new RMaxQTree(keys_a,n+1);
  RMaxQTree *tree_b = new RMaxQTree(keys_b,n+1);
  RMaxQTree *tree_c = new RMaxQTree(keys_c,n);
  RMaxQTree *tree_d = new RMaxQTree(keys_d,n);
  for(int ta = 0; ta < n; ta++) tree_a->update(keys_a[ta+1],ta,INT_MIN);
  for(int tb = 0; tb < n; tb++) tree_b->update(keys_b[tb+1],tb,INT_MIN);
  for(int tc = 0; tc < n; tc++) tree_c->update(keys_c[tc]  ,tc,INT_MIN);
  for(int td = 0; td < n; td++) tree_d->update(keys_d[td]  ,td,INT_MIN);

  tree_a->update(0,0,0);
  vector<pair<int,int>> E;
  for(int j = 0; j < n; j++){
    auto p1 = make_pair(Ap[j][0], j);
    auto p2 = make_pair(Ap[j][1], j);
    E.push_back(p1);
    E.push_back(p2);
  }
  sort(E.begin(), E.end());

  int ca[n], cb[n], cc[n], cd[n];
  //cout << "chain loop start" << endl;
  for(int i = 0; i < E.size(); i++){
    int j, e1;
    tie(e1, j) = E[i];
    //cout << "Ap[j] = (" << Ap[j][0] << "," << Ap[j][1] << "),("<< Ap[j][2] << "," << Ap[j][3]<< ")\t";
    auto rlfl = Ap[j][2]-Ap[j][0];
    //cout << "rlfl = " << rlfl << ", e1 = " << e1 << " j = "<< j <<endl;
    if(Ap[j][0] == e1){
      int ia, ib, ic, id, maxid = -102;
      //cout << "starting queries...";
#pragma omp parallel sections
      {
#pragma omp section
        {
          //a = T_a.rangeMax(0,			I.reverse.left); //I.reverse.left-1 causes failure when (reverse.left = forward.left = 0) with zeroth-index-indexing
          //cout << "a bounds" << minkey_ab << "," << Ap[j][2] << "\t";
          auto aret = tree_a->query(minkey_ab, Ap[j][2]-1);
          //cout << "a(" << aret.first << ", " << aret.second <<")..." << endl;
          ca[j] = aret.second;
          ia = aret.first;
        }
#pragma omp section
        {
          //b = T_b.rangeMax(I.reverse.left,	I.reverse.right);
          tie(ib, cb[j]) = tree_b->query(Ap[j][2],Ap[j][3]);
          cb[j] = Ap[j][2] + cb[j];
          //cout << "b(" << cb[j] << ", " << ib <<")...";
        }
#pragma omp section
        {
          //c = T_c.rangeMax(INT_MIN,		rlfl);
          tie(ic, cc[j]) = tree_c->query(minkey_cd,rlfl);
          cc[j] = Ap[j][0] + cc[j];
          //cout << "c(" << cc[j] << ", " << ic <<")...";
        }
#pragma omp section
        {
          //d = T_d.rangeMax(rlfl+1,		INT_MAX);
          tie(id,cd[j]) = tree_d->query(rlfl+1,maxkey_cd);
          cd[j] = Ap[j][2] + cd[j];
          //cout << "d(" << cd[j] << ", " << id <<")...";
        }
      }
      //cout << "done...";
      auto temp = {make_pair(ca[j],ia),
                   make_pair(cb[j],ib),
                   make_pair(cc[j],ic),
                   make_pair(cd[j],id)};
      tie(C[j],maxid) = max(temp);
      //cout << "C[j] = " << C[j];
      int maxSum = C[j] + Ap[j][1] - Ap[j][0] + 1;
      Cp[j][0] = maxSum;
      Cp[j][1] = maxid;
      //cout << "...max got...(" << maxSum << "," << maxid << ")..." << endl;
      //cout << "original retvector " << retVec[j].first << ",(" << retVec[j].second.first << "," << retVec[j].second.second << ")" << endl;
      retVec[j] = make_pair(maxSum, make_pair(maxid, j));
      //cout << "updated retvector " << retVec[j].first << ",(" << retVec[j].second.first << "," << retVec[j].second.second << ")" << endl;
      //cout << "got max: " << max << endl;

      //cout << "if-clause updates...";
      //cout << "upgrading c key " << rlfl << " to " << C[j]-Ap[j][0] << endl;
      tree_c->update(rlfl,j,C[j]-Ap[j][0]); //cout << "1..";
      //cout << "upgrading d key " << rlfl << " to " << C[j]-Ap[j][2] << endl;
      tree_d->update(rlfl,j,C[j]-Ap[j][2]); //cout << "2..Done, return loop" << endl << endl;
    }
    else{
      //cout << "else-clause updates...";
      //cout << "update a: " << Ap[j][3] << ", " << j << " into " << Cp[j][0] << endl;
      tree_a->update(Ap[j][3],j,Cp[j][0]); //cout << "1...";
      tree_b->update(Ap[j][3],j,C[j]-Ap[j][2]); //cout << "2...";

      tree_c->update(rlfl,j,INT_MIN); //cout << "3...";
      tree_d->update(rlfl,j,INT_MIN); //cout << "4...Done, return loop" << endl << endl;
    }
    //cout << endl;
  }
  //cout << "returning from new chaining" << endl;
  return retVec;
}
vector<pair<int,pair<int,int>>> chaining(vector<Interval_pair> A, int size){
  // cout << "beginning chaining for " << A.size() << " intervals...";

  rsa1d T_a = rsa1d(A.size()); 
  rsa1d T_b = rsa1d(A.size());
  rsa1d T_c = rsa1d(A.size());
  rsa1d T_d = rsa1d(A.size());
  vector<pair<int, int>> C_a(A.size(), make_pair(0,0));
  vector<pair<int, int>> C_b(A.size(), make_pair(0,0));
  vector<pair<int, int>> C_c(A.size(), make_pair(0,0));
  vector<pair<int, int>> C_d(A.size(), make_pair(0,0));
  vector<pair<int, pair<int,int>>> C_p(A.size(), make_pair(0,make_pair(0,0)));
  vector<pair<int, int>> C(A.size(), make_pair(0,0));
  vector<pair<int, int>> E_1;
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
  // cout << "chaining initializer done" << endl;
  for(int i = 0; i < E_1.size(); i++){
    // cout << "i = " << i << endl;
    auto e = E_1[i];
    int j = e.second;
    // cout << "j = " << j << endl;
    Interval_pair I = A[j];
    // cout << "Ipair = " << I.toString() << endl;
    auto rlfl = I.reverse.left - I.forward.left;
    if(I.forward.left == E_1[i].first){
      //  cout << "before rmax calls" <<  i << "/"  << E_1.size() <<  endl;
#pragma omp sections
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
      // cout << "after rmax calls" << endl;
       C_a[j].first = a.first.primary.second;	C_a[j].second = a.second;
       C_b[j].first = b.first.primary.second;	C_b[j].second = I.reverse.left + b.second;
       C_c[j].first = c.first.primary.second;	C_c[j].second = I.forward.left + c.second;
       C_d[j].first = d.first.primary.second;	C_d[j].second = I.reverse.left + d.second;

      auto max = chainingMax(C_a[j].second, C_b[j].second, C_c[j].second, C_d[j].second);
      if     (C_d[j].second == max) C[j] = C_d[j];
      else if(C_b[j].second == max) C[j] = C_b[j];
      else if(C_a[j].second == max) C[j] = C_a[j];
      else if(C_c[j].second == max) C[j] = C_c[j];


      auto cpsum = C[j].second+I.forward.right-I.forward.left+1;
      C_p[j] = make_pair(cpsum, make_pair(C[j].first,j));
      //cout << "updated retvector " << C_p[j].first << ",(" << C_p[j].second.first << "," << C_p[j].second.second << ")" << endl;

      auto upgsumL = (int)C[j].second-(int)I.forward.left;
      T_c.upgrade(make_pair(rlfl, j), upgsumL);
      //cout << "done upgrade 1" << endl;
      auto upgsumR = (int)C[j].second-(int)I.reverse.left;
      T_d.upgrade(make_pair(rlfl, j), upgsumR);
      //cout << "done upgrade 2" << endl;
    }
    else{
      T_a.upgrade(make_pair(I.reverse.right, j), C_p[j].first);
      T_b.upgrade(make_pair(I.reverse.right, j), C[j].second- I.reverse.left);
      //cout << "done upgrades in else 1" << endl;
      T_c.update(make_pair(rlfl, j), INT_MIN);
      T_d.update(make_pair(rlfl, j), INT_MIN);
      //cout << "done updates" << endl;
    }
    // cout << endl;
  }
  // cout << "returning from chaining" << endl;
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
  int i = maxIndex;
  //cout <<"Ret Chain["<< i << "]: " << chains[i].first << "," << chains[i].second.first << ":"<<chains[i].second.second << "\n";
  i = chains[maxIndex].second.first;
  for(int j = chains.size(); j >= 0; j--){
    if(i < 0 || last == 0){
      break;
    }
    if(chains.at(i).second.second >= 0){
      auto I = Ipairs.at(chains.at(i).second.second);
      if(chainIntervals.size() > 0 && chains.size() > 1
         //&&
         //&& chainIntervals.back().forward.left > I.forward.left
         //&& chainIntervals.back().reverse.left > I.reverse.left
         ){ //Ensuring (weak) precedence

        symcov.push_back(chains.at(i).first);
        chainIntervals.push_back(Ipairs.at(chains.at(i).second.second));
        last = i;
        //cout <<"Ret Chain["<< i << "]: " << chains[i].first << "," << chains[i].second.first << ":"<<chains[i].second.second << "\n";
        i = chains[i].second.first;
        if(i == 0){ //would print out same index again => chain is done.
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
