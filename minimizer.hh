#ifndef BDBWT_MEM_CHAIN_MINIMIZER_HH
#define BDBWT_MEM_CHAIN_MINIMIZER_HH
#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include "util.hh"
#include "sdsl/int_vector.hpp"
#include <omp.h>
using namespace std;
/** comparison function to determine the lexicographical order between two k-mers, or minimizers, given as pair<string,int>
 */
bool mimCompare (pair<string,int> first, pair<string,int> second);
/**  Computes all minimizers for input string t1. The computation is done by adding w k-mers into a bucket inheretently sorted by the set datastructure in cpp.
 For each window w, the lexicographically smallest k-mer is then selected as the minimizer.
*/
vector<pair<string,int>> minimizers(string t1, int k, int w);
/** Naive-ish way of computing which minimizers appear in both strings
 */
pair<vector<pair<string,int>>,vector<pair<string,int>>> mutualMinimizers(vector<pair<string,int>> m1, vector<pair<string,int>> m2);
/**  Translating minimizers into tuples tuple<i,j,d>, where:
 i = starting location of minimizer on text 1
 j = starting location of minimizer on text 2
 d = depth, or length of the minimizers.
*/
pair< vector<tuple<int,int,int>> ,pair< vector<pair<string,int>> , vector<pair<string,int>> > > minimizerTuples(vector<pair<string,int>> m1, vector<pair<string,int>> m2, Configuration conf, bool unique);

pair<vector<pair<string,int>>,vector<pair<string,int>>> mergeMinimizerPairs(vector<pair<string,int>> mini1, vector<pair<string,int>> mini2, Configuration conf);

/** Naive implementation to extend mutual minimizers left and right to obtain MEM match.
 */
vector<tuple<int,int,int>> memifyMinimizers(vector<tuple<int,int,int>> mini, Configuration conf);
/** Efficient computation of Permuted LCP array using BWT index and suffix array datastructures. As defined in paper "Permuted Longest-Common-Prefix Array"
 */
pair<vector<int>,vector<int>> createPLCP(BD_BWT_index<> index, int q, string t, bool lcp, vector<pair<int,int>> SA, bool direction);
/** Converting PLCP into LCP in linear time 
 */
vector<int> createLCPFromPLCP(vector<int> PLCP, vector<pair<int,int>> SA);
/** Creating an partition bitvector, such that B[i] = 1, iff: LCP[i] < k
 */
sdsl::int_vector<1> partitioning(int k, BD_BWT_index<> index, pair<vector<int>,vector<int>> PLCPLCP);
/** Using Suffix Array indexes to compute BWT intervals from minimizers, */
vector<pair<Interval_pair,string>> minimizerToBWTInterval(sdsl::int_vector<1> bv, sdsl::int_vector<1> bvr, vector<pair<string,int>> mini, vector<pair<int,int>> SA, vector<pair<int,int>> SAr,BD_BWT_index<> index, pair<vector<int>,vector<int>> lcp1, pair<vector<int>,vector<int>> lcp2, string text);


struct mimsort {
  bool operator() (const pair<string,int>& first, const pair<string,int>& second) const{
    return first.first < second.first;
  }
};
struct mimsortIndex {
  bool operator() (const pair<string,int>& first, const pair<string,int>& second) const{
    return first.second < second.second;
  }
};

bool mimCompare (pair<string,int> first, pair<string,int> second){
  if(first.first != second.first){
    return first.first < second.first;
  }
  else{
    return first.second < second.second;
  }
}
bool mimCompareIndex (pair<string,int> first, pair<string,int> second){
  return first.second < second.second;
}


struct minimizerStruct {
  string minimizer;
  vector<int> occ1;
  vector<int> occ2;
};
pair<vector<pair<string,int>>,vector<pair<string,int>>> minimizerBlobbingUnique(vector<pair<string,int>> m1, vector<pair<string,int>> m2){
  vector<pair<string,int>> retA;
  vector<pair<string,int>> retB;
  unordered_map<string, pair<vector<int>,vector<int>>> minimap;
  auto smaller = ((m1.size() < m2.size())? m1.size() : m2.size());
  //cout << "Smaller =" << smaller << endl;
  for(int i = 0; i < smaller; i++){
    minimap[m1[i].first].first.push_back(m1[i].second);
    minimap[m2[i].first].second.push_back(m2[i].second);
  }
  if(m1.size() > smaller){
    for(int i = smaller; i < m1.size(); i++){
      minimap[m1[i].first].first.push_back(m1[i].second);
    }
  }else if(m2.size() > smaller){
    for(int i = smaller; i < m2.size(); i++){
      minimap[m2[i].first].second.push_back(m2[i].second);
    }
  }
  for(const auto &c : minimap){ //O(|n|)
    auto minLen = c.first.length();
    if(!(c.second.first.size() > 0 && c.second.second.size() > 0)){ //O(1)
      continue;
    }else{
      // cout << c.first << ", " << c.second.first[0] << endl;
      //cout << c.first << ", " << c.second.second[0] << endl;
      //cout << endl;
      retA.push_back(make_pair(c.first,c.second.first[0]));
      retB.push_back(make_pair(c.first,c.second.second[0]));
    }
  }
  //cout << "Done blobs" << endl;
  return make_pair(retA,retB);
}


vector<tuple<int,int,int>> minimizerBlobbing(Configuration conf, vector<pair<string,int>> m1, vector<pair<string,int>> m2){
  vector<tuple<int,int,int>> retTuple;
  map<string,vector<int>> minimapF;
  map<string,vector<int>> minimapR;
  auto smaller = ((m1.size() < m2.size())? m1.size() : m2.size());
  //cout << "Smaller =" << smaller << endl;
  for(int i = 0; i < smaller; i++){
      minimapF[m1[i].first].push_back(m1[i].second);
      minimapR[m2[i].first].push_back(m2[i].second);

 /*     if(m1[i].first != conf.text1.substr(m1[i].second, m1[i].first.length())){
	      cout << "mismatch 1" << endl;
      }else cout << "Pushed 1 " << conf.text1.substr(m1[i].second,m1[i].first.length()) << " With key " << m1[i].first << endl;

      if(m2[i].first != conf.text2.substr(m2[i].second, m2[i].first.length())){
	      cout << "mismatch 2" << endl;
      }else cout << "Pushed 2 " << conf.text2.substr(m2[i].second,m2[i].first.length()) << " With key " << m2[i].first << endl;
  */
  }
  if(m1.size() > smaller){
    for(int i = smaller; i < m1.size(); i++){
      minimapF[m1[i].first].push_back(m1[i].second);
  /*    if(m1[i].first != conf.text1.substr(m1[i].second, m1[i].first.length())){
	      cout << "mismatch 3" << endl;
      }else cout << "Pushed 1 " << conf.text1.substr(m1[i].second,m1[i].first.length()) << " With key " << m1[i].first << endl;
   */}
  }else if(m2.size() > smaller){
    for(int i = smaller; i < m2.size(); i++){
      minimapR[m2[i].first].push_back(m2[i].second);
      /*if(m2[i].first != conf.text2.substr(m2[i].second, m2[i].first.length())){
	      cout << "mismatch 4" << endl;
      }else cout << "Pushed 2 " << conf.text2.substr(m2[i].second,m2[i].first.length()) << " With key " << m2[i].first << endl;
    */}
  }

  // for(auto a : m1){
  //   minimap[a.first].first.push_back(a.second);
  // }
  // for(auto b : m2){
  //   minimap[b.first].second.push_back(b.second);
  // }
 
  // for(const auto &c : minimap){
  //   //cout << "Blob: " << c.first;
  //   cout << "First: ";
  //   for(auto e : c.second.first){
  //     cout << e << ", ";
  //   }
  //   cout << "Second: ";
  //   for(auto e : c.second.second){
  //     cout << e << ", ";
  //   }
  //   cout << endl;
  // }
  // Blobs to tuples

  //vector<int> crossVector (m1.size()*m2.size(),-1);
  bool fancyType = false;
  if(fancyType){
/*
      for(const auto &c : minimap){
        int alpha = c.second.first.size();
        int beta = c.second.second.size();
        if(c.second.first.size() > 0 && c.second.second.size() > 0){
          //int size = alpha*beta;
          //cout << "size: " << alpha << ", " << beta << endl;
          int size = (alpha*beta > alpha+beta)? alpha*beta : alpha+beta;
          vector<int> crossVector ((size)*2,-1);
          int i = 0;
          int j = 1;
          int d = alpha-1+beta-1;
#pragma omp sections
    {
#pragma omp section
          {
            for(auto i1 : c.second.first){
              for(int x = 0; x < beta; x++){
                //cout << "i = " << i << endl;
                crossVector[i] = i1;
                i = i+2;
              }
              //cout << "---" << endl;
            }
          }
#pragma omp section
          {
            for(auto i2 : c.second.second){
              for(int x = 0; x < alpha; x++){
                //cout << "j+(x*d)= " << j+(x*d) << endl;
                crossVector[j + ((x*beta))] = i2;
              }
              j = j+2;
              //cout << "---" << endl;
            }
          }
    }
          //cout << "vector done" << endl;
          int k = 0;
          while(k < crossVector.size()){
            retTuple.push_back(make_tuple(crossVector[k], crossVector[k+1],c.first.length()));
            k = k+2;
          }
        }
      }
      */
  }else{
    for(auto c : minimapF){ //O(|n|)
      auto minLen = c.first.length();
      //if(!(c.second.first.size() > 0 && c.second.second.size() > 0)){ //O(1)
      //  continue;
      //}
      auto vec1 = c.second;
      auto vec2 = minimapR[c.first];
      int li1 = 0;
      int li2 = 0;
      int lminlen = 0;
      for(int a = 0; a < vec1.size();a++){ //O(k log n)
        for(int b = 0; b < vec2.size();b++){ 
	  auto i1 = vec1[a];
	  auto i2 = vec2[b];
   	  if(li1 == i1 && li2 == i2 && lminlen == minLen){continue;}
         // cout << "BlobTuple: " << i1 << ", " << i2 << ", " << minLen <<" Min: " << c.first << endl;
          retTuple.push_back(make_tuple(i1,i2,minLen-1)); //O(1)
	  li1 = i1;
	  li2 = i2;
	  lminlen = minLen;
	  auto str1 = conf.text1.substr(i1,minLen);
	  auto str2 = conf.text2.substr(i2,minLen);

	  /*cout << str1 << "\n <> \n" << str2 << endl;
	  if(str1 != str2){
		  cout << "Strings are not the same" << endl;
	  }else{
		  cout << "Strings are the same" << endl;
	  }cout << endl;
	  */
        }
      }
    }
  }
  //cout << "Done blobs" << endl;
  return retTuple;
}
vector<pair<string,int>> minimizers(string t1, int k, int w){
  //  cout << "enter minimizer";
  vector<pair<string,int>> kmers((t1.size()));
  vector<pair<string,int>> ret;
  set<pair<string,int>,mimsortIndex> minimizers;
  vector<pair<string,int>>::iterator minimizer;

  for(int i = 0; i < t1.size()-k; i++){
    kmers.push_back(make_pair(t1.substr(i,k),i));
    //minimizers.insert(kmers.back());
    if(i >= w){
      minimizer = min_element(kmers.end()-w,kmers.end(),mimCompare);
      auto mini = *minimizer;
      minimizers.insert(mini);
      }
  }
  // cout << "minimizers: " << endl;
  for(auto km : minimizers){
    //  cout << km.first << ", " << km.second << endl;
    ret.push_back(km);
  }
  //cout << "found " << ret.size() << " minimizers" << endl;
  //sort(ret.begin(), ret.end(), minimizerLexSort);
  return ret;
}
pair<vector<pair<string,int>>,vector<pair<string,int>>> mutualMinimizers(vector<pair<string,int>> m1, vector<pair<string,int>> m2){
  vector<pair<string,int>> m3, m4;
  int b = 0;
  int c = 0;
  auto newtype = true;
  if(!newtype){
    for(int a = 0; a < m1.size(); a++){
      auto first = m1[a];
      while(m2.at(b).first < first.first && b < m2.size()-1){
        b++;
      }
      c = 0;
      while(b+c < m2.size() && m2.at(b+c).first == first.first){
        if(m2.at(b+c).first == first.first){
          m3.push_back(first);
          m4.push_back(m2.at(b+c));
          c++;
        }else{
          c++;
        }
	    }
    }
    return make_pair(m3,m4);
  }else{
    int i = 0;
    int j = 0;
    while(i < m1.size()-1 || j < m2.size()-1){
      if(i < m1.size()-1 && j < m2.size()-1 && m1.at(i).first == m2.at(j).first){
        m3.push_back(m1.at(i));
        m4.push_back(m2.at(j));
        i++; j++;
      }else if(j < m2.size() && m1.at(i) > m2.at(j)){
        j++;
      }else if(i < m1.size() && m1.at(i) < m2.at(j)){
        i++;
      }
    }
    return(make_pair(m3,m4));
  }
 }
pair<vector<pair<string,int>>,vector<pair<string,int>>> minimizerAnchors(vector<pair<string,int>> m1, vector<pair<string,int>> m2){
  vector<pair<string,int>> m3, m4;
  int b = 0;
  int i = 0;
  auto newtype = true;
  if(newtype){
    for(auto first : m1){
      if(m3.size() > 0){
        if(m3.back() == first){
          continue;
        }
      }
    if(b < m2.size()-1){
        while(m2.at(b).first < first.first){
          b++;
        }
        if(b+i < m2.size()-1 && m2.at(b+i).first == first.first){
          m3.push_back(first);
          m4.push_back(m2.at(b));
          i++;
        }else{
          i = 0;
        }
	    }
    }
    return make_pair(m3,m4);

  }else{
  unordered_set<string> m1Kmer(m1.size()-1);
  unordered_set<string> m2Kmer(m2.size()-1);
  unordered_set<string> anchorMers;
  vector<pair<string,int>> m3, m4;
  for(auto first : m1){
    m1Kmer.insert(first.first);
  }
  for(auto first : m2){
    m2Kmer.insert(first.first);
  }
  for(auto m : m1Kmer){
    if(m2Kmer.find(m) != m2Kmer.end()){
      anchorMers.insert(m);
    }
  }
  auto temp = anchorMers;
  for(auto m : m1){
    if(anchorMers.find(m.first) != anchorMers.end()){
      m3.push_back(m);
      anchorMers.erase(m.first);
    }
  }
  anchorMers = temp;
  for(auto m : m2){
    if(anchorMers.find(m.first) != anchorMers.end()){
      m4.push_back(m);
      anchorMers.erase(m.first);
    }
  }
  //auto ret = mergeMinimizerPairs(m3,m4);
  return make_pair(m3,m4);
  }
}

pair< vector<tuple<int,int,int>> ,pair< vector<pair<string,int>> , vector<pair<string,int>> > > minimizerTuples(vector<pair<string,int>> m1, vector<pair<string,int>> m2, Configuration conf, bool unique = false){
  vector<tuple<int,int,int>> retTuple;
  pair<vector<pair<string,int>>,vector<pair<string,int>>> retRaw;

  sort(m1.begin(), m1.end(), mimCompare);
  sort(m2.begin(), m2.end(), mimCompare);

  sdsl::int_vector<1> bv(conf.text1.length()+1, 0, 1);
  sdsl::rank_support_v<1> b_rank(&bv);
  sdsl::select_support_mcl<1> b_select(&bv);

  sdsl::int_vector<1> bv2(conf.text2.length()+1, 0, 1);
  sdsl::rank_support_v<1> b_rank2(&bv2);
  sdsl::select_support_mcl<1> b_select2(&bv2);
  pair<vector<pair<string,int>>,vector<pair<string,int>>>  muts;
  if(unique){
    muts = minimizerAnchors(m1,m2);
  }else{
    muts = mutualMinimizers(m1,m2);
  }
  m1 = muts.first;
  m2 = muts.second;
  if(conf.miniMergerCount > 0){
    sort(m1.begin(),m1.end(), mimCompareIndex);
    sort(m2.begin(),m2.end(), mimCompareIndex);
    pair<vector<pair<string,int>>,vector<pair<string,int>>> merger = make_pair(m1,m2);
    for(int i = 0; i < conf.miniMergerCount; i++){
      merger = mergeMinimizerPairs(merger.first,merger.second,conf);
    }
    m1 = merger.first;
    m2 = merger.second;
    if(conf.verbosity > 2) cout << "sorted minimems" << endl;
    sort(m1.begin(),m1.end(), mimCompare);
    sort(m2.begin(),m2.end(), mimCompare);
  }
  for(int i = 0; i < m1.size()-1; i++){
    auto x = m1[i];
    auto y = m2[i];
    if(x.first.compare(y.first) == 0){
      if(unique){
        int end = ((conf.text1.size()-1 < x.second + x.first.length()-1+20)? 1 : x.first.length()-1+20);
        int end2 = ((conf.text2.size()-1 < y.second + y.first.length()-1+20)? 1 : y.first.length()-1+20);
        auto r1 = b_rank(x.second+end) - b_rank(x.second+1);
        auto r2 = b_rank2(y.second+end2) - b_rank2(y.second+1);
        if(r1 > 0 || r2 > 0){
          continue;
        }
      }
      auto tup = make_tuple(x.second, y.second, x.first.length());
      retRaw.first.push_back(x);
      retRaw.second.push_back(y);
      retTuple.push_back(tup);

      if(unique){
        for(int a = x.second; a < x.second + x.first.length(); a++){
          bv[a] = 1;
        }
        for(int b = y.second; b < y.second + y.first.length(); b++){
          bv2[b] = 1;
        }
      }
    }
  }
  // cout << "ret muts" << endl;
  return make_pair(retTuple,retRaw);
}
pair<vector<pair<string,int>>,vector<pair<string,int>>> mergeMinimizerPairs(vector<pair<string,int>> mini1, vector<pair<string,int>> mini2, Configuration conf){
  int a,b,c,i,j,k;
  set<pair<string,int>> retMini1;
  set<pair<string,int>> retMini2;
  if(conf.verbosity > 2) cout << "merger" << endl;
  int count = 0;
  for(int x = 0; x < mini1.size(); x++){
    for(int y = 0; y < mini2.size(); y++){
      if(mini1.at(x).first.compare(mini2.at(y).first) != 0) continue;
      if(x > 0 && y > 0){
        a = mini1.at(x-1).second;
        b = mini2.at(y-1).second;
        c = mini1.at(x-1).first.length();

        i = mini1.at(x).second;
        j = mini2.at(y).second;
        k = mini1.at(x).first.length();
      }
      if((x > 0 && y > 0) && (i-a == j-b && a+c >= i && b+c >= j && a+c <= i+k && b+c <= j+k)){
        string f = mini1.at(x-1).first;
        string e = mini1.at(x).first;
        int d = a+c-i;
        int l = abs(d)+1;
        string merged = conf.text1.substr(a,c+k-l);
        if(merged.compare(conf.text2.substr(b,c+k-l)) != 0) continue;
        count++;
        retMini1.insert(retMini1.end(),make_pair(merged,a));
        retMini2.insert(retMini2.end(),make_pair(merged,b));
      }
      else{
        retMini1.insert(retMini1.end(),mini1.at(x));
        retMini2.insert(retMini2.end(),mini2.at(y));
      }
    }
  }
  vector<pair<string,int>> retMiniVec1, retMiniVec2;
  for(auto v : retMini1){
    retMiniVec1.push_back(v);
  }
  for(auto v : retMini2){
    retMiniVec2.push_back(v);
  }
  cout << "merged " << count << " times" << endl;
  return make_pair(retMiniVec1,retMiniVec2);
}

vector<tuple<int,int,int>> memifyMinimizers(vector<tuple<int,int,int>> mini, Configuration conf){
  set<tuple<int,int,int>> miniMemsSet;
  vector<tuple<int,int,int>> miniMems;
  vector<vector<tuple<int,int,int>>> miniMemThreads(omp_get_max_threads());
#pragma omp parallel for
  for(int a = 0; a < mini.size(); a++){
    int i,j,k;
    tie(i,j,k) = mini.at(a);
    int i2,j2,k2;
    i2 = i;
    j2 = j;
    k2 = k;
    while(i2-1 >= 0 && j2-1 >= 0 && conf.text1.at(i2-1) == conf.text2.at(j2-1)){
      i2 = i2-1;
      j2 = j2-1;
      k2 = k2+1;
    }
    i = i2; j = j2; k = k2;
    while(i+k+1 < conf.text1.length() && j+k+1 < conf.text2.length() && conf.text1.at(i+k+1) == conf.text2.at(j+k+1)){
      k++;
    } 
    miniMemThreads[omp_get_thread_num()].push_back(make_tuple(i,j,k));
  }
  for(auto mm : miniMemThreads){
    for(auto m2 : mm){
      miniMemsSet.insert(m2);
    }
  }
  miniMems.assign(miniMemsSet.begin(), miniMemsSet.end());
  return miniMems;
}
vector<int> createLCPFromPLCP(vector<int> PLCP, vector<pair<int,int>> SA){
  vector<int> LCP(PLCP.size());
  for(int i = 0; i < PLCP.size(); i++){
    LCP[i] = PLCP[SA[i].first];
  }
  return LCP;
}
/** Works correctly, not optimal for space as per reference paper.
    By definition: PLCP[SA[j]] = LCP[j].
*/
pair<vector<int>,vector<int>> createPLCP(BD_BWT_index<> index, int q, string text, bool lcp, vector<pair<int,int>> SA, bool direction){
  vector<int> phi(SA.size());
  vector<int> PLCP(SA.size());
  auto LF = mapLF(index, direction);
  string t = text;
  t.append("$");
  for(int j = 1; j < SA.size(); j++){
    if(SA[j].first % q == 0){
      phi[SA[j].first] = (SA[j-1].first);
    }
  }
  int ell = 0;
  for(int i = 0; i < floor((SA.size()-1)/q); i++){
    int j = phi[i];
    while(t.at(i*q+ell) == t.at(j+ell)){
      ell = ell+1;
    }
    PLCP[i] = ell;
    ell = (((ell-q) > 0)? (ell-q) : 0);
  }
  if(lcp){
    auto LCP = createLCPFromPLCP(PLCP, SA);
    return (make_pair(PLCP, LCP));
  }
  vector<int> dummy(1,-1);
  return make_pair(PLCP, dummy);
}

sdsl::int_vector<1> partitioning(int k, BD_BWT_index<> index, pair<vector<int>,vector<int>> PLCPLCP){
  auto PLCP = PLCPLCP.first;
  auto LCP = PLCPLCP.second;

  sdsl::int_vector<1> B(PLCP.size()+1, 0, 1);
  B[0] = 1;
  B[PLCP.size()] = 1;
  for(int i = 0; i < PLCP.size(); i++){
    if(LCP[i] < k){
      B[i] = 1;
    }
  }
  return B;
}
vector<pair<Interval_pair,string>> minimizerToBWTInterval(sdsl::int_vector<1> bv, sdsl::int_vector<1> bvr, vector<pair<string,int>> mini, vector<pair<int,int>> SA, vector<pair<int,int>> SAr,BD_BWT_index<> index, pair<vector<int>,vector<int>> lcp1, pair<vector<int>,vector<int>> lcp2, string text = ""){
  // cout << "Enter minimizerToBWTInterval" << endl;
  sdsl::rank_support_v<1> b_rank(&bv);
  sdsl::rank_support_v<1> b_rankr(&bvr);
  sdsl::select_support_mcl<1> b_select(&bv);
  sdsl::select_support_mcl<1> b_selectr(&bvr);

  vector<pair<Interval_pair,string>> P;
  unordered_set<int> stored;
  for(int j = 0; j < mini.size(); j++){
    if(mini[j].second + mini[j].first.length() > text.length()) continue;
    int i = mini[j].second; //mini[j].second denotes the index of the k-mer on the original k-mer listing (all k-mers of text).
    int ir = (text.length()-1)-i;
    int revOffset (mini[j].first.length()-1);
    auto s  = SA[i].second;
    auto s2 = SAr[ir-revOffset].second;

    int r1,r2,r3,r4;
    if(bv[s] == 1){
      r1 = s;
    }else{
      r1 = b_select(b_rank(s));
    }
    r2 = b_rank(r1)+1;

    if(bvr[s2] == 1){
      r3 = s2;
    }else{
      r3 = b_selectr(b_rankr(s2));
    }
    r4 = b_rankr(r3)+1;

    int a,b,c,d;
    a = r1;
    b = b_select(r2+1)-1;
    c = r3;
    d = b_selectr(r4+1)-1;
    P.push_back(make_pair(Interval_pair(a,b,c,d),mini[j].first));
    stored.insert(a);
  }
  return P;
}
#endif
 
