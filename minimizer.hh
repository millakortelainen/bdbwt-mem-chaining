#ifndef BDBWT_MEM_CHAIN_MINIMIZER_HH
#define BDBWT_MEM_CHAIN_MINIMIZER_HH
#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include "util.hh"
#include "sdsl/int_vector.hpp"
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
pair< vector<tuple<int,int,int>> ,pair< vector<pair<string,int>> , vector<pair<string,int>> > > minimizerTuples(vector<pair<string,int>> m1, vector<pair<string,int>> m2, bool unique, string text1, string text2);

vector<tuple<int,int,int>> mergeMinimizers(vector<tuple<int,int,int>> mini, string text1, string text2);
pair<vector<pair<string,int>>,vector<pair<string,int>>> mergeMinimizerPairs(vector<pair<string,int>> mini1, vector<pair<string,int>> mini2,string text1, string text2);

/** Naive implementation to extend mutual minimizers left and right to obtain MEM match.
 */
vector<tuple<int,int,int>> memifyMinimizers(vector<tuple<int,int,int>> mini, string text, string text2);
/** Efficient computation of Permuted LCP array using BWT index and suffix array datastructures. As defined in paper "Permuted Longest-Common-Prefix Array"
 */
pair<vector<int>,vector<int>> createPLCP(BD_BWT_index<> index, int q, string t, bool lcp, vector<pair<int,int>> SA, bool direction);
/** Converting PLCP into LCP in linear time 
 */
vector<int> createLCPFromPLCP(vector<int> PLCP, vector<pair<int,int>> SA);
/** Creating an partition bitvector, such that B[i] = 1, iff: LCP[i] < k
 */
sdsl::int_vector<1> partitioning(int k, BD_BWT_index<> index, pair<vector<int>,vector<int>> PLCPLCP);
/** Using Suffix Array indexes to compute BWT intervals from minimizers, thus allowing use of BDBWT to finalize MEM matching.
 */
vector<Interval_pair> minimizerToBWTIntervalV2(vector<pair<string,int>> mini, vector<pair<int,int>> SA, vector<pair<int,int>> SAr,BD_BWT_index<> index, string text);
/** Using Suffix Array indexes to compute BWT intervals from minimizers, thus allowing use of BDBWT to finalize MEM matching. Unused, left for completeness. Implementation contains issues on correctness of results.
 */
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
  return first.first < second.first;
}
bool mimCompareIndex (pair<string,int> first, pair<string,int> second){
  return first.second < second.second;
}

vector<pair<string,int>> minimizers(string t1, int k, int w){
  //  cout << "enter minimizer";
  vector<pair<string,int>> kmers((t1.size()));
  vector<pair<string,int>> ret;
  set<pair<string,int>,mimsortIndex> minimizers;
  vector<pair<string,int>>::iterator minimizer;

  for(int i = 0; i < t1.size()-k; i++){
    kmers.push_back(make_pair(t1.substr(i,k),i));
    if(i >= w){
      minimizer = min_element(kmers.end()-w,kmers.end(),mimCompare);
      auto mini = *minimizer;
        minimizers.insert(mini);
      }
  }
  // cout << "minimizers: " << endl;
  for(auto km : minimizers){
  //  cout << km.first << endl;
    ret.push_back(km);
  }
  cout << "found " << ret.size() << " minimizers" << endl;
  sort(ret.begin(), ret.end(), minimizerLexSort);
  return ret;
}

pair<vector<pair<string,int>>,vector<pair<string,int>>> mutualMinimizers(vector<pair<string,int>> m1, vector<pair<string,int>> m2){
  // unordered_set<string> m1Kmer(m1.size()-1);
  // unordered_set<string> m2Kmer(m2.size()-1);
  // set<string> anchorMers;
  vector<pair<string,int>> m3, m4;
  int b = 0;
  int i = 0;
  auto newtype = true;
  if(newtype){
    for(auto first : m1){
	    if(b < m2.size()-1){
      while(m2.at(b).first < first.first){
        b++;
      }
      if(b+i < m2.size() && m2.at(b+i).first == first.first){
        m3.push_back(first);
        m4.push_back(m2.at(b));
        i++;
      }else{
        i = 0;
      }
	    }	
    }
    //cout << "hue" << endl;
    //auto ret = mergeMinimizerPairs(m3,m4);
    return make_pair(m3,m4);
  }else{
    // for(auto first : m2){
    //   m2Kmer.insert(first.first);
    // }
    // for(auto m : m1Kmer){
    //   if(m2Kmer.find(m) != m2Kmer.end()){
    //     anchorMers.insert(m);
    //   }
    // }
    // // for(auto m: anchorMers){
    // //   cout << "anchorMer: "<< m << endl;
    // // }
    // int maxSize = (m1.size()-1 > m2.size()-1)? m1.size()-1 : m2.size()-1;
    // for(int i = 0; i < maxSize; i++){
    //   if((m1.size()-1 > maxSize) && (anchorMers.find(m1.at(i).first) == anchorMers.end())){
    //     m1.erase(m1.begin()+i);
    //   }
    //   if((m2.size()-1 > maxSize) && (anchorMers.find(m2.at(i).first) == anchorMers.end())){
    //     m2.erase(m2.begin()+i);
    //   }
    // }
    // return make_pair(m1,m2);
  }
}
pair<vector<pair<string,int>>,vector<pair<string,int>>> minimizerAnchors(vector<pair<string,int>> m1, vector<pair<string,int>> m2){
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

pair< vector<tuple<int,int,int>> ,pair< vector<pair<string,int>> , vector<pair<string,int>> > > minimizerTuples(vector<pair<string,int>> m1, vector<pair<string,int>> m2, bool unique = false, string text1 = "", string text2 = ""){
  vector<tuple<int,int,int>> retTuple;
  pair<vector<pair<string,int>>,vector<pair<string,int>>> retRaw;

  sort(m1.begin(), m1.end(), mimCompare);
  sort(m2.begin(), m2.end(), mimCompare);
  pair<vector<pair<string,int>>,vector<pair<string,int>>>  muts;
  if(unique){
    muts = minimizerAnchors(m1,m2);
  }else{
    muts = mutualMinimizers(m1,m2);
  }
  m1 = muts.first;
  m2 = muts.second;
  bool doMerger = false;
  if(doMerger){
    cout << "test" << endl;
    //sort(retTuple.begin(),retTuple.begin(),memSort);
    //retTuple = mergeMinimizers(retTuple,text1,text2);
    sort(m1.begin(),m1.end(), mimCompareIndex);
    sort(m2.begin(),m2.end(), mimCompareIndex);
    auto merger = mergeMinimizerPairs(m1,m2,text1,text2);
    cout << "merger done" << endl;
    m1 = merger.first;
    m2 = merger.second;
    cout << "sorted minimems" << endl;
    sort(m1.begin(),m1.end(), mimCompare);
    sort(m2.begin(),m2.end(), mimCompare);
  }
  for(auto x : m1){
    int i = 0;
    for(auto y : m2){
      i++;
      //cout << "Comparing: " << x.first << " & " << y.first;
      if(x.first.compare(y.first) == 0){
        // m2.erase(m2.begin(), m2.begin()+i);
        i = 0;
        // if(retRaw.first.size() > 0){
        //   if(retRaw.first.back() == retRaw.second.back()){
        //     int a,b,c;
        //     int l = x.first.length();
        //     tie(a,b,c) = retTuple.back();

        //     if(x.second - l <= a+c && x.second-l >= a &&
        //        y.second - l <= b+c && y.second-l >= b){
        //       //retTuple.pop_back();
        //       //retTuple.push_back(make_tuple(a,b,(l+c)-((a+c+1)-(x.second))));
        //     }
        //   }
        // }
        auto tup = make_tuple(x.second, y.second, x.first.length());;
        //cout << "pushed";
        retRaw.first.push_back(x);
        retRaw.second.push_back(y);
        retTuple.push_back(tup);
        //break;
      }
      //      cout << endl;
      // if(x.first.length() > 0 && y.first.length() > 0 && x.first.at(0) != y.first.at(0)){
      //	break;
      //}
    }
  }

  cout << "ret muts" << endl;
  return make_pair(retTuple,retRaw);
}
vector<tuple<int,int,int>> mergeMinimizers(vector<tuple<int,int,int>> mini, string text1, string text2){  
  int a,b,c,i,j,k;
  for(int x = 1; x < mini.size(); x++){
    tie(a,b,c) = mini.at(x-1);	
    tie(i,j,k) = mini.at(x);
    if(i-a == j-b && a+c >= i && b+c >= j){
      mini.at(x) = make_tuple(a,b,c+abs(c-k));
      mini.erase(mini.begin()+x-1);
      x--;
    }
  }
  return mini;
}
pair<vector<pair<string,int>>,vector<pair<string,int>>> mergeMinimizerPairs(vector<pair<string,int>> mini1, vector<pair<string,int>> mini2, string text1, string text2){
  int a,b,c,i,j,k;
  cout << "merger" << endl;
  for(auto m : mini1){
	  cout << m.first << "," << m.second << endl;
 }
  for(int x = 1; x < mini1.size(); x++){
    if(mini1.at(x).first.compare(mini2.at(x).first) != 0) continue;
    if(mini1.at(x-1).first.compare(mini2.at(x-1).first) != 0) continue;

    a = mini1.at(x-1).second;
    b = mini2.at(x-1).second;
    c = mini1.at(x-1).first.length()-1;

    i = mini1.at(x).second;
    j = mini2.at(x).second;
    k = mini1.at(x).first.length()-1;
    int xlen = c-(i-(a+c))+k+1;
    int ylen = c-(j-(b+c))+k+1;

    if(a+c < i && b+c< j || c < 1 || k < 1)  continue;
    //cout << "dist = " << i-(a+c+1) << ",(" << a << "," << i << ")" << endl; 
    if(i-a == j-b && a+c >= i && b+c >= j){
	    string f = mini1.at(x-1).first;
	    string e = mini1.at(x).first;
	    int d = a+c-i;
	    int l = abs(d)+1;
	    string merged = text1.substr(a,c+abs(c-k));
	    cout << "merged: " << merged << " (" << a << "," << i << ")" <<endl;

      mini1.at(x) = make_pair(merged,a);
      mini1.erase(mini1.begin()+x-1);
      mini2.at(x) = make_pair(merged,b);
      mini2.erase(mini2.begin()+x-1);
      x--;
    }
  }
  return make_pair(mini1,mini2);
}
vector<tuple<int,int,int>> memifyMinimizers(vector<tuple<int,int,int>> mini, string text, string text2){
  set<tuple<int,int,int>> miniMemsSet;
  vector<tuple<int,int,int>> miniMems;
  mini = mergeMinimizers(mini, text, text2);
  for(auto m : mini){
    int i,j,k;
    tie(i,j,k) = m;
    int i2,j2,k2;
    i2 = i;
    j2 = j;
    k2 = k;
    while(i2-1 >= 0 && j2-1 >= 0 && text.at(i2-1) == text2.at(j2-1)){
      i2 = i2-1;
      j2 = j2-1;
      k2 = k2+1;
    }
    i = i2; j = j2; k = k2;
    while(i+k < text.length() && j+k < text2.length() && text.at(i+k) == text2.at(j+k)){
      k++;
    }
    miniMemsSet.insert(make_tuple(i,j,k));
  }

  miniMems.assign(miniMemsSet.begin(), miniMemsSet.end());
  // for(auto m : miniMems){
  //   int i,j,k;
  //   tie(i,j,k) = m;
  //   cout << "minimem: (" << i << "," << j << "," << k << ") "<<endl;
  // }
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
  //cout << "SA size: " << SA.size() << endl;
  auto LF = mapLF(index, direction);
  string t = text;
  t.append("$");
  //cout << "enter plcp" << endl;
  for(int j = 1; j < SA.size(); j++){
    if(SA[j].first % q == 0){
      phi[SA[j].first] = (SA[j-1].first);
    }
  }
  int ell = 0;
  for(int i = 0; i < floor((SA.size()-1)/q); i++){
    int j = phi[i];
    //    cout << "j = " << j << ", i = " << i << endl;
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
  //  cout << "Enter partitioning" << endl;
  auto PLCP = PLCPLCP.first;
  auto LCP = PLCPLCP.second;

  sdsl::int_vector<1> B(PLCP.size()+1, 0, 1);
  B[0] = 1;
  B[PLCP.size()] = 1;
  for(int i = 0; i < PLCP.size(); i++){
    if(LCP[i] < k){
      B[i] = 1;
    }else{
      //B[i] = 0;
    }
  }
  return B;
}

Interval suffixBinary(vector<pair<int,int>> SA, string text, string p, int d, vector<int64_t> C){
  Interval ret;
  int lowerBound = C[p.at(0)];
  int upperBound = lowerBound;
  int cmp = text.substr(SA[lowerBound].first,d).compare(p);
  cout << "binary...";
  while(cmp != 0){
    if(cmp < 0){
      lowerBound = (SA.size()-1 < lowerBound * 2)? SA.size()-1 : lowerBound*2;
    }else if(cmp > 0){
      //upperBound = lowerBound;
      lowerBound = lowerBound - floor(lowerBound/4)-1;
    }
    cmp = text.substr(SA[lowerBound].first,d).compare(p);
  }
  cout << "done";
  cmp = text.substr(SA[lowerBound].first,d).compare(p);
  while(cmp == 0){
    lowerBound--;
    cmp = text.substr(SA[lowerBound-1].first,d).compare(p);
  }
  upperBound = lowerBound;
  cmp = text.substr(SA[upperBound+1].first,d).compare(p);
  while(cmp == 0){
    upperBound++;
    cmp = text.substr(SA[upperBound+1].first,d).compare(p);
  }
  cout << "...done" << p << endl;
  return Interval(lowerBound, upperBound);
}

vector<Interval_pair> minimizerToBWTIntervalV2(vector<pair<string,int>> mini,int minimumDepth, vector<pair<int,int>> SA, vector<pair<int,int>> SAr, std::vector<int64_t> gc_arr, string text = ""){
  vector<Interval_pair> P;
  auto C = gc_arr;
  unordered_set<int> stored;
  string text2 = string(text.rbegin(), text.rend());
  for(int j = 0; j < mini.size(); j++){
    //    cout << "C[j] = " << C[mini[j].first.at(0)] << ", " << mini[j].first.at(0) << endl;
    /* This can be optimized */
    string revmini = string(mini[j].first.rbegin(), mini[j].first.rend());
    // Interval i1 = suffixBinary(SA, text, mini[j].first, minimumDepth, C);
    //Interval i2 = suffixBinary(SAr, text2, revmini, minimumDepth, C);
    auto a = C[mini[j].first.at(0)];
    if(a > 0){
      a--;
    }
    while(text.substr(SA[a].first,minimumDepth).compare(mini[j].first) != 0 && a < text.size()-minimumDepth){
      a++;
    }
    auto b = a;
    while(text.substr(SA[b].first,minimumDepth).compare(mini[j].first) == 0 && b < text.size()-minimumDepth){
      b++;
    }
    auto c = C[revmini.at(0)];
    if(c > 0){
      c--;
    }
    while(text2.substr(SAr[c].first,minimumDepth).compare(revmini) != 0 && c < text2.size()-minimumDepth){
      c++;
    }

    auto d = c;
    while(text2.substr(SAr[d].first,minimumDepth).compare(revmini) == 0 && d < text2.size()-minimumDepth){
      d++;
    }
    if((P.size() == 0 || stored.count(a) == 0)) {
      P.push_back(Interval_pair(a,b,c,d));
      //P.push_back(Interval_pair(i1,i2));
      //      cout << Interval_pair(a,b-1,c,d-1).toString() << endl;
      stored.insert(a);
    }
  }
  cout << "return p";
  return P;
}
vector<pair<Interval_pair,string>> minimizerToBWTInterval(sdsl::int_vector<1> bv, sdsl::int_vector<1> bvr, vector<pair<string,int>> mini, vector<pair<int,int>> SA, vector<pair<int,int>> SAr,BD_BWT_index<> index, pair<vector<int>,vector<int>> lcp1, pair<vector<int>,vector<int>> lcp2, string text = ""){
  // cout << "Enter minimizerToBWTInterval" << endl;
  sdsl::rank_support_v<1> b_rank(&bv);
  sdsl::rank_support_v<1> b_rankr(&bvr);
  //Returns the amount of ones
  size_t ones = sdsl::rank_support_v<1>(&bv)(bv.size());
  size_t onesr = sdsl::rank_support_v<1>(&bvr)(bvr.size());
  sdsl::bit_vector::select_1_type b_select(&bv);
  sdsl::bit_vector::select_1_type b_selectr(&bvr);
  // for(auto a : bv){
  //   cout << a << ",";
  // }cout << endl;
  // for(auto a : bvr){
  //   cout << a << ",";
  // }
  // cout << endl;
  //     cout << endl;
  //rank_support_v<1> returns the amount of 1's up to, but NOT including x in b_rank(x), count(x \in [0...x[ )
  //cout << b_rank(18) << endl;
  //bit_vector::select_0_type returns the proper index of x:th 0 in given vector.
  //cout << b_select(0+1) << endl;
  cout << endl;
  // for(int i = 0; i < SAr.size(); i++){
  //   cout << SA[i].first << "\t" << SA[i].second << "\t" << bv[i] << "\t" << lcp1.second[i] << "\t"<<i<<"\t" <<SAr[i].first << "\t" << SAr[i].second << "\t" << bvr[i] << "\t" << lcp2.second[i] << endl;
  // }
  // cout << "ones: " << ones << endl;
  // cout << "onesr: " << onesr << endl << endl;
  vector<pair<Interval_pair,string>> P;
  //auto C = index.get_global_c_array();
  unordered_set<int> stored;
  for(int j = 0; j < mini.size(); j++){
    //    cout << mini[j].second << ", " << mini[j].first << endl;
    if(mini[j].second + mini[j].first.length() > text.length()) continue;
    int i = mini[j].second; //mini[j].second denotes the index of the k-mer on the original k-mer listing (all k-mers of text).
    int ir = (text.length()-1)-i;
    //i = (i < bv.size()-1)? i : bv.size()-1;
    //    i = (i > 1)? i : 1;
    // cout << " i = " << i << endl;
    // auto s  = SA[i].second;
    // auto s2 = SAr[i].second;
    // if(i > SA.size()-1){
    //   i = 0;
    // }
    // if(i == 0){
    //   i = SA.size()-1;
    // }

    int revOffset (mini[j].first.length()-1);
    auto s  = SA[i].second;
    auto s2 = SAr[ir-revOffset].second;
    //cout << "bv  at s(" <<s<<") = " << bv[s] << endl;
    //cout << "bvr at s2("<<s2<<")= " << bvr[s2] << endl;
    //auto s  = i;
    //auto s2 = ir-revOffset;

    string text2 = string(text.rbegin(), text.rend());
    //cout << "Minimizer = " << mini[j].first << endl;
    //cout << "Minimizer from text fwd("<<i<<"): "  << text.substr(i, mini[j].first.length()) << endl;
    //cout << "Minimizer from text rev("<<ir<<"): " << text2.substr(ir-revOffset, mini[j].first.length()) << endl;
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

    //cout << "b_rank(SA[i].second) r1,r2= " << r1 << "," << r2 << "\t s = " << s << endl;
    //cout << "b_rank(SA[i].second) r3,r4= " << r3 << "," << r4 << "\t s2= " << s2 << endl;
    //
    // if(r1 < ones){
    //   r1 = r1;
    // }
    // if(r3 < ones){
    //   r3 = r3;
    // }
     int a,b,c,d;
    if(r1 == 0){
      a = r1;
    }else{
     	a = r1;
    }
    b = b_select(r2+1)-1;
    if(bv[b] == 1) //b--;
    if(r3 == 0){
      c = r3;
    }else{
      c = r3;
    }
    d = b_selectr(r4+1)-1;
    if(bvr[d] == 1) //d--;
    //  cout << "Minimizer from select fwd("<<a<<","<<b<<"): "  << text.substr(SA[a].first) << endl;
    //cout << "Minimizer from select rev("<<c<<","<<d<<"): "  << text2.substr(SAr[c].first) << endl;
    // cout << "interval: [" << a << ", " << b  <<"]" << ",[" << c << "," << d << "] >>" << mini[j].first << endl;

    // for(int i = a; i <= b; i++){
    //   cout << "\t SA[" << i << "]" << text.substr(SA[i].first,3) <<endl;
    // }
    // for(int i = c; i <= d; i++){
    //   cout << "\t SAr[" << i <<"]" << text2.substr(SAr[i].first,3) << endl;
    // }
   // if((P.size() == 0) || true) {
      // if(b > SA.size()-1){
      // 	b = b-SA.size()-1;
      // 	a = a-SA.size()-1;
      // }
      // if(d > SAr.size()-1){
      // 	d = d-SAr.size()-1;
      // 	c = c-SAr.size()-1;
      // }
      //      cout << " pushed" << endl;
      P.push_back(make_pair(Interval_pair(a,b,c,d),mini[j].first));
      //cout << "pushed: " << P.back().first.toString() << endl;
      stored.insert(a);
   // }
  }
  return P;
}

#endif
