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
pair< vector<tuple<int,int,int>> ,pair< vector<pair<string,int>> , vector<pair<string,int>> > > minimizerTuples(vector<pair<string,int>> m1, vector<pair<string,int>> m2);
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
vector<Interval_pair> minimizerToBWTInterval(sdsl::int_vector<1> bv, sdsl::int_vector<1> bvr, vector<pair<string,int>> mini, vector<pair<int,int>> SA, vector<pair<int,int>> SAr,BD_BWT_index<> index, string text);


struct mimsort {
  bool operator() (const pair<string,int>& first, const pair<string,int>& second) const{
    return first.first < second.first;
  }
};


bool mimCompare (pair<string,int> first, pair<string,int> second){
  return first.first < second.first;
}


vector<pair<string,int>> minimizers(string t1, int k, int w){
  //  cout << "enter minimizer";
  vector<set<pair<string,int>,mimsort>> kmers((t1.size()));
  vector<pair<string,int>> ret;

  for(int i = 0; i < t1.size()-k; i++){
    //inserting k-mers into w buckets, which are kept in lexicographic order.
    kmers.at((int)(i/w)).insert(make_pair(t1.substr(i,k),i));
    //    cout << ".";
  }
  for(auto km : kmers){
    if(!km.empty()){
      auto x = *km.begin();
      //    cout << x.first << endl;
      ret.push_back(x);
    }
  }
  cout << "found " << ret.size() << " minimizers" << endl;
  sort(ret.begin(), ret.end(), minimizerLexSort);
  return ret;
}

pair<vector<pair<string,int>>,vector<pair<string,int>>> mutualMinimizers(vector<pair<string,int>> m1, vector<pair<string,int>> m2){
  unordered_set<string> m1Kmer(m1.size()-1);
  unordered_set<string> m2Kmer(m2.size()-1);
  set<string> anchorMers;
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
  // for(auto m: anchorMers){
  //   cout << "anchorMer: "<< m << endl;
  // }
  int maxSize = (m1.size()-1 > m2.size()-1)? m1.size()-1 : m2.size()-1;
  for(int i = 0; i < maxSize; i++){
    if((m1.size()-1 > maxSize) && (anchorMers.find(m1.at(i).first) == anchorMers.end())){
      m1.erase(m1.begin()+i);
    }
    if((m2.size()-1 > maxSize) && (anchorMers.find(m2.at(i).first) == anchorMers.end())){
      m2.erase(m2.begin()+i);
    }
  }
  return make_pair(m1,m2);
}
pair< vector<tuple<int,int,int>> ,pair< vector<pair<string,int>> , vector<pair<string,int>> > > minimizerTuples(vector<pair<string,int>> m1, vector<pair<string,int>> m2){
  //  pair<vector<tuple<int,int,int>>,pair<vector<pair<string,int>>,vector<pair<string,int>>>> ret;
  vector<tuple<int,int,int>> retTuple;
  pair<vector<pair<string,int>>,vector<pair<string,int>>> retRaw;
  
  auto muts = mutualMinimizers(m1,m2);
  m1 = muts.first;
  m2 = muts.second;
  sort(m1.begin(), m1.end(), mimCompare);
  sort(m2.begin(), m2.end(), mimCompare);
  //cout << "sorted minimems" << endl;
  for(auto x : m1){
    int i = 0;
    for(auto y : m2){
      i++;
      //cout << "Comparing: " << x.first << " & " << y.first;
      if(x.first.compare(y.first) == 0){
	m2.erase(m2.begin(), m2.begin()+i);
	i = 0;
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
  sort(retTuple.begin(),retTuple.begin(),memSort);
  return make_pair(retTuple,retRaw);
}
vector<tuple<int,int,int>> memifyMinimizers(vector<tuple<int,int,int>> mini, string text, string text2){
  set<tuple<int,int,int>> miniMemsSet;
  vector<tuple<int,int,int>> miniMems;
  auto last = mini.at(0);
  int a,b,c;
  tie(a,b,c) = mini.at(0);
  for(int x = 1; x < mini.size(); x++){
    int i,j,k;
    tie(i,j,k) = mini.at(x);
    if(i-a == j-b && a+c >= i && b+c >= j){
      mini.erase(mini.begin()+x);
      mini.at(x-1) = make_tuple(a,b,c+abs(c-k));
    }
    a = i;
    b = j;
    c = k;
  }
  for(auto m : mini){
    int i,j,k;
    tie(i,j,k) = m;
    int i2,j2,k2;
    //while(i > 0 && j > 0 && i < text.length() && j < text2.length() && text.at(i-1) == text2.at(j-1)){
    i2 = i;
    j2 = j;
    k2 = k;
    while(i2-k2 >= 0 && j2-k2 >= 0 && text.substr(i2,k).compare(text2.substr(j2,k)) == 0){
      i2 = i2-k;
      j2 = j2-k;
      k2 = k2+k;
    }
    while(text.substr(i2,k2) != text2.substr(j2,k2) && (i2+k2) < text.length() && (j2+k2) < text2.length() && k2 > 0){
      //      cout << "minimem: (" << i2 << "," << j2 << "," << k2 << ") "<< endl;
      i2 = i2+1;
      j2 = j2+1;
      k2 = k2-1;
    }
    //cout << "second done" << endl;
    i = i2; j = j2; k = k2;
    //}
    while(i+k < text.length() && j+k < text2.length() && text.at(i+k) == text2.at(j+k)){
      k++;
    }
    miniMemsSet.insert(make_tuple(i,j,k));
    tie(a,b,c) = m;
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
pair<vector<int>,vector<int>> createPLCP(BD_BWT_index<> index, int q, string t, bool lcp, vector<pair<int,int>> SA, bool direction){
  vector<int> phi(SA.size());
  vector<int> PLCP(SA.size());
  //cout << "SA size: " << SA.size() << endl;
  auto LF = mapLF(index, direction);
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
    }
  }
  return B;
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
    if((P.size() == 0) || stored.count(a) == 0) {
      P.push_back(Interval_pair(a,b,c,d));
      //      cout << Interval_pair(a,b-1,c,d-1).toString() << endl;
      stored.insert(a);
    }
  }
  cout << "return p";
  return P;
}
  
vector<Interval_pair> minimizerToBWTInterval(sdsl::int_vector<1> bv, sdsl::int_vector<1> bvr, vector<pair<string,int>> mini, vector<pair<int,int>> SA, vector<pair<int,int>> SAr,BD_BWT_index<> index, string text = ""){
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
  // for(int i = 0; i < SAr.size(); i++){
  //   cout << SA[i].first << "\t" << SA[i].second << "\t\t\t" <<SAr[i].first << "\t" << SAr[i].second << endl;
  // }
  // cout << "ones: " << ones << endl;
  // cout << "onesr: " << onesr << endl << endl;
  vector<Interval_pair> P;
  auto C = index.get_global_c_array();
  unordered_set<int> stored;
  for(int j = 0; j < mini.size(); j++){
    //    cout << mini[j].second << ", " << mini[j].first << endl;
    int i = mini[j].second; //mini[j].second denotes the index of the k-mer on the original k-mer listing (all k-mers of text).
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
    //auto s  = SA[i].second;
    //auto s2 = SAr[i].second;

    // auto s  = i;
    // auto s2 = i;
    string text2 = string(text.rbegin(), text.rend());
    //    cout << "Minimizer from text: "     << text.substr((SA[SA[i].second].first)) << endl;
    //cout << "Minimizer from text rev: " << text2.substr((SAr[SAr[i].second].first)) << endl;
    //int r1 = b_rank(s);
    //int r3 = b_rankr(s2);
    // cout << "r1 = " << r1 << "\t r3 = " << r3 << endl;
    // cout << "SA[i].first = " << SA[i].first << "\t SA^r[i].first = " << SAr[i].first << endl;
    // cout << "SA[i].second = " << SA[i].second << "\t SA^r[i].second = " << SAr[i].second << endl;
    // cout << "b_rank(SA[i].second) r1= " << r1 << endl;
    // cout << "b_rank(SA[i].second) r2= " << r2 << endl;
    // cout << "b_rank(SA[i].second) r3= " << r3 << endl;
    // cout << "b_rank(SA[i].second) r4= " << r4 << endl; 
    // cout <<"i: "<< i << endl;
											    
    // if(r1 < ones){
    //   r1 = r1;
    // }
    // if(r3 < ones){
    //   r3 = r3;
    // }
    // auto a = b_select(r1);
    // auto b = b_select(r1+1)-1;
    // auto c = b_selectr(r3);
    // auto d = b_selectr(r3+1)-1;

    /* This can be optimized */
    string revmini = string(mini[j].first.rbegin(), mini[j].first.rend());
    auto a = C[mini[j].first.at(0)];
    while(text.substr(SA[a].first,3).compare(mini[j].first) != 0){
      a++;
    }
    auto b = a;
    while(text.substr(SA[b+1].first,3).compare(mini[j].first) == 0){
      b++;
    }
    auto c = C[revmini.at(0)];
    while(text2.substr(SAr[c].first,3).compare(revmini) != 0){
      c++;
    }
    auto d = c;
    while(text2.substr(SAr[d+1].first,3).compare(revmini) == 0){
      d++;
    }

    // cout << "interval: [" << a << ", " << b  <<"]" << ",[" << c << "," << d << "] >>" << mini[j].first << endl;
    
    // for(int i = a; i <= b; i++){
    //   cout << "\t SA[" << i << "]" << text.substr(SA[i].first,3) <<endl;
    // }
    // for(int i = c; i <= d; i++){
    //   cout << "\t SAr[" << i <<"]" << text2.substr(SAr[i].first,3) << endl;
    // }
    if((P.size() == 0) || true) {
      // if(b > SA.size()-1){
      // 	b = b-SA.size()-1;
      // 	a = a-SA.size()-1;
      // }
      // if(d > SAr.size()-1){
      // 	d = d-SAr.size()-1;
      // 	c = c-SAr.size()-1;
      // }
      //      cout << " pushed" << endl;
      P.push_back(Interval_pair(a,b,c,d));
      stored.insert(a);
    }
  }
  return P;
}

#endif
