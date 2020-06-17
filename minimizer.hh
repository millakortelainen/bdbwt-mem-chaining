#include <iostream>
#include <vector>
#include <set>
#include <unordered_set>
#include <algorithm>
#include "util.hh"
using namespace std;
struct mimsort {
  bool operator() (const pair<string,int>& first, const pair<string,int>& second) const{
    return first.first < second.first;
  }
};

vector<pair<string,int>> minimizers(string t1, int k, int w){
  vector<set<pair<string,int>,mimsort>> kmers((t1.size()/w));
  vector<pair<string,int>> ret;
  for(int i = 0; i <= t1.size()-k; i++){
    kmers.at((int)(i/w)).insert(make_pair(t1.substr(i,k),i));
  }
  for(auto km : kmers){
    if(!km.empty()){
      auto x = *km.begin();
      //cout << x.first << endl;
      ret.push_back(x);
    }
  }
  cout << "found " << ret.size() << " minimizers" << endl;
  sort(ret.begin(), ret.end(), minimizerLexSort);
  return ret;
}

vector<tuple<int,int,int>> minimizerTuples(vector<pair<string,int>> m1, vector<pair<string,int>> m2){
  unordered_set<string> m1Kmer(m1.size()-1);
  unordered_set<string> m2Kmer(m2.size()-1);
  set<string> anchorMers;
  vector<tuple<int,int,int>> ret;
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
  // for(auto j : anchorMers){
  //   cout <<"anchor: "<< j << endl;
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
  int i = 0;
  for(auto x : m1){
    for(auto y : m2){
      if(x.first == y.first){
	auto tup = make_tuple(x.second, y.second, x.first.length());
	//	cout << "mini tuple (" << x.second << "," <<  y.second << "," << x.first.length() << ")" << endl;
	ret.push_back(tup);
      }
    }
  }
  sort(ret.begin(),ret.begin(),memSort);
  return ret;
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
    
    while(i > 0 && j > 0 && i < text.length() && j < text2.length() && text.at(i-1) == text2.at(j-1)){
      i--;
      j--;
      k++;
    }
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
