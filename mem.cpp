#include <iostream>
#include <bitset>
#include "include/BD_BWT_index.hh"
#include "include/Iterators.hh"
#include "rsa1d.hh"
#include "rsa2d.hh"
#include <string>
#include <set>
#include <stdio.h>
#include <inttypes.h>
#include <tuple>
#include <stack>

using namespace std;

bool verboseEnumeration = false;
bool verboseSigma = false;
bool verboseI = false;
bool verboseMaximum = false;
bool verboseSubroutine = false;
bool verboseElement = false;
set<char> alphabet;

struct suffix{
  int index;
  string suffix;
};
int suffcmp(struct suffix first, struct suffix last){
  return (first.suffix < last.suffix)? 1:0;
}
int* build_suffix_array(string text){
  text += BD_BWT_index<>::END;
  struct suffix suffix_array[text.size()+1];
  int *suffix_index_array = new int[text.size()+1];

  for(int i = 0; i < text.size(); i++){
    suffix_array[i].index = i;
    suffix_array[i].suffix = text.substr(i,text.size()+1);
  }
  sort(suffix_array,suffix_array+text.size(),suffcmp);

  for(int i = 0; i < text.size(); i++){
    cout << suffix_array[i].index << "\t" << suffix_array[i].suffix  << "\n";
    suffix_index_array[i] = suffix_array[i].index;
  }
  return suffix_index_array;
}

/** Map LF values for the given index.
    @return pair<map<int,int>,int> such that .first is the actual mapping, and .second is the zeroth index so that we don't need to search for it again.
*/
pair<map<int,int>,int> mapLF(BD_BWT_index<>& index){
  map<char,int> m;
  map<int,int> lfmapping;
  auto C = index.get_global_c_array();
  int zeroth_index = -1;
  for(auto i : alphabet){ 
    m[i] = 0;
  }
  for(int i = 0; i < index.size(); i++){
    char l = char(index.forward_bwt_at(i));
    lfmapping[i] = C[l]+m[l];
    if(lfmapping[i] == 0){
      zeroth_index = i;
    }
    m[l] += 1;    
    
  }
  return make_pair(lfmapping,zeroth_index);
}
/** Print out BWT, LF and SA indices.
param index BD_BWT_index<>& Burrows Wheeler index.
param text1 string original text that the BWT is based on.
*/
void pretty_print_all(BD_BWT_index<>& index, string text1){
  string separator = "+------------------------------------------------+\n";
  cout << separator;
  auto SA = build_suffix_array(text1);
  cout << separator;
  cout << "Text: " << text1 << "\n";
  cout << separator;
  cout << "| bwd\t" << "ln\t" << "fwd\t" << "LF\t" << "SA[i]\t" << "\t |\n";
  auto mapping = mapLF(index).first;
  for(int i = 0; i < index.size(); i++){
    char t = (char)index.forward_bwt_at(i);
    char tr= (char)index.backward_bwt_at(i);
    cout << "| "<< tr << "\t("<<i<<")\t" << t << "\t" << mapping[i] << "\t" << SA[i] << "\t\t |\n";
  }
  cout << separator;
}
/** Enumerates the unique characters on the forward index of the BWT.
    param idx BD_BWT_index<> Burrows-wheeler transform
    param ip Interval_pair range to find the unique characters in.
    return std::vector<uint8_t> array containing each unique character.
*/
std::vector<uint8_t> enumerateLeft(BD_BWT_index<> idx, Interval_pair ip){
  set<uint8_t> s;
  vector<uint8_t> ret;
  for(int i = ip.forward.left; i <= ip.forward.right; i++){
    auto c = idx.forward_bwt_at(i);
    s.insert(c); //Set is always sorted
    //    cout << "enumerated left symbol: " << c << ", from forward interval of: " << ip.forward.toString() <<"\n";
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
  set<uint8_t> s;
  vector<uint8_t> ret;
  for(int i = ip.reverse.left; i <= ip.reverse.right; i++){
    auto c = idx.backward_bwt_at(i);
    s.insert(c); //Set is always sorted
    //cout << "enumerated right symbol: " << c  << ", from backward interval of: "<< ip.reverse.toString() << "\n";
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
  for(auto i : A){
    a = i.first;
    b = i.second;
    for(auto j: B){
      c = j.first;
      d = j.second;
      if(a != c && b != d){ // Pseudocode gives the definition of C = ((a,b,c,d) | (a,b) \in A, (c,d) \in B, a != c, b != d). This however is leading to inconsistent behaviour on returning MEM's.

	cross.push_back(make_tuple(a,b,c,d));
	if(verboseSubroutine){
	  if(a == BD_BWT_index<>::END) a = '$';
	  if(b == BD_BWT_index<>::END) b = '$';
	  if(c == BD_BWT_index<>::END) c = '$';
	  if(d == BD_BWT_index<>::END) d = '$';
	}
      }
    }
    if(verboseSubroutine) cout << "---" << "\n";
  }
  
  for(int k = 0; k < cross.size(); k++){
    tie(a,b,c,d) = cross[k];
    Interval_pair i_1 = idxS.left_extend(pr.first,a);  //if(i_1.reverse.size() == 0) i_1 = pr.first;
    if(verboseSubroutine) cout << "Extended interval " << pr.first.toString() << " to the left and got: " << i_1.toString() << "\n";
    
    Interval_pair i_2 = idxS.right_extend(i_1,b); //if(i_2.reverse.size() == 0) i_2 = i_1; 
    if(verboseSubroutine) cout << "Extended interval " << i_1.toString() << " to the right and got: " << i_2.toString() << "\n";
    
    Interval_pair i_3 = idxT.left_extend(pr.second,c); //if(i_3.reverse.size() == 0) i_3 = pr.second;
    if(verboseSubroutine) cout << "Extended interval " << pr.second.toString() << " to the left and got: " << i_3.toString() << "\n";
    
    Interval_pair i_4 = idxT.right_extend(i_3,d); //if(i_4.reverse.size() == 0) i_4 = i_3;
    if(verboseSubroutine) cout << "Extended interval " << i_3.toString() << " to the right and got: " << i_4.toString() << "\n";

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

  //  auto lftext1 = mapLF(idxS).second;
  //  auto lftext2 = mapLF(idxT).second;
  // cout << "size " << itrl1Size << "\n";
  S.push(make_tuple(Interval_pair(0,itrl1Size,0,itrl1Size), Interval_pair(0,itrl2Size,0,itrl2Size),0));
  
  while(!S.empty()){
    if(verboseElement) cout << "\n\n---Popping new element:" << "";
    tie(ip0,ip1,depth) = S.top(); S.pop();
    if(verboseElement) cout << ip0.toString() << ip1.toString() << ", current depth is: " << depth << " \nStack has: " << S.size() << " elements left."  << "\n\n";
    
    if((ip0.forward.right - ip0.forward.left+1) < 1 || (ip1.forward.right - ip1.forward.left+1) < 1){
      if(verboseElement) cout << "Invalid Element" << ip0.toString() << ip1.toString() << "\n";
      //continue;
    }
    
    if(idxS.is_left_maximal(ip0) || idxT.is_left_maximal(ip1) || (enumerateLeft(idxS,ip0) != enumerateLeft(idxT,ip1))){ //had to take forward and reverse unlike pseudo
      if(depth > 0){
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
      
      if(verboseSigma){
	if(cp == BD_BWT_index<>::END){
	  cp = '$';
	}
	cout << "Sigma: " << cp << "\n";
	cout << "\tExtending interval left on IdxS" << ip0.toString() << " with character: " << cp << " ... -> " << i1.toString()<<"\n";
	cout << "\tExtending interval left on IdxT" << ip1.toString() << " with character: " << cp << " ... -> " << i2.toString()<<"\n";
      }
      if(i1.forward.left < 0 || i2.forward.left < 0){ //left_extend() returns [-1,-2] interval if extending with character c is not possible.
	continue; //We don't need to bother with invalid index.
      }
      I.insert(make_pair(i1,i2));
      if(verboseI) cout << "\t\tI size: " << I.size() << " Added: " << i1.toString() << i2.toString() << "\n";
    }
    
    if(verboseSigma) cout << "---\n";
    pair<Interval_pair,Interval_pair> x;
    if(I.size() != 0){
      x = *I.begin(); //Initialize first item for comparison
    }else{
      continue;
    }

    int xForwardDelta = x.first.forward.right - x.first.forward.left;
    int xForwardDelta2 = x.second.forward.right - x.second.forward.left;
    int maxDelta = xForwardDelta + xForwardDelta2;
    if(verboseMaximum) cout << "\n";
    
    for(auto y : I){
      if(y.first.forward.left < 0 || y.second.forward.left < 0 || y.first.reverse.left < 0 || y.second.reverse.right < 0){
	continue;
      }
      if(verboseMaximum) cout << "Comparing x:\t " << x.first.toString() << "\t" << x.second.toString() << "\nagainst y:\t " << y.first.toString() << "\t" << y.second.toString() << "\n";

      int yForwardDelta = y.first.forward.right - y.first.forward.left;
      int yForwardDelta2 = y.second.forward.right - y.second.forward.left;

      int zDelta 	= yForwardDelta + yForwardDelta2;
      
      if(zDelta > maxDelta){
	x = y;
	if(verboseMaximum) cout << "Found greater:\t " << y.first.toString() << "\t" << y.second.toString() << "\n";
	maxDelta = zDelta;
      }
    }
    if(verboseMaximum) cout << "\nMaximum interval computed was: x = " << x.first.toString() << x.second.toString() << " Proceeding to remove from I..."<< "\n";
    ip0 = x.first;
    ip1 = x.second;
    
    if(I.size() == 0){
      if(verboseI) cout << "\t" << "No values in I" << "\n";
    }else{
      I.erase(x);
      if(verboseI) cout << "Removed x from I, I is now: " << "\n";
      if(I.size() == 0){
	if(verboseI) cout << "\t" << "No values in I" << "\n";
      }else{
	for(auto i : I){
	  if(verboseI) cout << "\t" << i.first.toString() << i.second.toString() << "\n";	  
	}
	for(auto y : I){
	  if(idxS.is_right_maximal(ip0) || idxT.is_right_maximal(ip1) || (enumerateRight(idxS, ip0) != enumerateRight(idxT, ip1))){
	    if(verboseElement) cout << "Pushed y:" << y.first.toString() << y.second.toString() << "\n";
	    S.push(make_tuple(y.first,y.second, depth+1));
	  }
	}
      }
    }
    //cout << "Interval x: " <<ip0.toString() << ip1.toString() << "rmax1: " << idxS.is_right_maximal(ip0) << "rmax2: " << idxS.is_right_maximal(ip1) << "\n";
    // for(auto i : enumerateRight(idxS, ip0)){
    //   cout << "enum1: " << i << "\n";
    // }
    // for(auto i : enumerateRight(idxS, ip0)){
    //   cout << "enum2: " << i << "\n";
    // }
    if(idxS.is_right_maximal(ip0) || idxT.is_right_maximal(ip1) || (enumerateRight(idxS, ip0) != enumerateRight(idxT, ip1))){ //had to take forward and reverse unlike pseudo
      if(verboseElement) cout << "Pushed x:" << x.first.toString() << x.second.toString() << "\n";
      S.push(make_tuple(x.first,x.second, depth+1));
    }
  }
  return ret;
}

bool pairSort(const pair<int,int> &first, const pair<int,int> &second){
  return (first.first < second.first);
}
bool pairSort2(const pair<int,tuple<int,int,int>> &first, const pair<int,tuple<int,int,int>> &second){
  return (first.first < second.first);
}
int maxFour(int &a, int &b, int &c, int &d){
  int e = (a > b)? a:b;
  int f = (c > d)? c:d;
  return (e > f)? e:f;
}

vector<int> chaining(vector<Interval_pair> A, int size){
  rsa1d T_a = rsa1d(size); 
  rsa1d T_b = rsa1d(size);
  rsa2d T_c = rsa2d(size);
  rsa2d T_d = rsa2d(size);
  vector<int> C_a(A.size(),-1);
  vector<int> C_b(A.size(),-1);
  vector<int> C_c(A.size(),-1);
  vector<int> C_d(A.size(),-1);
  vector<int> C_p(A.size(),-1);
  vector<int> C(A.size(),-1);
  vector<pair<int,int> > E_1;
  vector<pair<int,int> > E_2;
  T_a.update(0,-999);
  T_b.update(0,-999);
  for(int j = 0; j < A.size(); j++){
    T_a.update(A[j].reverse.right, -999);
    T_b.update(A[j].reverse.right, -999);
    if((A[j].reverse.left - A[j].forward.left) >= 0){
      T_c.update(A[j].reverse.left - A[j].forward.left,A[j].forward.right,-999);
      T_d.update(A[j].reverse.left - A[j].forward.left,A[j].reverse.right,-999);
      // cout << A[j].forward.left << "\t" << A[j].reverse.left <<"\t" << A[j].reverse.left - A[j].forward.left << "\n";
    }
    auto p1 = make_pair(A[j].forward.left, j);
    auto p2 = make_pair(A[j].forward.right, j);
    // cout << p2.first << "\t" << p2.second << "\n";
    E_1.push_back(p1);
    E_2.push_back(p2);

  }
  for(auto k : E_2){
    E_1.push_back(k);
  }
  sort(E_1.begin(), E_1.end(), pairSort);

  for(int i = 0; i < E_1.size(); i++){
    int j = E_1[i].second;
    Interval_pair I = A[j];

    if(I.forward.left == E_1[i].first){
      // cout << "round: " << j << "\n";
      C_a[j] = T_a.rangeMax(0, I.reverse.left-1).second;
      C_b[j] = I.reverse.left + T_b.rangeMax(I.reverse.left,I.reverse.right).second; 
      C_c[j] = I.forward.left + T_c.rangeMax(0,I.reverse.left - I.forward.left,0, I.forward.right).second;
      C_d[j] = I.reverse.left + T_d.rangeMax((I.reverse.left - I.forward.left)+1,999,0, I.reverse.right).second;

      // cout << "error?: " <<  T_d.rangeMax((I.reverse.left - I.forward.left),999,0, I.reverse.right).second << "\n";

      // cout << C_a[j] << "\t" << C_b[j] << "\t" <<  C_c[j] << "\t"  <<  C_d[j] << "\n" ;

      C[j] = maxFour(C_a[j], C_b[j],C_c[j],C_d[j]);
      C_p[j] = C[j]+I.forward.right-I.forward.left+1;
      
      // cout << (int)C[j]-(int)I.forward.left << "\n";

      
      T_c.upgrade(I.reverse.left - I.forward.left, I.forward.right, (int)C[j]-(int)I.forward.left);
      T_d.upgrade(I.reverse.left - I.forward.left, I.forward.right, (int)C[j]-(int)I.reverse.left);
    }else{
      T_a.upgrade(I.reverse.right,C_p[j]);
      T_b.upgrade(I.reverse.right,C[j]- I.reverse.left);
      T_c.update(I.reverse.left - I.forward.left, I.forward.right, -999);
      T_d.update(I.reverse.left - I.forward.left, I.forward.right, -999);
    }
  }
  return C_p;
}
/** Creating SA with use of recursive LF mapping.
*/
vector<int> int_ret_recurse(BD_BWT_index<> idxS, map<int,int> LFI, int i, int currIndex, int k, vector<int> retSA){
  if(k < idxS.size()-1){
    retSA[currIndex] = idxS.size()-1 - k;
    return int_ret_recurse(idxS, LFI,i,LFI[currIndex], k+1, retSA);
  }else{
    retSA[currIndex] = idxS.size()-k-1;
    return retSA;
  }
  
}
/** Creating SA with use of recursive LF mapping.
*/
vector<int> returnIntervals(BD_BWT_index<> idxS, BD_BWT_index<> idxT, int i){
  auto lfS = mapLF(idxS);
  vector<int> retSA(idxS.size(),-9);
  retSA[lfS.second] = idxS.size()-1;
  retSA = int_ret_recurse(idxS, lfS.first, i, lfS.first[lfS.second],  0, retSA);

  return retSA;
    
}
/** TODO RadixSort
 */
void radixSort(vector<pair<int,int>> list, int r){
  sort(list.begin(), list.end(), pairSort);
}
/** TODO RadixSort
 */
void radixSort(vector<pair<int,tuple<int,int,int>>> list, int r){
  sort(list.begin(), list.end(), pairSort2);
}
/** TODO batch locate
 */	       
vector<pair<int,tuple<int,int,int>>> batchLocate(vector<pair<int,tuple<int,int,int>>>  pairs, vector<bool> marked, BD_BWT_index<> bwt){
  int i = 1;
  int n = bwt.size()-1;
  vector<pair<int,int>> translate;
  auto LFindex = mapLF(bwt).first;
  
  for(int j = n; j > 1; j--){
    if (marked[i]){
      translate.push_back(make_pair(i,j));
    }
    i = LFindex[i];
  }
  radixSort(translate, 1);
  radixSort(pairs, 1);
  i = 1;
  int j = 1;

  while(i <= pairs.size()){
    if(pairs[i].first == translate[j].first){
      pairs[i].first = translate[j].second;
      i++;
    }else{
      j++;
    }
  }
  return pairs;
}

int main(){
  string text =  "GTGCGTGATCATCATTT";
  string text2 = "GTGCAAAGTGATTACCA";
  //string text2 = "ASDKASSATTKATI";

  //500 bases.
  //string text = "gttcaccatttaaataatcttcaatatcaacacgcgaagctcgcttgcagggatgaactgaatagacctgtttactccggaaaagcaagactatcctggtgctgatgctacggtacattgttcttggcacgattacggactattcacactgaatccgggtggggagggccttatggacacgtaatatgcgcgtactggttggcgttgtagacgcgcaacttcatcgataatctgactgcctgacaagctaccagcaatacgttactccatcccgctatcctcggtactgcttgcggtgtcaccccgttaagtgacgtcctgttcgcggctaggctacgagttgcgttaatgcactctgaatcagaattccgcagcgttaagctggcttcaccagcgtcttcggtctgacttaaacctactcccgacatttctacagtgactactgtgtacgccccacgaagtcaaccccgagctacacctaaccggcctccagcactgcc";
  
  //string text2 = "aattcgaatagttagctgacgtacgacatgttaccttaataatataactggtgtccgcgactgagtgctctcctacctcccacgagcctcaggaaaaacgtctttaaatctctacccggagctgtttaaggggaagccaactcgaacctagcagggcattaaatttgtattgcaccaaaacgaccggcttaacattccgtgtctcactggacggaaaaccaacctaagcagtatttggcctcctggtaggcgaaccatctacggtggaccgtataatcggactaaccggcaggtttacacttcgcaatgctacgctgcccagggccgggcccccagtaggtttgcactgtagagggagggccggagtgtatcccccatcggtaactctacatatgcgcaagccgccctgggcaagatcccatcccactcgtgtggctctcgcgccgggtggattgtacgatcggaatcctctggggacgcgcgttcagtaacttcgctta";
  
  std::vector<Interval_pair> Ipairs;
  
  BD_BWT_index<> index((uint8_t*)text.c_str());
  BD_BWT_index<> index2((uint8_t*)text2.c_str());
  
  auto mems  = bwt_mem2(index, index2);
  set<tuple<int,int,int> > removedupes;
  for(auto a: mems ){
    removedupes.insert(a);
  }

  auto retSA = returnIntervals(index,index2,0);
  auto retSA2 = returnIntervals(index2,index,0);
  //int ind = 0;
  //for(auto r : retSA){
  //  cout << "retSA: " << ind<< ", "<< r << "\n";
  //  ind++;
  //}
  
  for(auto a : removedupes){
    int i, j, depth;
    tie(i,j,depth) = a;
    if(depth > 2 && retSA[i]+depth+1 < index.size() && retSA2[j]+depth+1 < index2.size()){
      // cout << "subroutine > i: " << i << " j: " << j << " depth: " << depth << "\n";
      int begin_i= retSA[i];
      int end_i  = begin_i+depth-1;
      int begin_j= retSA2[j];
      int end_j  = begin_j+depth-1;
      
      //Printing with the depth information, currently depth is wrong.
      cout << "Depth given in tuple (i,j,d): "<< depth << "\n";
      cout << "S: ["<< begin_i <<","<< end_i <<"]" << "-->\t";

      for(int b = begin_i+1; b <= end_i+1; b++){
      	cout << text[b];
      }
      cout << "\n";
      cout << "T: ["<< begin_j <<","<< end_j <<"]" << "-->\t";
      for(int b = begin_j+1; b <= end_j+1; b++){
      	cout << text2[b];
      }
      cout << "\n";

      //Printing by checking equality, not ideal but works for ensuring correctness until depth works properly.
      int d = 0;
      while(text[begin_i+d] == text2[begin_j+d] && begin_i+d < text.size() && begin_j+d < text2.size()){
	cout << text[begin_i+d];
	d++;
      }
      if(begin_j >= begin_i){
	Ipairs.push_back(Interval_pair(begin_i, begin_i+d,begin_j, begin_j+d));
      }
      cout << ", final depth: " << d<< "\n";
      
      cout << "\n";
    }
  }

  //pretty_print_all(index,text);
  //pretty_print_all(index2,text2);

  auto chains = chaining(Ipairs, text2.size());
  for(int i = 0; i < chains.size(); i++){
   cout <<"Chain["<< i << "]: " << chains[i] << "\n";
  }
}

