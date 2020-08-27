#ifndef BDBWT_MEM_CHAIN_RSA1D_HH
#define BDBWT_MEM_CHAIN_RSA1D_HH
#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>
#include <bits/stdc++.h>
using namespace std;
class rsa1d{
private:
  bool verbose = false;
  unordered_map<int,int> hash;
public:
  //<value, identifier>
  struct cell1d{
    pair<int,int> primary;
    int value;
  };
  static int cellcmp1d(struct cell1d first, struct cell1d last){
    return (first.primary.first < last.primary.first)? 1:0;
  }
  std::vector<struct cell1d> array;
  rsa1d(int size){
    array.reserve(size);
  }

  int insertCell(pair<int,int> prim, int val){
    struct cell1d c;
    c.primary = prim;
    c.value = val;
    array.push_back(c);
    return 0;
  }
  int sortArray(){
    sort(array.begin(), array.end(), cellcmp1d);

    for(int i = 0; i < array.size(); i++){
      auto a = array[i];
      hash.insert(make_pair(a.primary.first,i));
    }
    return 0;
  }

  int getArrayIndex(pair<int,int> i){
    int h = hash.at(i.first);
    return h;
  }
  /** Finds the integer corresponding to the hash-map, or the closest corresponding integer, rounding down in the case that the desired value does not exist in the hash map.
      Worst case linear in size of hash array, constant on average case.
  */
  int hash_back(int i, bool fwd){
    if(i == array.size()-1){
      return i;
    }
    if(i == 0){
      return i;
    }
    int h;
    auto eq = hash.equal_range(i); // constant on average
    auto index = distance(hash.begin(),eq.first); //constant
    return index   ;
  }
  /**
     Finds the maximum value in the given range of keys.
     @return pair<struct cell1d, int>
  */
  std::pair<struct cell1d,int> rangeMax(int min1, int max1){

    int maxVal = INT_MIN;
    struct cell1d maxCell;
    maxCell.primary = make_pair(-101,-101);
    if(min1 == -1 && max1 == -1){
      return make_pair(maxCell, maxVal);
    }
    int count = 0;
    int m,m2;
    m = hash_back(min1,false);
    m2 = array.size()-1;

    if(verbose) cout << "rangemax bounds: "<< min1 << ", " << max1 << "...";
    if(verbose) cout << "hashed bounds: "<< m << ", " << m2 << "...";
    if(verbose) cout << "rangemax: array.size()=" << array.size() << endl;
    for(int i = m; i <= m2; i++){
      count++;
      if(array[i].primary.first > max1){
        break;
      }
      if(array[i].primary.first >= min1 && array[i].primary.first <= max1){
        int tempmax = array[i].value;
	
        if(tempmax > maxVal){
          maxVal = tempmax;
          maxCell = array[i];	  
        }
      }
    }
    if(verbose) cout <<"1D rangemax: " << maxCell.primary.second << " value: " << maxVal << " after: "<<count<< " queries." << endl << endl;;
    return std::make_pair(maxCell, maxVal);
  }

  void printKey(pair<int,int> primary){
    auto i = getArrayIndex(primary);
    std::cout << "Primary key: "   << array[i].primary.first
              << " Value: "        << array[i].value
              << "\n";
  }

  void update(pair<int,int> primary, int value){
    auto i = getArrayIndex(primary);
    array[i].primary = primary;
    array[i].value = value;
  }

  void upgrade(pair<int,int> primary, int value){
    auto i = getArrayIndex(primary);
    int oldValue = array[i].value;
    array[i].value = (oldValue > value)? oldValue : value;
  }
};
#endif
