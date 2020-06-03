#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>
#include <bits/stdc++.h>
using namespace std;
class rsa1d{
private:
bool verbose = false;
  
public:
  //<value, identifier>
  struct cell{
    pair<int,int> primary;
    int value;
  };
  std::vector<cell> array;
  
  rsa1d(int size){
    array.reserve(size);
  }

  int insertCell(pair<int,int> prim, int val){
    struct cell c;
    c.primary = prim;
    c.value = val;
    array.push_back(c);
    return 0;
  }
  int getArrayIndex(pair<int,int> i){
    for(int a = 0; a < array.size(); a++){
      if(array[a].primary.first == i.first && array[a].primary.second == i.second){
	return a;
      }
    }
    return -1;
  }

  /**
     Finds the maximum value in the given range of keys.
     @return pair<struct cell, int>
  */
  std::pair<struct cell,int> rangeMax(int min1, int max1){
    if(verbose) cout << "rangemax bounds: "<< min1 << ", " << max1 << endl;
    if(verbose) cout << "rangemax: array.size()=" << array.size() << endl;
    
    int maxVal = INT_MIN;
    struct cell maxCell;
    maxCell.primary = make_pair(-1,-1);
    
    for(int i = 0; i < array.size(); i++){
      if(array[i].primary.first >= min1 && array[i].primary.first <= max1){
	int tempmax = array[i].value;
	
	if(tempmax > maxVal){
	  maxVal = tempmax;
	  maxCell = array[i];	  
	}
      }
    }
    if(verbose) cout <<"1D rangemax: " << maxCell.primary.second << "value: " << maxVal << endl;
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
  
