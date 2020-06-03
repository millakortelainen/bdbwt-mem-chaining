#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>
#include <unordered_map>
#include <bits/stdc++.h>
using namespace std;
class rsa2d{
private:
bool verbose = false;

public:
  struct cell{
    pair<int,int> primary;
    vector<int> secondary;
    vector<int> value;
  };
  
  std::vector<cell> array;
  //key, location in vector 'array'
  //unordered_map<int, int> hash;
  rsa2d(int size){
    array.reserve(size);
    // vector<cell> temp;
    // for(int i = 0; i < size; ++i){
    //   struct cell c;
    //   c.primary = make_pair(i,-1);
    //   for(int j = 0; j < size; ++j){
    // 	c.secondary.push_back(j);
    // 	c.value.push_back(INT_MIN);
    // 	temp.push_back(c);
    //   }
    //   array.insert(array.end(),temp);
    // }
  }

  int insertCell(pair<int,int> prim, int sec, int val){
    struct cell c;
    if(verbose) cout << "inserting prim: <" << prim.first << "," << prim.second << ">, sec: " << sec << ", val: " << val << "..."; 
    c.primary = prim;
    c.secondary.push_back(sec);
    c.value.push_back(val);
    if(verbose) cout << "inserting, Secondary lenght ="<< c.secondary.size() << "...";
    array.push_back(c);
    if(verbose) cout << "success" << endl;
    return 0;
  }
  pair<int,int> getArrayIndex(pair<int,int> i, int s){
    for(int a = 0; a < array.size(); a++){
      if(array[a].primary.first == i.first && array[a].primary.second == i.second){
	for(int b = 0; b < array[a].secondary.size(); b++){
	  if(b == s){
	    return make_pair(a,b);
	  }
	}
	array.at(a).secondary.push_back(s);
	array.at(a).value.push_back(INT_MIN);
	return getArrayIndex(i,s);
      }
    }
    return make_pair(-1,-2);
  }
  /**
     Finds the maximum value in the given range of keys.
     @return pair<struct cell, int>
  */
  std::pair<struct cell, int> rangeMax(int min1, int max1, int min2, int max2){
    int maxVal = INT_MIN;
    struct cell maxCell;
    maxCell.primary = make_pair(-1,-1);
    
    if(verbose) cout << "rangemax bounds: "<< min1 << ", " << max1 << "::" << min2 << ", " << max2 << endl;
    if(verbose) cout << "rangemax: array.size()=" << array.size() << endl;
    for(int i = 0; i < array.size(); i++){
      if(array[i].primary.first >= min1 && array[i].primary.first <= max1){
	for(int k = 0; k < array[i].secondary.size(); k++){
	  int sec = array[i].secondary[k];
	  if(sec >= min2 && sec <= max2){
	    int tempMax = array[i].value[k];
	    	  
	    if(tempMax > maxVal){
	      maxVal = tempMax;
	      maxCell = array[i];
	    }
	  }
	}
      }
    }
    if(verbose) cout <<"2D rangemax: " << maxCell.primary.second << " value: " << maxVal << endl;
    return std::make_pair(maxCell, maxVal);
  }
  
  void printKey(pair<int,int> primary, int secondary){
    auto i = getArrayIndex(primary, secondary);
    std::cout << "Primary key: "   << array[i.first].primary.first
	      << " Secondary key: "<< array[i.first].secondary[i.second]
	      << " Value: "        << array[i.first].value[i.second]
	      << "\n";
  }

  void update(pair<int,int> primary, int secondary, int value){
    auto i = getArrayIndex(primary, secondary);
    array[i.first].secondary[i.second] = secondary;
    array[i.first].value[i.second] = value;
    if(verbose) cout << "updated primary(" << i.first <<", "<< i.second <<"), secondary: " << secondary << " with value: " << value << endl;
  }

  void upgrade(pair<int,int> primary, int secondary, int value){
    auto i = getArrayIndex(primary, secondary);
    int oldValue = array[i.first].value.at(i.second);
    array[i.first].value.at(i.second) = (oldValue > value)? oldValue : value;
  }
};
  
