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
  unordered_map<int,int> hash;
public:
struct cell2d{
  pair<int,int> primary;
  vector<pair<int,int>> secondary;
};
  static int cellcmp2d(struct cell2d first, struct cell2d last){
  return (first.primary.first < last.primary.first)? 1:0;
}
  static bool secondValueSort(const pair<int,int> &first, const pair<int,int> &second){
  return (first.first < second.first);
}
  
  std::vector<struct cell2d> array;
  //key, location in vector 'array'
  //unordered_map<int, int> hash;
  rsa2d(int size){
    array.reserve(size);
    // vector<cell> temp;
    // for(int i = 0; i < size; ++i){
    //   struct cell2d c;
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
    struct cell2d c;
    if(verbose) cout << "inserting prim: <" << prim.first << "," << prim.second << ">, sec: " << sec << ", val: " << val << "..."; 
    c.primary = prim;
    c.secondary.push_back(make_pair(sec,val));
    //    c.value.push_back(val);
    if(verbose) cout << "inserting, Secondary lenght ="<< c.secondary.size() << "...";
    array.push_back(c);
    if(verbose) cout << "success" << endl;
    //int n = array.size()-1;
    //hash.insert(make_pair(prim.first, n));
    return 0;
  }
  int sortArray(){
    sort(array.begin(), array.end(), cellcmp2d);

    int min = array[0].primary.first;
    int max = min;
    hash.insert(make_pair(array[0].primary.first, 0));
    for(int i = 1; i < array.size(); i++){
      auto a = array[i];
      if(a.primary.first > min){
	hash.insert(make_pair(a.primary.first, i));
	min = a.primary.first;
      }
    }
    return 0;
  }
  pair<int,int> getArrayIndex(pair<int,int> i, int s){
    int h = hash.at(i.first);
    for(int a = h; a < array.size(); a++){
      if(array[a].primary.first == i.first && array[a].primary.second == i.second){
	for(int b = 0; b < array[a].secondary.size(); b++){
	  if(array[a].secondary[b].first == s){
	    return make_pair(a,b);
	  }
	}
	array.at(a).secondary.push_back(make_pair(s, INT_MIN));
	return getArrayIndex(i,s);
      }
    }
    return make_pair(-1,-2);
  }
  /**
     Finds the maximum value in the given range of keys.
     @return pair<struct cell2d, int>
  */
  std::pair<struct cell2d, int> rangeMax(int min1, int max1, int min2, int max2){
    int m, m2;
    int count = 0;
    if(min1 == INT_MIN){
      m = 0;
    }else{
      m = hash.at(min1-1); //since input is r.left - l.left +1
    }
    if(max1 == INT_MAX){
      max1 = array.size()-1;
    }else{
      m2 = hash.at(max1);
    }
    
    int maxVal = INT_MIN;
    struct cell2d maxCell;
    maxCell.primary = make_pair(-1,-1);
    
    if(verbose) cout << "rangemax2D bounds: "<< min1 << ", " << max1 << "::" << min2 << ", " << max2 << "...";
    if(verbose) cout << "hashed   bounds: "<< m << ", " << m2 << "...";

    for(int i = m; i <= m2; i++){
      if(array[i].primary.first >= min1 && array[i].primary.first <= max1){
	for(int k = 0; k < array[i].secondary.size(); k++){
	  count++;
	  int sec = array[i].secondary.at(k).first;
	  if(array[i].secondary.size() > 1){
	    cout << "secondary: "<< sec;
	  }
	  if(sec >= min2 && sec <= max2){
	    int tempMax = array[i].secondary[k].second;
	    
	    if(tempMax > maxVal){
	      maxVal = tempMax;
	      maxCell = array[i];
	    }
	  }
	}
      }
    }
    if(verbose) cout <<"2D rangemax: " << maxCell.primary.second << " value: " << maxVal << "after: " << count << " queries." << endl;
    return std::make_pair(maxCell, maxVal);
  }
  
  void printKey(pair<int,int> primary, int secondary){
    auto i = getArrayIndex(primary, secondary);
    std::cout << "Primary key: "   << array[i.first].primary.first
	      << " Secondary key: "<< array[i.first].secondary[i.second].first
	      << " Value: "        << array[i.first].secondary[i.second].second
	      << "\n";
  }

  void update(pair<int,int> primary, int secondary, int value){
    auto i = getArrayIndex(primary, secondary);
    array[i.first].secondary[i.second].first = secondary;
    array[i.first].secondary[i.second].second = value;
    if(verbose) cout << "updated primary(" << i.first <<", "<< i.second <<"), secondary: " << secondary << " with value: " << value << endl;
  }

  void upgrade(pair<int,int> primary, int secondary, int value){
    auto i = getArrayIndex(primary, secondary);
    int oldValue = array[i.first].secondary.at(i.second).second;
    array[i.first].secondary.at(i.second).second = (oldValue > value)? oldValue : value;
  }
};
  
