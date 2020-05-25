#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>

class rsa2d{
public:
  struct cell{
    int primary;
    std::vector<int> secondary;
    std::vector<int> value;
  };
  
  std::vector<cell> array;
  rsa2d(int size){
    for(int i = 0; i < size; ++i){
      struct cell c;
      c.primary = i;
      for(int j = 0; j < size; ++j){
	c.secondary.push_back(j);
	c.value.push_back(-1);
    }
      	array.insert(array.end(),c);
    }
  }

  int insertCell(int prim, int sec, int val){
    struct cell c;
    c.primary = prim;
    c.secondary.push_back(sec);
    c.value.push_back(val);
    array.insert(array.end(),c);
    return 0;
  }

  /**
     Finds the maximum value in the given range of keys.
     @return pair<struct cell, int>
  */
  std::pair<struct cell, int> rangeMax(int min1, int max1, int min2, int max2){
    if(min1 > array.size() || min2 > array.size()){
      throw "Minimum Before array beginning, throwing an error out.";
     }
    if(max1 > array.size()-1){
      max1 = array.size()-1;
    }
    if(max2 > array.size()-1){
      max2 = array.size()-1;
    }
    //    std::cout << "bounds: "<< min1 << ", " << max1 << "::" << min2 << ", " << max2 << "\t array size: " << array.size() << "\n";

    
    int maxVal = -1;
    struct cell maxCell;
    
    for(int i = min1; i <= max1; i++){
      if(i > array.size()-1){
	break;
      }
      for(int j = min2; j <= max2; j++){
	int tempmax = array[i].value[j];
	//	std::cout <<"tempmax: " << tempmax << "\n";

	if(tempmax > maxVal){
	  maxVal = tempmax;
	  maxCell = array[i];
	}
      }
    }
    //    std::cout << "maxVal = " << maxVal << "\n";
    return std::make_pair(maxCell, maxVal);
  }

  void printKey(int primary, int secondary){
    std::cout << "Primary key: "   << array[primary].primary
	      << " Secondary key: "<< array[primary].secondary[secondary]
	      << " Value: "        << array[primary].value[secondary]
	      << "\n";
  }

  void update(int primary, int secondary, int value){
    array[primary].primary = primary;
    array[primary].secondary[secondary] = secondary;
    array[primary].value[secondary] = value;
  }

  void upgrade(int primary, int secondary, int value){
    int oldValue = array[primary].value[secondary];
    array[primary].value[secondary] = (oldValue > value)? oldValue : value;
  }
};
  
