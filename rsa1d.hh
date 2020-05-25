#include <stdio.h>
#include <iostream>
#include <vector>
#include <utility>

class rsa1d{
public:
  struct cell{
    int primary;
    int value;
  };
  std::vector<cell> array;
  
  rsa1d(int size){
    for(int i = 0; i < size; ++i){
      struct cell c;
      c.primary = i;
      c.value = -1;
      array.insert(array.end(),c);
    }
  }

  int insertCell(int prim, int val){
    struct cell c;
    c.primary = prim;
    c.value = val;
    array.insert(array.end(),c);
    return 0;
  }

  /**
     Finds the maximum value in the given range of keys.
     @return pair<struct cell, int>
  */
  std::pair<struct cell,int> rangeMax(int min1, int max1){
    if(min1 > array.size()){
      throw "Minimum Before array beginning, throwing an error out.";
    }
    int maxVal = -1;
    struct cell maxCell;
    
    for(int i = min1; i <= max1; i++){
      if(i > array.size()-1){
	break;
      }
      int tempmax = array[i].value;
      
      if(tempmax > maxVal){
	maxVal = tempmax;
	maxCell = array[i];
      }
    }
    return std::make_pair(maxCell, maxVal);
  }

  void printKey(int primary){
    std::cout << "Primary key: "   << array[primary].primary
	      << " Value: "        << array[primary].value
	      << "\n";
  }

  void update(int primary, int value){
    array[primary].primary = primary;
    array[primary].value = value;
  }

  void upgrade(int primary, int value){
    int oldValue = array[primary].value;
    array[primary].value = (oldValue > value)? oldValue : value;
  }
};
  
