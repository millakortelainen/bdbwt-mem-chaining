#ifndef BDBWT_MEM_CHAIN_IO_HH
#define BDBWT_MEM_CHAIN_IO_HH
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <iostream>
using namespace std;

vector<string> readInputFromFasta(string filename){
  vector<string> texts;
  string line;
  string text;
  ifstream fa(filename);
  int index = -1;
  while(getline(fa,line)){
    if(line == "" || line.at(0) == ';' || line.at(0) == '\n'){
      cout << line << endl;
      continue;
    }if(line.at(0) == '>'){
      cout << line << endl;
      index++;
      if(text.length() > 0){
	texts.push_back(text);
	text = "";
      }
      continue;
    }
    if(index == -1){
      cout << "no line starting with '>' found." << endl;
      break;
    }
    //    cout << line << endl;
    text.append(line);
  }
  texts.push_back(text);
  return texts;
}
void writeStringToFasta(string filename, string fasta){
  ofstream f (filename, ofstream::out);
  f << ">unnamed";
  for(int i = 0; i < fasta.length(); i++){
    if((i % 70) == 0){
      f << endl;
    }
    f << fasta[i];
  } 
}
#endif
