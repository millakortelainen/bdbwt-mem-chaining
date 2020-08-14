#ifndef BDBWT_MEM_CHAIN_IO_HH
#define BDBWT_MEM_CHAIN_IO_HH
#include <fstream>
#include <vector>
#include <stdlib.h>
#include <iostream>
#include "driver.hh"
#include <regex>
using namespace std;

vector<string> readInputFromFasta(string filename, int specificIndex = -1){
  vector<string> texts;
  string line;
  string text;
  filename = regex_replace(filename, regex("^\\s+"), ""); //trim whitespace
  filename = regex_replace(filename, regex("\\s+$"), "");
  ifstream fa(filename);
  int index = -1;
  while(getline(fa,line)){
    if(line == "" || line.at(0) == ';' || line.at(0) == '\n'){
      cout << line << endl;
      continue;
    }if(line.at(0) == '>'){
      index++;
      if(specificIndex < 0 || (specificIndex == index && specificIndex >= 0)){
	cout << line << endl;
      }
      if(text.length() > 0 && specificIndex < 0){
	texts.push_back(text);	
      }
      else if(text.length() > 0 && specificIndex >= 0 && specificIndex == index-1){
	texts.push_back(text);
	break;
      }
      text = "";
    }else{
      text.append(line);
    }
  }
  if(index == -1){
    cout << "no line starting with '>' found." << endl;
  } 
  texts.push_back(text);
  //  cout << text << endl;
  
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
Configuration readConfiguration(string filename){
  Configuration conf;
  string line;
  ifstream fa(filename);
  getline(fa,line);
  string delimiter = ">";
  string parse = line.substr(line.find(delimiter)+1, line.length());
  string filename1 = parse;
  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  int index1 = strtol(parse.c_str(),NULL,10);
  cout << "Fasta1: " << filename1 << " [" << index1 << "]" << endl;
  conf.text1 = readInputFromFasta(filename1, index1).at(0);
  conf.text1 = conf.text1.substr(0,conf.text1.length()-1);
  cout <<">";
  cout <<conf.text1.substr(0,30);
  cout << "..." << endl;

  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  string filename2 = parse;
  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  int index2 = strtol(parse.c_str(),NULL,10);
  cout << "Fasta2: " << filename2 << " [" << index2 << "]" << endl;
  conf.text2 = readInputFromFasta(filename2, index2).at(0);
  conf.text2 = conf.text2.substr(0,conf.text2.length()-1);
  cout <<">";
  cout << conf.text2.substr(0,30);
  cout << "..." << endl;
  
  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  int mode = strtol(parse.c_str(),NULL,10);
  cout << "Computation mode: " << mode << endl;
  conf.mode = mode;

  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  int size = strtol(parse.c_str(),NULL,10);
  cout << "Minimum Depth / K-mer size: " << size << endl;
  conf.minimumDepth = size;

  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  int minimizerWinSize = strtol(parse.c_str(),NULL,10);
  cout << "Window Size: " << minimizerWinSize << endl;
  conf.minimizerWindowSize = minimizerWinSize;

  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  int minimizerMergerCount = strtol(parse.c_str(),NULL,10);
  cout << "Minimizer merger count: " << minimizerWinSize << endl;
  conf.miniMergerCount = minimizerMergerCount;
  return conf;
  
}
#endif
