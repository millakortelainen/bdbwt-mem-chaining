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

  string delimiter = ">";

  getline(fa,line);
  string parse = line.substr(line.find(delimiter)+1, line.length());
  int verbosity = strtol(parse.c_str(),NULL,10);
  conf.verbosity = verbosity;

  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  string filename1 = parse;
  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  int index1 = strtol(parse.c_str(),NULL,10);
  if(verbosity > 2) cout << "Fasta1: " << filename1 << " [" << index1 << "]" << endl;
  conf.text1 = readInputFromFasta(filename1, index1).at(0);
  conf.text1 = conf.text1.substr(0,conf.text1.length()-1);
  if(verbosity > 2)cout <<">";
  if(verbosity > 2) cout <<conf.text1.substr(0,30);
  if(verbosity > 2)cout << "..." << endl;

  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  string filename2 = parse;
  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  int index2 = strtol(parse.c_str(),NULL,10);
  if(verbosity > 2)cout << "Fasta2: " << filename2 << " [" << index2 << "]" << endl;
  conf.text2 = readInputFromFasta(filename2, index2).at(0);
  conf.text2 = conf.text2.substr(0,conf.text2.length()-1);
  if(verbosity > 2)cout <<">";
  if(verbosity > 2)cout << conf.text2.substr(0,30);
  if(verbosity > 2)cout << "..." << endl;
  
  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  int mode = strtol(parse.c_str(),NULL,10);
  if(verbosity > 2) cout << "Computation mode: " << mode << endl;
  conf.mode = mode;

  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  int size = strtol(parse.c_str(),NULL,10);
  if(verbosity > 2)cout << "Minimum Depth / K-mer size: " << size << endl;
  conf.minimumDepth = size;

  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  int minimizerWinSize = strtol(parse.c_str(),NULL,10);
  if(verbosity > 2) cout << "Window Size: " << minimizerWinSize << endl;
  conf.minimizerWindowSize = minimizerWinSize;

  getline(fa,line);
  parse = line.substr(line.find(delimiter)+1, line.length());
  int minimizerMergerCount = strtol(parse.c_str(),NULL,10);
  if(verbosity > 2)cout << "Minimizer merger count: " << minimizerMergerCount << endl;
  conf.miniMergerCount = minimizerMergerCount;

  getline(fa,line);
  bool threadedBWT = true;
  parse = line.substr(line.find(delimiter)+1, line.length());
  int threadedBWTint = strtol(parse.c_str(),NULL,10);
  if(threadedBWTint == 0){
    threadedBWT = false;
  }
  if(verbosity > 2) cout << "BWT Threading: " << threadedBWT << endl;
  conf.threadedBWT = threadedBWT;

  getline(fa,line);
  bool printAbsentAndChains = false;
  parse = line.substr(line.find(delimiter)+1, line.length());
  int printAbsentAndChainsint = strtol(parse.c_str(),NULL,10);
  if(printAbsentAndChainsint == 1){
    printAbsentAndChains = true;
  }
  if(verbosity > 2)cout << "Verbose Chain and Absent sections: " << printAbsentAndChains << endl;
  conf.printAbsentAndChains = printAbsentAndChains;

  getline(fa,line);
  bool VerboseEditDistances = false; 
  parse = line.substr(line.find(delimiter)+1, line.length());
  int VerboseEditDistancesint = strtol(parse.c_str(),NULL,10);
  if(VerboseEditDistancesint == 1){
    VerboseEditDistances = true;
  }
  if(verbosity > 2)cout << "Verbose Chain and Absent sections: " << VerboseEditDistances << endl;
  conf.verboseEditDistances = VerboseEditDistances;

  getline(fa,line);
  bool rawChains = false;
  parse = line.substr(line.find(delimiter)+1, line.length());
  int rawChainsint = strtol(parse.c_str(),NULL,10);
  if(rawChainsint == 1){
    rawChains = true;
  }
  if(verbosity > 2) cout << "Printing raw chains: " << rawChains << endl;
  conf.rawChains = rawChains;

  getline(fa,line);
  bool chainStringSegments = false;
  parse = line.substr(line.find(delimiter)+1, line.length());
  int chainStringSegmentsint = strtol(parse.c_str(),NULL,10);
  if(chainStringSegmentsint == 1){
    chainStringSegments = true;
  }
  if(verbosity > 2) cout << "Printing chains as strings: " << chainStringSegments << endl;
  conf.chainStringSegments = chainStringSegments;
  return conf;

}
#endif
