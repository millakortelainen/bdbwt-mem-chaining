#include <random>
#include <vector>
#include <set>
#include <string>
#include <algorithm>
#include <iostream>
#include <regex>
#include <fstream>
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

string mutateText(string text, int blockCounts, double targetSim){
  int tlen = text.length();
  random_device rd;
  mt19937 gen(rd());
  set<char> alphabet;
  vector<pair<int,int>> blocks;
  cout << "Text lenght = " << tlen << endl;
  cout << "Target similarity = " << targetSim << endl;
  int editsLeft = (tlen-(tlen*targetSim));
  cout << "Doing "<< editsLeft << " edits" << endl;
  for(auto t : text){
    alphabet.insert(t);
  }
  if(blockCounts < 1){
    blockCounts = 1;
  }
  auto blockLen = tlen / blockCounts;

  for(int i = 0; i < blockCounts; i++){
    blocks.push_back(make_pair(i*blockLen, (i+1)*(blockLen-1)));
  }

  shuffle(blocks.begin(), blocks.end(), gen);
  while(editsLeft > 0){
    for(auto b : blocks){
      auto l = b.first;
      auto r = b.second;
      int max = (editsLeft < (r-l))? editsLeft : (r-l);
      set<int> editIndices;
      uniform_int_distribution<> distrib(l, r);
      for(int i = 0; i < max/blocks.size()+1; i++){
	editIndices.insert(distrib(gen));
      }
      uniform_int_distribution<> distrib2(0, alphabet.size()-2);
      for(auto e : editIndices){
	char tc = text.at(e);
	vector<char> changes;
	for(auto a : alphabet){
	  if(a != tc){
	    changes.push_back(a);
	  }
	}

	auto ind = distrib2(gen);
	text[e] = changes.at(ind);
	editsLeft--;
      }
    }
  }
  return text;
}

int main(int argc, char *argv[]){
  auto fileinput = readInputFromFasta(argv[1]);
  auto mutateSimil = stod(argv[2]);
  auto blockCounts = stod(argv[3]);

  auto mutated = mutateText(fileinput.at(0), blockCounts, mutateSimil);
  writeStringToFasta(argv[4], mutated);
}
