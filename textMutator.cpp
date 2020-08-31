#include <random>
#include <vector>
#include <set>
#include <string>
#include "io.hh"
#include <algorithm>
using namespace std;
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
