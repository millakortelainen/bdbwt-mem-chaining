#include <fstream>
using namespace std;

vector<string> readInputFromFasta(string filename){
  vector<string> texts;
  string line;
  string text;
  ifstream fa(filename);
  int index = -1;
  while(getline(fa,line)){
    if(line.at(0) == ';' || line.at(0) == '\n'){
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
