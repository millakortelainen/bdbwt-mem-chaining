#include "mem.hh"
using namespace std;
void naiveOutput(BD_BWT_index<> index, BD_BWT_index<> index2, vector<tuple<int,int,int>> memVector, string text1 = "", string text2 = "",bool verbose = false){
  auto retSA = buildSAfromBWT(index); //RetSA builds SA array for the given text from it's BWT transform without having to use the extra space from permutating whole original text.
  auto retSA2 = buildSAfromBWT(index2);

  //Naive returning of the intervals
  for(auto a : memVector){
    int i, j, depth;
    tie(i,j,depth) = a;
    int begin_i= retSA[i]+1;
    int begin_j= retSA2[j]+1;
    
    //Handling the specific special case when SA^S[i] = text.size(); We use index.size()-1 instead so wouldn't need to keep the text in memory at all.
    if(begin_i >= index.size()-1){
      begin_i = index.size()-retSA[i]-1; //Index.size() will always be text.size()+1 due to the added END marker. Furthermore, we need to minus one to get proper offset from general case of begin_i = retSA[i]+1;
    }
    //Handling the specific special case when SA^T[j] = text2.size(); We use index2.size()-1 instead so wouldn't need to keep the text in memory at all.
    if(begin_j >= index2.size()-1){
      begin_j = index2.size()-retSA2[j]-1; //Analogously to above. 
    }
    int end_i  = begin_i+depth-1;
    int end_j  = begin_j+depth-1;
    
    cout << "Given tuple (i,j,d): "<< "(" << i << ","<<j <<"," << depth << ")" << "\n";
    
    cout << "S: ["<< begin_i <<","<< end_i <<"]" << "-->\t";
    if(text1.size() > 1 && verbose){
      for(int b = begin_i; b <= end_i; b++){
	cout << text1[b]; //For verificiation of results
      }
      cout << "\n";
    }
    if(text2.size() > 1 && verbose){
      cout << "T: ["<< begin_j <<","<< end_j <<"]" << "-->\t";
      for(int b = begin_j; b <= end_j; b++){
	cout << text2[b]; //For verification of results
      }
      cout << "\n";
    }
  }
}
vector<tuple<int,int,int>> batchOutput(BD_BWT_index<> index, BD_BWT_index<> index2, vector<tuple<int,int,int>> memVector, bool verbose = false){
  std::vector<struct occStruct> Ipairs;
  std::vector<struct occStruct> Ipairs2;
  std::vector<bool> marked1(index.size(),false);
  std::vector<bool> marked2(index2.size(),false);
  vector<tuple<int,int,int>> retVector;
  int p = 0;
  
  for(auto m : memVector){
    struct occStruct newOcc;
    struct occStruct newOcc2;
    int i,j,d;
    tie(i,j,d) = m;
    
    newOcc.key  = i;
    newOcc2.key = j;
    newOcc.primary  = p;
    newOcc2.primary = p;
    p++;
    
    marked1[i] = true;
    marked2[j] = true;
    Ipairs.push_back(newOcc);
    Ipairs2.push_back(newOcc2);
  }
  
  auto bl1 = batchLocate(Ipairs,marked1,index);
  auto bl2 = batchLocate(Ipairs2,marked2,index2);
  int maxkey = -1;
  for(int k = 0; k < memVector.size(); k++){
    int i,j,d;
    tie(i,j,d) = memVector[bl1[k].primary];
    int ik = bl1[k].key+1;
    int jk = bl2[k].key+1;
    if(ik+1 >= index.size()-1){ //Special case when interval begins at the beginning, SA[i]=text.size()
      ik = index.size()-ik;
    }
    if(jk+1 >= index2.size()-1){ //Special case when interval begins at the beginning, SA[i]=text.size()
      jk = index2.size()-jk;
    }
    retVector.push_back(make_tuple(ik, jk, d));
  }

  if(verbose){
    for(auto k : retVector){
      int i,j,d;
      tie(i,j,d) = k;
      cout << "triple: (" << i <<","<< j <<","<< d <<")\n";
    }
  }
  return retVector;
}
  
int main(){
  string text;
  string text2;
  switch(4){
  case 1: {
    text  = "GTGCGTGATCATCATTT";
    text2 = "AGTGCAAAGTGATTACC";
    break;
  }
  case 2:{
    text  = "ASDKISSAIKALAS";
    text2 = "ASDKASSAIAKALA";
    break;
  }
  case 3: {
    text  = "GTGCGTGTTCATCATTT";
    text2 = "GTGCGTGATCATCATTT";
    break;
  }
  case 4: {
    text  = "GTGCGTGATCATCATTTA";
    text2 = "AGTGCGCGTGACATCTTT";
    break;
  }
  case 500: {
    text = "gttcaccatttaaataatcttcaatatcaacacgcgaagctcgcttgcagggatgaactgaatagacctgtttactccggaaaagcaagactatcctggtgctgatgctacggtacattgttcttggcacgattacggactattcacactgaatccgggtggggagggccttatggacacgtaatatgcgcgtactggttggcgttgtagacgcgcaacttcatcgataatctgactgcctgacaagctaccagcaatacgttactccatcccgctatcctcggtactgcttgcggtgtcaccccgttaagtgacgtcctgttcgcggctaggctacgagttgcgttaatgcactctgaatcagaattccgcagcgttaagctggcttcaccagcgtcttcggtctgacttaaacctactcccgacatttctacagtgactactgtgtacgccccacgaagtcaaccccgagctacacctaaccggcctccagcactgcc";
    text2 = "aattcgaatagttagctgacgtacgacatgttaccttaataatataactggtgtccgcgactgagtgctctcctacctcccacgagcctcaggaaaaacgtctttaaatctctacccggagctgtttaaggggaagccaactcgaacctagcagggcattaaatttgtattgcaccaaaacgaccggcttaacattccgtgtctcactggacggaaaaccaacctaagcagtatttggcctcctggtaggcgaaccatctacggtggaccgtataatcggactaaccggcaggtttacacttcgcaatgctacgctgcccagggccgggcccccagtaggtttgcactgtagagggagggccggagtgtatcccccatcggtaactctacatatgcgcaagccgccctgggcaagatcccatcccactcgtgtggctctcgcgccgggtggattgtacgatcggaatcctctggggacgcgcgttcagtaacttcgctta";
    break;
  }
  }
  minimumDepth = 3;  
  BD_BWT_index<> index((uint8_t*)text.c_str());
  BD_BWT_index<> index2((uint8_t*)text2.c_str());
  
  auto mems  = bwt_mem2(index, index2);//Find MEMS between two BDBWT indexes.
  sort(mems.begin(), mems.end(), memSort); //Proper sorting of the tuples with priority order of i --> d --> j
  auto filtered = filterMems(mems);
  
  pretty_print_all(index,text);
  pretty_print_all(index2,text2);

  naiveOutput(index,index2,filtered,text,text2, true);
  cout << endl;
  auto bo = batchOutput(index,index2,filtered, true);

  cout << endl;
  sort(bo.begin(), bo.end(), memSort);
  vector<Interval_pair> Ipairs;
  for(auto b : bo){
    int i,j,d;
    tie(i,j,d) = b;
    if(i <= j){
      Interval_pair temp(i,i+d-1, j,j+d-1);
      cout << temp.toString() << endl;
      Ipairs.push_back(temp);
    }
  }
  sort(Ipairs.begin(), Ipairs.end(), intervalSort);
  auto chains = chaining(Ipairs, text2.size());
  for(int i = 0; i < chains.size(); i++){
    cout <<"Chain["<< i << "]: " << chains[i] << "\n";
  }
  vector<pair<Interval_pair, int>> chainPairs;
  int lf, lr, lc;
  cout << Ipairs.size() << endl;
  chainPairs.push_back(make_pair(Ipairs[0],chains[0]));
  lf = Ipairs[0].forward.right;
  lr = Ipairs[0].reverse.right;
  lc = chains[0];
  for(int i = 1; i < Ipairs.size(); i++){
    if((Ipairs[i].forward.left < lf || Ipairs[i].reverse.left < lr)){
      if(lc < chains[i]){
	cout << Ipairs[i].forward.left <<" < " << lf <<" && "<< lc<< " < "<< chains[i] << endl;;
	chainPairs.pop_back();
	chainPairs.push_back(make_pair(Ipairs[i],chains[i]));
	lf = Ipairs[i].forward.right;
	lr = Ipairs[i].reverse.right;
      }
      continue;      
    }
    // if(Ipairs[i].reverse.left < lr && lc < chains[i]){
    //   cout << Ipairs[i].reverse.left <<" < " << lr <<" && "<< lc<< " < "<< chains[i] << endl;;
    //   chainPairs.pop_back();
    //   chainPairs.push_back(make_pair(Ipairs[i],chains[i]));
    //   lr = Ipairs[i].reverse.right;
    //   lc = chains[i];
    //   continue;
    //    }
    chainPairs.push_back(make_pair(Ipairs[i],chains[i]));
    continue;
  }
  int totalScore = 0;
  for(auto c : chainPairs){
    auto i = c.first;
    auto cp = c.second;
    totalScore += cp;
    cout << i.toString() << " score: " << cp << ", total score = " << totalScore << endl;
  }
}
  

