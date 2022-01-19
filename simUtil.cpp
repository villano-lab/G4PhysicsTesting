#include <memory>
#include <fstream>
#include <iostream>
#include <sstream>

#include "TTree.h"
#include "TFile.h"

int *getList(int n,std::string filename)
{
  //make a list of ints, default contents zero
  int *out = (int*) std::malloc(n*sizeof(int));

  for(int i=0;i<n;i++){
    out[i]=0;
  }

  //read a file, and for each entry modify the contents
  //remember event numbers start at 1!
  std::ifstream in(filename.c_str(),std::ios::in);

  if(!in){
    std::cout << "File not present" << std::endl;
    return out;
  }

  int ev;
  in>>ev;
  while(!in.eof()){
    out[ev-1]=1;
    in>>ev;
  }

  return out;
}
void getListNew(int n,std::string filename,int *e,int *d,int *dI)
{
  //make a list of ints, default contents zero
  int *out = (int*) std::malloc(n*sizeof(int));

  for(int i=0;i<n;i++){
    out[i]=0;
  }

  //read a file, and for each entry modify the contents
  //remember event numbers start at 1!
  std::ifstream in(filename.c_str(),std::ios::in);

  if(!in){
    std::cout << "File not present" << std::endl;
    return;//out; // wasn't sure what the correct fix was here.
  }

  //allocate 
  //cout << "Allocating: " << n << endl;
  //int *e = (int*) malloc(n*sizeof(int));
  //int *d = (int*) malloc(n*sizeof(int));
  //int *dI = (int*) malloc(n*sizeof(int));

  int ev,det,didInteract;
  //in>>ev;
  std::string line;
  while(!getline(in,line).eof()){
    //cout << line << endl;
    std::istringstream iss(line);
    iss >> ev >> det >> didInteract;
    //cout << ev << endl;
    e[ev-1] = ev;
    d[ev-1] = det;
    dI[ev-1] = didInteract;
    //out[ev-1]=1;
    //in>>ev;
  }

  return;
}
bool insertEscapeInfo(std::string rootfile,std::string filename,bool force=false)
{
  //read the root file, find a Tree called metaData, overwrite it with new meta data
  //open the file and get the old tree
  TFile *f = new TFile(rootfile.c_str(),"UPDATE");

  bool havetree = (bool)f->GetListOfKeys()->FindObject("escapeInfo");

  if(havetree && !force) return true;

  TTree *oldt = (TTree*)f->Get("cascade");
  if(!oldt) return false;
  //make a tree
  TTree *t = new TTree("escapeInfo","escapeInfo");

  int ev,det,didInteract;
  t->Branch("ev",&ev,"ev/I");
  t->Branch("det",&det,"det/I");
  t->Branch("didInteract",&didInteract,"didInteract/I");

  //get escape info data 
  //extract the shift from the filename
  int *e,*d,*dI; 
  //MUST allocate
  e = (int*) malloc(oldt->GetEntries()*sizeof(int));
  d = (int*) malloc(oldt->GetEntries()*sizeof(int));
  dI = (int*) malloc(oldt->GetEntries()*sizeof(int));
  getListNew(oldt->GetEntries(),filename,e,d,dI);
  for(int i=0;i<oldt->GetEntries();i++){
    ev = e[i];
    det = d[i];
    didInteract = dI[i];
    t->Fill();
  }
  
  //write it back to the file
  t->Write("",TObject::kOverwrite);
  f->Close();

  return true;
}
