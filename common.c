#include "common.h"


//encodes A,C,G,T in binary 00, 01, 10, 11
int encodeACGT(char c){
  c=c&0xDF; // 11011111 sets letter to uppercase
  switch(c){
  case 'A':
    return 0;
  case 'C':
    return 1;
  case 'G':
    return 2;
  case 'T':
    return 3;
  default:
    return -1;
  }
}

unsigned long int packKmer(char *s,int k){
  int Base,j;
  unsigned long int kmer=0;
  //  assert((k<<1)<sizeof(unigned long int));
  for(j=0;j<k;j++){
    Base=encodeACGT(s[j]);
    if(Base==-1)
      continue; // This will ignore if not ACGT and continue. Check if that is the desired behaviour. 
      //return 0xFFFFFFFFFFFFFFFF; //
    kmer=(kmer<<2)+Base;
  }
  return kmer;
}

/*
// It don't think this is faster than the option above. but I wanted to try. 
unsigned long int packKmer(char *s, int k) {
  static const int encoding[128] = {
    ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3,
    ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3
  };
  unsigned long int kmer = 0;
  for (int j = 0; j < k; ++j) {
    char c = s[j] & 0xDF; // Convert to uppercase
    if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
      //return 0xFFFFFFFFFFFFFFFF;  // I don't hink I can't comment this thing out. 
    }
    kmer = (kmer << 2) | encoding[c]; // Use bitwise OR for slight optimization
  }
  return kmer;
}
*/

void unPackKmer(unsigned long int kmer,int k,char *s){
  int j;
  char ACGT[]={"ACGT"};
  for(j=k-1;j>=0;j--){
    s[j]=ACGT[(kmer & 0x03)];
    kmer=(kmer>>2);
  }
  s[k]='\0'; //string end. 
}

