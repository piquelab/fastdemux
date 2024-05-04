
#ifndef COMMON_H
#define COMMON_H


#define VERSION "0.1.1"

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))


// Function prototypes
int encodeACGT(char c);

unsigned long int packKmer(char *s,int k);

void unPackKmer(unsigned long int kmer,int k,char *s);


#endif // COMMON_H


