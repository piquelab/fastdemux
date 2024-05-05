
#ifndef COMMON_H
#define COMMON_H


#include <htslib/sam.h>


#define VERSION "0.1.1"

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))


// Function prototypes
int encodeACGT(char c);

unsigned long int packKmer(char *s,int k);

void unPackKmer(unsigned long int kmer,int k,char *s);

int find_read_index_for_ref_pos(const bam1_t *aln, int ref_pos);


#endif // COMMON_H


