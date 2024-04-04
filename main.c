#include <stdio.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include "khash.h"


#include <omp.h>

// Define the hash table for storing barcodes (CB tags) with an index to array
KHASH_MAP_INIT_STR(barcode,int)

#define SQRT2 1.41421356237
#define BAM_BUFF_SIZE 3600
#define OMP_THREADS 6


// Need to make a new hash table with values that contain info to store barcode or index to bc, 
// bci index of the barcode. 
// rc  reference allele count.
// ac  alternate allele count. 
// Define the structure to hold the barcode index, reference allele count, and alternate allele count
typedef struct {
  int bc_index;
  int ref_count;
  int alt_count;
} bc_info_t;

// Define the hash map type
KHASH_MAP_INIT_INT(bc_info_hash_t, bc_info_t)


// Define the hash table for storing barcodes (CB tags) with an index to array, 
KHASH_MAP_INIT_INT(bc_hash_t,int) //converting bc to integer, 32 bits should be enough but usign 64 for longer BC.

/*
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
*/


unsigned long int packKmer(char *s, int k) {
  static const int encoding[128] = {
    ['A'] = 0, ['C'] = 1, ['G'] = 2, ['T'] = 3,
    ['a'] = 0, ['c'] = 1, ['g'] = 2, ['t'] = 3
  };
  unsigned long int kmer = 0;
  for (int j = 0; j < k; ++j) {
    char c = s[j] & 0xDF; // Convert to uppercase
    if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
      //return 0xFFFFFFFFFFFFFFFF;
    }
    kmer = (kmer << 2) | encoding[c]; // Use bitwise OR for slight optimization
  }
  return kmer;
}


void unPackKmer(unsigned long int kmer,int k,char *s){
  int j;
  char ACGT[]={"ACGT"};
  for(j=k-1;j>=0;j--){
    s[j]=ACGT[(kmer & 0x03)];
    kmer=(kmer>>2);
  }
  s[k]='\0';
}


// Function to load valid barcodes from a file into a hash table
khash_t(bc_hash_t) *load_barcodes(char *filename) {
    khash_t(bc_hash_t) *hbc = kh_init(bc_hash_t);
    FILE *fp = fopen(filename, "r");
    char line[256];
    int i=0;
    while (fgets(line, sizeof(line), fp)) {
        int ret;
        khint_t k;
        line[strcspn(line, "\n")] = 0; // Remove newline character
        k = kh_put(bc_hash_t, hbc, packKmer(line,16), &ret); // Need to add logic for different size BC
	assert(ret>0); // Making sure the key is not there. No repeats allowed.
	kh_val(hbc, k)=i++;  //This is the bc_index for later. 
    }
    fclose(fp);
    return hbc;
}


// Function to find the read sequence index corresponding to a reference position
// based on the CIGAR string of the alignment. Returns -1 if the position is not covered.
/* 
int find_read_index_for_ref_pos(const bam1_t *aln, int ref_pos) {
  uint32_t *cigar = bam_get_cigar(aln);
  int read_pos = 0; // Position in the read
  int ref_pos_aligned = aln->core.pos; // Aligned position in the reference

  for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
    int op = cigar[i] & BAM_CIGAR_MASK;
    int len = cigar[i] >> BAM_CIGAR_SHIFT;

    switch (op) {
    case BAM_CMATCH:
    case BAM_CEQUAL:
    case BAM_CDIFF:
      if (ref_pos_aligned + len > ref_pos) {
	return read_pos + (ref_pos - ref_pos_aligned);
      }
      read_pos += len;
      ref_pos_aligned += len;
      break;
    case BAM_CINS:
    case BAM_CSOFT_CLIP:
      read_pos += len;
      break;
    case BAM_CDEL:
    case BAM_CREF_SKIP:
      if (ref_pos_aligned + len > ref_pos) {
	return -1; // Reference position falls within a deletion
      }
      ref_pos_aligned += len;
      break;
    default:
      // Other CIGAR operations are not considered here
      break;
    }
  }
  return -1; // Reference position is not covered by the read
}

*/

int find_read_index_for_ref_pos(const bam1_t *aln, int ref_pos) {
  uint32_t *cigar = bam_get_cigar(aln);
  int read_pos = 0; // Position in the read
  int ref_pos_aligned = aln->core.pos; // Aligned position in the reference

  for (uint32_t i = 0; i < aln->core.n_cigar; ++i) {
    int op = cigar[i] & BAM_CIGAR_MASK;
    int len = cigar[i] >> BAM_CIGAR_SHIFT;

    // Optimization for deletion and reference skip
    if ((op == BAM_CDEL || op == BAM_CREF_SKIP) && ref_pos_aligned + len > ref_pos) {
      return -1; // Reference position falls within a deletion or skip
    }

    switch (op) {
    case BAM_CMATCH:
    case BAM_CEQUAL:
    case BAM_CDIFF:
      if (ref_pos_aligned + len > ref_pos) {
	return read_pos + (ref_pos - ref_pos_aligned);
      }
      read_pos += len;
      ref_pos_aligned += len;
      break;
    case BAM_CINS:
    case BAM_CSOFT_CLIP:
      read_pos += len;
      break;
    case BAM_CDEL:
    case BAM_CREF_SKIP:
      ref_pos_aligned += len;
      break;
      // Other CIGAR operations are not considered here
    default:
      break;
    }
  }
  return -1; // Reference position is not covered by the read
}




// Main function to process BAM and VCF files
int main(int argc, char *argv[]) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <BAM file> <VCF file> <Barcode file> <Output file>\n", argv[0]);
        return 1;
    }

    char *bam_filename = argv[1];
    char *vcf_filename = argv[2];
    char *barcode_filename = argv[3];
    char *output_filename = argv[4];

    int nbcs=0;

    
    omp_set_num_threads(OMP_THREADS); //Adjust this to env variable.... or option


 
    // Load barcodes
    fprintf(stderr, "Reading the barcode file %s\n", barcode_filename);
    khash_t(bc_hash_t) *hbc = load_barcodes(barcode_filename);
    nbcs = kh_size(hbc);
    fprintf(stderr, "Number of barcodes is: %d\n", nbcs);

    // Open BAM file
    fprintf(stderr, "Opening BAM file %s\n", bam_filename);
    samFile *sam_fp = sam_open(bam_filename, "r");
    if (sam_fp == NULL) {
        fprintf(stderr, "Failed to open BAM file %s\n", bam_filename);
        return 1;
    }

    // multithreaded for sam, does not seem to make it faster, but slower
    //hts_set_threads(sam_fp, 6); // 6 threads?

    bam_hdr_t *bam_hdr = sam_hdr_read(sam_fp);
    hts_idx_t *bam_idx = sam_index_load(sam_fp, bam_filename); // Load BAM index
    if (bam_idx == NULL) {
        fprintf(stderr, "Failed to open BAM index for %s\n", bam_filename);
        sam_close(sam_fp);
        return 1;
    }


    // Open VCF file
    fprintf(stderr, "Open VCF file %s\n", vcf_filename);
    htsFile *vcf_fp = bcf_open(vcf_filename, "r");
    if (vcf_fp == NULL) {
        fprintf(stderr, "Failed to open VCF file %s\n", vcf_filename);
        return 1;
    }
    bcf_hdr_t *vcf_hdr = bcf_hdr_read(vcf_fp);
    bcf1_t *rec = bcf_init();

    if (!vcf_hdr) {
      fprintf(stderr, "Error: Could not read header from VCF file %s\n", vcf_filename);
      return 1;
    }

    // Process VCF and BAM files...
    // Assuming bam_hdr, sam_fp, vcf_fp, vcf_hdr are already defined and initialized
    int nsmpl = bcf_hdr_nsamples(vcf_hdr); // Number of samples

    // Iterate over sample names in the VCF header
    fprintf(stderr,"There are %d samples in the vcf file: \n",bcf_hdr_nsamples(vcf_hdr));
    for (int i = 0; i < nsmpl; i++) {
      const char *sample_name = vcf_hdr->samples[i]; // Get sample name
      fprintf(stderr,"%s\t",sample_name);
    }
    fprintf(stderr,"\n");


    int *bc_count_snps = (int *)malloc(nbcs * sizeof(int));  // Counts SNPs per barcode with >0 reads/UMIs.
    int *bc_count_umis = (int *)malloc(nbcs * sizeof(int)); // Counts UMIs in SNPs per barcode. 
    float *bc_mat_dlda = (float *)malloc(nsmpl * nbcs * sizeof(float)); //Matrix nbcs rows and nsmpl columns, to store DLDA results.  


    int snpcnt=0;

    for (int i = 0; i < nbcs; i++) bc_count_snps[i]=0;
    for (int i = 0; i < nbcs; i++) bc_count_umis[i]=0;
    for (int i = 0; i < nbcs*nsmpl; i++) bc_mat_dlda[i]=0.0;

    while (bcf_read(vcf_fp, vcf_hdr, rec) == 0) { // Iterate over VCF records
      if (bcf_unpack(rec, BCF_UN_STR) < 0) continue; // Unpack the current VCF record

      // Filter bi-allelic SNPs
      if (rec->n_allele != 2 || strlen(rec->d.allele[0]) != 1 || strlen(rec->d.allele[1]) != 1) continue;

      // Retrieve SNP position and alleles
      int pos = rec->pos + 1; // Convert to 1-based position
      char ref_allele = rec->d.allele[0][0];
      char alt_allele = rec->d.allele[1][0];

      // Initialize allele counts
      int ref_count = 0, alt_count = 0;

      //Extract dosages from vcf. 
      float *dosages=NULL;
      float alf=1E-5; //Instead of 0 I put this so variance does not become to small and avoid division by zero. 
      int nds_arr = 0;

      snpcnt++;

      if(snpcnt > 20000) break;  // for debugging just a small number of SNPs

      if (bcf_get_format_float(vcf_hdr, rec, "DS", &dosages, &nds_arr) < 0) continue;

      for (int i = 0; i < nsmpl; i++) {
	//	fprintf(stderr,"%f ",dosages[i]);
	alf +=dosages[i];
      }
      alf=alf/(nsmpl*2);

      if(snpcnt%1000==1)
	fprintf(stderr, "Processing %d: %d %s %d %c %c %f\n",snpcnt,rec->rid, bcf_hdr_id2name(vcf_hdr, rec->rid), rec->pos, ref_allele, alt_allele,alf);      

      if(alf * (1-alf) > 0.2) continue; 

      //using dosage array to store the DLDA weights. 
      for (int i = 0; i < nsmpl; i++)    
	dosages[i] = (dosages[i] - 2.0 * alf) / (SQRT2 * alf * (1-alf));

      /* GENOTYPE based DOSAGE */
      /*
      int *dosages = (int *)malloc(nsmpl * sizeof(int)); // Allocate memory for dosages
      // Temporary storage for genotype data
      int32_t *gt_arr = NULL, ngt_arr = 0;
      // Extract the GT field
      int ngt = bcf_get_genotypes(vcf_hdr, rec, &gt_arr, &ngt_arr);
      if (ngt <= 0){
	free(gt_arr);
	continue;
      }
      // Process each sample
      for (int i = 0; i < nsmpl; i++) {
        // Genotypes are stored in the format: allele1 | allele2, with -1 indicating missing data
        // Calculate the dosage: 0 for reference allele, 1 for each alternate allele
        int dosage = 0;
        for (int j = 0; j < ngt / nsmpl; j++) { // Loop over each allele (diploid: 2 alleles)
	  if (gt_arr[i * (ngt / nsmpl) + j] == bcf_int32_vector_end) break; // Reached the end of alleles for this sample
	  int allele_index = bcf_gt_allele(gt_arr[i * (ngt / nsmpl) + j]);
	  if (allele_index > 0) { // Count only alternate alleles
	    dosage += allele_index; // Assuming biallelic. For multiallelic sites, this needs adjustment
	  }
        }
        dosages[i] = dosage;
	fprintf(stderr,"%d ",dosage);
      }
      fprintf(stderr,"\n");
      // Cleanup
      free(gt_arr);
      */

      //Need to chr are the same 
      // Setup pileup
      hts_itr_t *iter = sam_itr_queryi(bam_idx, rec->rid, rec->pos, rec->pos + 1);
      if (iter != NULL){

      //           fprintf(stderr,"iter not null\n");
      // Assuming `barcodes` is a khash_t(barcode)* containing valid CB tags
      // and `reads_hash` is a khash_t(barcode)* used to track unique CB+UB combinations
      khash_t(bc_info_hash_t) *cb_h = kh_init(bc_info_hash_t); // Hash table for counting reads for CB with unique UMI
      khash_t(barcode) *cb_ub_h = kh_init(barcode); // Hash table for counting unique CB+UB combos as STRING think to convert to integer as well

      bam1_t *b[BAM_BUFF_SIZE];
      int j=0;
      for(j=0;j<BAM_BUFF_SIZE; j++)
	b[j] = bam_init1(); // Initialize a container for a BAM record.

      while(1){
      for(j=0;j<BAM_BUFF_SIZE; j++)
	if(sam_itr_next(sam_fp, iter, b[j])<0) break;
    
      if(j==0)
	break;
      //schedule(dynamic)
#pragma omp parallel for 
      for(int jj=0;jj<j;jj++){
	//while (sam_itr_next(sam_fp, iter, b) >= 0) {
	// Extract CB tag (cell barcode)
	uint8_t *cb_ptr = bam_aux_get(b[jj], "CB");
	if (cb_ptr==NULL) continue; // Skip if no CB tag
	char* cb = bam_aux2Z(cb_ptr); // Convert to string
	khint_t k;

	// Follow the CIGAR and get index to the basepair in the read
	int read_index = find_read_index_for_ref_pos(b[jj], rec->pos);
	if (read_index < 0) continue; // Skip if the position is not covered by the read

	unsigned long int kmer=packKmer(cb,16);
	// Check if the CB tag is in the list of valid barcodes
	khint_t kcb = kh_get(bc_hash_t, hbc, kmer);
	if (kcb == kh_end(hbc)) continue; // Skip if CB not found among valid barcodes 

  
	// Extract UB tag (unique molecular identifier)
	uint8_t *ub_ptr = bam_aux_get(b[jj], "UB");
	char ub[100]; // Assuming UB won't exceed this length
	if (ub_ptr) strcpy(ub, bam_aux2Z(ub_ptr)); // Convert to string if present
	else strcpy(ub, ""); // Use an empty string if no UB tag
    
	// Construct a unique key for CB+UB combination
	char cb_ub_key[200]; // Ensure this is large enough to hold CB+UB+delimiter
	sprintf(cb_ub_key, "%s-%s", cb, ub);

    
	// Check if we have already processed this CB
	int ret=0; 

#pragma omp critical
	{
	k = kh_put(bc_info_hash_t, cb_h, kmer , &ret); // Add it to the hash table
	
	assert(ret>=0);
	if(ret>0){ 
	  kh_val(cb_h,k).bc_index=kh_val(hbc,kcb);
	  kh_val(cb_h,k).ref_count=0;
	  kh_val(cb_h,k).alt_count=0;
	}
	}
	// only if unique UMI	  kh_val(cb_h,k)++;	

	// Check if we have already processed this CB+UB combination
#pragma omp critical
	{
	int k2 = kh_put(barcode, cb_ub_h, strdup(cb_ub_key), &ret);
	}
	if (ret>=0) { // This is a new, unique CB+UB combination
	  //  I could do majority voting for multiple cominations, but I will pick the first one for now. 
        
	  uint8_t *seq = bam_get_seq(b[jj]); // Get the sequence
	  // Get base at the read's position corresponding to the SNP
	  char base = seq_nt16_str[bam_seqi(seq, read_index)]; 

	  // Assuming ref_allele and alt_allele are char variables holding the reference and alternate alleles respectively
	  if (base == ref_allele) { 
	    ref_count++;
#pragma omp critical
	    kh_val(cb_h,k).ref_count++; 
	  }
	  else if (base == alt_allele) { 
	    alt_count++;
#pragma omp critical
	    kh_val(cb_h,k).alt_count++; 
	  }
	}// if ret
	    //	}//omp critical
       }//For omp
      }//while(1)

      for(int jj=0;jj<BAM_BUFF_SIZE;jj++)
	bam_destroy1(b[jj]); // Clean up BAM record container
      //}//pragma parallel 
     
      //I could use a parallel here too, perhaps. Barcodes are independent. 
#pragma omp parallel for schedule(dynamic)
      for (khint_t k = kh_begin(cb_h); k != kh_end(cb_h); ++k)  // traverse
	if (kh_exist(cb_h, k))            // test if a bucket contains data
          if ((kh_val(cb_h,k).ref_count > 0) || (kh_val(cb_h,k).alt_count > 0)) {
	    // If we want to write a pileup. 
	    // You may want to write this information to a file instead of printing it. This will allow for ASE calculations down the line. 
	    // Here we should do the demuxlet calculations using dosage and store results in a matrix.
	    //    printf("--%s\t%d\t%s\t%f\t%d\t%d\t%d\t%d\t%d\n", bcf_hdr_id2name(vcf_hdr, rec->rid), pos, rec->d.id, alf, ref_count, alt_count, 
	    // 	   kh_val(cb_h,k).bc_index, kh_val(cb_h,k).ref_count,kh_val(cb_h,k).alt_count);
    	    
	    assert(kh_val(cb_h,k).bc_index<nbcs);
	    
	    bc_count_snps[kh_val(cb_h,k).bc_index]++; // Increase number of SNPs seen for this barcode. 
	    bc_count_umis[kh_val(cb_h,k).bc_index]+= kh_val(cb_h,k).ref_count + kh_val(cb_h,k).alt_count; // Increase number of UMIs in SNPs seen for this barcode. 
	    
	    float val = kh_val(cb_h,k).alt_count - alf * bc_count_umis[kh_val(cb_h,k).bc_index];

	    //Loop across individuals for a barcode. 
	    for(int i = 0; i < nsmpl; ++i)
	      //The weight could be pre-calculated based on alf and dosage for each sample. dosage becomes the wight
	      // weigths: dosage[i] = ((dosages[i] - 2 * alf) / (SQRT2 * alf * (1-alf)))
	      bc_mat_dlda[ (kh_val(cb_h,k).bc_index * nsmpl) + i ] += dosages[i] * val;
	  }
      
      // Clear the reads hash table for the next SNP
      kh_destroy(barcode, cb_ub_h);
      kh_destroy(bc_info_hash_t, cb_h);
      }


      //Cleaning up. 
      sam_itr_destroy(iter);
      // When you're done with dosages, remember to free it to avoid memory leaks
      free(dosages);

      // Output to file if there are reads covering the SNP
      if (ref_count > 0 || alt_count > 0) {
        //printf("++%s\t%d\t%s\t%f\t%d\t%d\n", bcf_hdr_id2name(vcf_hdr, rec->rid), pos, rec->d.id, alf, ref_count, alt_count);
	//        fprintf(stderr,"%s\t%d\t%s\t%d\t%d\n", bcf_hdr_id2name(vcf_hdr, rec->rid), pos, rec->d.id, ref_count, alt_count);
        // You may want to write this information to a file instead of printing it
      }

    }

    // Print barcode DLDA matrix and BC/SNP/UMI information. 

    FILE *fp = fopen(output_filename, "w");
    for(int i = 0; i<nbcs; i++){
      fprintf(fp,"%d\t%d\t%d\n",i,bc_count_snps[i],bc_count_umis[i]);
	}
    fclose(fp);


    fp = fopen("out.dlda.txt", "w");
    for(int i = 0; i<nbcs; i++){
      for(int j = 0;j<nsmpl; j++)
        fprintf(fp,"%f\t",bc_mat_dlda[i*nsmpl + j]);
      fprintf(fp,"%f\n");
    }
    fclose(fp);


    // Clean up
    fprintf(stderr, "Closing everything \n");

    free(bc_count_snps);
    free(bc_count_umis); 
    free(bc_mat_dlda); 


    bcf_destroy(rec);
    bcf_hdr_destroy(vcf_hdr);
    bcf_close(vcf_fp);
    hts_idx_destroy(bam_idx);
    bam_hdr_destroy(bam_hdr);
    sam_close(sam_fp);
    kh_destroy(bc_hash_t, hbc);

    return 0;
}

