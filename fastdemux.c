#include <stdio.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>

#include "khash.h"
#include "ketopt.h"
#include "ksort.h"

#include <math.h>
#include <zlib.h>

#include <omp.h>

#include "common.h"

// Define the hash table for storing barcodes (CB tags) with an index to array
// KHASH_MAP_INIT_STR(barcode,int)  //  decomissioning this

#define SQRT2 1.41421356237
#define BAM_BUFF_SIZE 4
#define OMP_THREADS 12


// Define a structure to hold the value with its index for sorting
typedef struct {
  float value;
  int index;
} pair_t;

// Comparator macro for descending sort
#define pair_gt(a, b) ((a).value > (b).value)
KSORT_INIT(pair, pair_t, pair_gt)

// Define the structure to hold the barcode index, reference allele count, and alternate allele count
typedef struct {
  unsigned long int kmer; //the kmer
  int bc_index;  // index of the barcode as it appears in the barcode file (see below). 
  int ref_count; // reference allele count.   
  int alt_count; // alternate allele count.       
} bc_info_t;

// Define the hash map type for CB, value counting number of UMIs, it could be used to count reads too. 
KHASH_MAP_INIT_INT64(bc_info_hash_t, bc_info_t)


typedef struct {
  unsigned long int kmer; //the kmer
  int bc_index;  // index of the barcode as it appears in the barcode file (see below). 
} bc_anno_t;


// Define the hash table for storing barcodes (CB tags) with an index to array, 
KHASH_MAP_INIT_INT64(bc_hash_t,bc_anno_t) //converting bc to integer, 32 bits should be enough but usign 64 for longer BC, and for the combined. 


// Function to load valid barcodes from a file into a hash table
khash_t(bc_hash_t) *load_barcodes(char *filename) {
    khash_t(bc_hash_t) *hbc = kh_init(bc_hash_t);
    //    FILE *fp = fopen(filename, "r");
    gzFile fp = gzopen(filename, "r");
    char line[256];
    int i=0;
    //    while (fgets(line, sizeof(line), fp)) {
    while (gzgets(fp, line, sizeof(line))) { 
       int ret;
        khint64_t k;
        line[strcspn(line, "\n")] = 0; // Remove newline character
        k = kh_put(bc_hash_t, hbc, packKmer(line,16), &ret); // Need to add logic for different size BC
	assert(ret>0); // Making sure the key is not there. No repeats allowed.
	kh_val(hbc, k).bc_index=i++;  //This is the bc_index for later. 
	kh_val(hbc, k).kmer= packKmer(line,16);
	if(i<10){
	  fprintf(stderr,"%d,%lu,%lu,%s\n",i,k,packKmer(line,16),line);
	}
    }
    //    fclose(fp);
    gzclose(fp);
    return hbc;
}




// Setting for long options. 
static ko_longopt_t longopts[] = {
  { "min-qual-score", ko_required_argument, 301 },
  { "min-base-qual",  ko_required_argument, 302 },
  { "num-threads",  ko_required_argument, 303 },
  { "max-var", ko_required_argument, 304 },
  { "max-snps", ko_required_argument, 305 },
  //  { "output-dlda", ko_no_argument, 306},
  { NULL, 0, 0 } // End marker
};

// Main function to process BAM and VCF files
int main(int argc, char *argv[]) {

    ketopt_t opt = KETOPT_INIT;
    float max_var = 0.1;
    int c, min_qual_score = 20, min_base_qual = 13, num_threads = OMP_THREADS, max_snps=0; // default values
    // int output_dlda =0;

    fprintf(stderr, "fastdemux version %s\n",VERSION);
    while ((c = ketopt(&opt, argc, argv, 1, "q:b:t:", longopts)) >= 0) {
      switch(c){
      case 'q':
      case 301:
	min_qual_score = atoi(opt.arg);
	break;
      case 'b':
      case 302:
	min_base_qual = atoi(opt.arg);
	break;
      case 't':
      case 303:
        num_threads = atoi(opt.arg);
	break;
      case 304:
        max_var = atof(opt.arg);
	break;
      case 305:
        max_snps = atoi(opt.arg);
	break;
      case '?':
	fprintf(stderr,"unknown opt: -%c\n", opt.opt? opt.opt : ':');
	return 1;
      case ':':
	fprintf(stderr,"missing arg: -%c\n", opt.opt? opt.opt : ':');
	return 1;
      default:
	break;
      }
    }

    // Check for required files and other positional arguments if necessary
    if (argc - opt.ind < 4) {
      fprintf(stderr, "Usage: %s [options] <BAM file> <VCF file> <Barcode file> <Output prefix>\n", argv[0]);
      fprintf(stderr, "Options:\n");
      fprintf(stderr, " -q, --min-qual-score NUM  Minimum quality score (default: %d)\n",min_qual_score);
      fprintf(stderr, " -b, --min-base-qual NUM   Minimum base quality (default: %d)\n",min_base_qual);
      fprintf(stderr, " --max-var REAL  Threshold on max 2*q*(1-q) to select more informative SNPs  (default: %f)\n",max_var);	
      fprintf(stderr, " --max-snps NUM  Max number of SNPs  (useful for debug, default: %d All)\n",max_snps);
      fprintf(stderr, " -t, --num-threads NUM  Number of OpenMP threads (default: %d)\n",num_threads);
      return 1;
    }

    fprintf(stderr, "Options settings:\n");
    fprintf(stderr, " -q, --min-qual-score: %d\n",min_qual_score);
    fprintf(stderr, " -b, --min-base-qual: %d\n",min_base_qual);
    fprintf(stderr, " --max-var: %f\n",max_var);	
    fprintf(stderr, " --max-snps %d (0=All)\n",max_snps);
    fprintf(stderr, " -t, --num-threads: %d\n",num_threads);
 

    char *bam_filename = argv[opt.ind];
    char *vcf_filename = argv[opt.ind + 1];
    char *barcode_filename = argv[opt.ind + 2];
    char *output_prefix = argv[opt.ind + 3];
    char filename[256];

    int nbcs=0;
    int i,j;
    

    omp_set_num_threads(num_threads); //Adjust this to env variable.... or option

    // Load barcodes, 
    fprintf(stderr, "Reading the barcode file %s\n", barcode_filename);
    khash_t(bc_hash_t) *hbc = load_barcodes(barcode_filename);
    nbcs = kh_size(hbc);
    fprintf(stderr, "Number of barcodes is: %d\n", nbcs);

    // Open BAM file
    samFile *sam_fp[BAM_BUFF_SIZE * num_threads];
    bam_hdr_t *bam_hdr;  //I only need one hdr in multithreaded. 
    hts_idx_t *bam_idx[BAM_BUFF_SIZE * num_threads];
    
    fprintf(stderr, "Opening BAM file %s\n", bam_filename);

    for(j=0; j < BAM_BUFF_SIZE*num_threads; j++){
      sam_fp[j] = sam_open(bam_filename, "r");
      if (sam_fp[j] == NULL) {
        fprintf(stderr, "Failed to open BAM file %s\n", bam_filename);
        return 1;
      }
      bam_idx[j] = sam_index_load(sam_fp[j], bam_filename); // Load BAM index
      if (bam_idx[j] == NULL) {
        fprintf(stderr, "Failed to open BAM index for %s\n", bam_filename);
        sam_close(sam_fp[j]);
        return 1;
      }
    }
    bam_hdr = sam_hdr_read(sam_fp[0]);

    // multithreaded for sam, does not seem to make it faster, but slower
    //hts_set_threads(sam_fp, 6); // 6 threads?

    // Open VCF file
    fprintf(stderr, "Open VCF file %s\n", vcf_filename);
    htsFile *vcf_fp = bcf_open(vcf_filename, "r");
    if (vcf_fp == NULL) {
        fprintf(stderr, "Failed to open VCF file %s\n", vcf_filename);
        return 1;
    }
    bcf_hdr_t *vcf_hdr = bcf_hdr_read(vcf_fp);

    if (!vcf_hdr) {
      fprintf(stderr, "Error: Could not read header from VCF file %s\n", vcf_filename);
      return 1;
    }

    // Process VCF and BAM files...
    // Assuming bam_hdr, sam_fp, vcf_fp, vcf_hdr are already defined and initialized
    int nsmpl = bcf_hdr_nsamples(vcf_hdr); // Number of samples

    // Iterate over sample names in the VCF header
    fprintf(stderr,"There are %d samples in the vcf file: \n",bcf_hdr_nsamples(vcf_hdr));
    for (i = 0; i < nsmpl; i++) {
      const char *sample_name = vcf_hdr->samples[i]; // Get sample name
      fprintf(stderr,"%s\t",sample_name);
    }
    fprintf(stderr,"\n");

    // Verify chromosome names match between bam/vcf
    // This assumes that are ordered the same way between the two files and it could be ext 
    fprintf(stderr,"There are %d contigs in the bam file\n",bam_hdr->n_targets);
    fprintf(stderr,"There are %d contigs in the vcf file\n",vcf_hdr->n[BCF_DT_CTG]); 
    fprintf(stderr,"Verifying chromosome order and naming for the vcf file is the same as in the bam file: ");
    for (i = 0; i < vcf_hdr->n[BCF_DT_CTG]; ++i) {
      const char *bam_chr = bam_hdr->target_name[i];
      const char *vcf_chr = bcf_hdr_id2name(vcf_hdr, i); //vcf_hdr->id[BCF_DT_CTG][i].key;
      // Find corresponding chromosome ID in BAM header
      //      int bam_contig_id = bam_get_tid(bam_header, vcf_contig);
      if (strcmp(bam_chr, vcf_chr) != 0) {
	fprintf(stderr,"NOT OK\n");
	fprintf(stderr,"Mismatch found: BAM[%d] = %s, VCF[%d] = %s\n", i, bam_chr, i, vcf_chr);
	return 1;
      }
    }
    fprintf(stderr,"OK\n");


    int *bc_count_snps = (int *)malloc(nbcs * sizeof(int));  // Counts SNPs per barcode with >0 reads/UMIs.
    assert(bc_count_snps != NULL);
    int *bc_count_umis = (int *)malloc(nbcs * sizeof(int)); // Counts UMIs in SNPs per barcode. 
    assert(bc_count_umis != NULL);
    float *bc_mat_dlda = (float *)malloc(nsmpl * nbcs * sizeof(float)); //Matrix nbcs rows and nsmpl columns, to store DLDA results.  
    assert(bc_mat_dlda != NULL);

    unsigned long int *bc_kmer = (unsigned long int *) malloc(nbcs * sizeof(unsigned long int)); 
    assert(bc_kmer != NULL);

    //To calculate correlation on genetic dosage between samples. We could consider calculating the kinship coefficient instead. 
    int npairs = nsmpl * (nsmpl+1) / 2; 
    float *corpairs = (float *)malloc(npairs * sizeof(float));
    float **cormat = (float **)malloc(nsmpl * sizeof(float *)); // To be defined as upper triangular matrix... Keep always j<=i
    assert(corpairs != NULL);
    assert(cormat != NULL);
    int kk=0;
    for(int ii=0;ii<nsmpl;ii++){
      cormat[ii]=&corpairs[kk];
      for(int jj=0;jj<=ii;jj++)
	corpairs[kk++] = 0; 
    }

    //Traverse the hash to save the kmer of CB?
    for (khint64_t k = kh_begin(hbc); k != kh_end(hbc); ++k)  // traverse
      if (kh_exist(hbc, k))            // test if a bucket contains data
	bc_kmer[kh_val(hbc,k).bc_index]= kh_val(hbc, k).kmer;

    // Verifying first 10 barcodes for debugging. 
    char kmerstr[20];
    for(i=0;i<10;i++){
      unPackKmer(bc_kmer[i],16,kmerstr);
      fprintf(stderr,"%d,%lu,%s\n",i,bc_kmer[i],kmerstr);
    }

    int snpcnt=0;
    int snpcntused=0;

    //Initialize snp/umi count per barcode, and dlda matrix, and correlation matrix.
    for (i = 0; i < nbcs; i++) bc_count_snps[i]=0;
    for (i = 0; i < nbcs; i++) bc_count_umis[i]=0;
    for (i = 0; i < nbcs*nsmpl; i++) bc_mat_dlda[i]=0.0;


    //buffer rec for multithreaded. 
    bcf1_t *rec[BAM_BUFF_SIZE*num_threads];
    for(j=0;j<BAM_BUFF_SIZE*num_threads; j++)
      rec[j] = bcf_init();


    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    //////////          VCF LOOP                    ////////////
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////

    //while (bcf_read(vcf_fp, vcf_hdr, rec) == 0) { // Iterate over VCF records
    while(1){ // Iterate over VCF records;
      // Buffering records for multithreaded.
      for(j = 0; j < BAM_BUFF_SIZE*num_threads; j++){
	snpcnt++;
	if(bcf_read(vcf_fp, vcf_hdr, rec[j]) != 0) break; 
      }
      if(j==0) // If no more things in the buffer. 
	break;

      // Useful for debugging a few snps; e.g. using --max-snps 10000 option.  
      if(max_snps>0 && snpcnt>max_snps) break;  

      //schedule(dynamic)
#pragma omp parallel for 
      for(int jj=0;jj<j;jj++){

	if (bcf_unpack(rec[jj], BCF_UN_STR) < 0) continue; // Unpack the current VCF record

	// Filter bi-allelic SNPs
	if (rec[jj]->n_allele != 2 || strlen(rec[jj]->d.allele[0]) != 1 || strlen(rec[jj]->d.allele[1]) != 1) continue;

	// Retrieve SNP position and alleles
	int pos = rec[jj]->pos + 1; // Convert to 1-based position
	char ref_allele = rec[jj]->d.allele[0][0];
	char alt_allele = rec[jj]->d.allele[1][0];

	// Initialize allele counts
	int ref_count = 0, alt_count = 0;

	//Extract dosages from vcf. 
	float *dosages=NULL;
	float alf=1E-5; //Instead of 0 I put this so variance does not become to small and avoid division by zero. 
	int nds_arr = 0;

	if (bcf_get_format_float(vcf_hdr, rec[jj], "DS", &dosages, &nds_arr) < 0) continue;
	
	for (int ii = 0; ii < nsmpl; ii++) {
	  //	fprintf(stderr,"%f ",dosages[i]);
	  alf +=dosages[ii];
	}
	alf=alf/(nsmpl*2);

	if(snpcnt%1000==0 && jj==0) // verbose progress every 1000 SNPs approximately. 
	  fprintf(stderr, "Processing %d: %d %s %d %c %c %f\n", snpcnt,rec[jj]->rid, bcf_hdr_id2name(vcf_hdr, rec[jj]->rid), rec[jj]->pos, ref_allele, alt_allele,alf);      

	if((2 * alf * (1 - alf) >= max_var) && (max_var > 0)){  // low freq. alleles are more informative. 
	  free(dosages);
	  continue; 
	}
	

	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	//////////        Iterating  bam records        ////////////
	////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////
	// Setup pileup looping over bam records overlapping SNP position. 
	// May need to convert rid to bamid. if not ordered the same, but for now I'm checking if same exact order.  
	hts_itr_t *iter = sam_itr_queryi(bam_idx[jj], rec[jj]->rid, rec[jj]->pos, rec[jj]->pos + 1);
	if (iter != NULL){
	  //           fprintf(stderr,"iter not null\n");

	  //Moved this inside iter, as no need to calculate DLDA weights if iter is NULL. 
	  //using dosage array to store the DLDA weights. 
	  float *weights=(float *)malloc(nsmpl * sizeof(float));
	  assert(weights!=NULL);
	  for (int ii = 0; ii < nsmpl; ii++)
	    weights[ii] = (dosages[ii] - 2.0 * alf) / (SQRT2 * alf * (1-alf));


	  // Assuming `barcodes` is a khash_t(barcode)* containing valid CB tags
	  // and `reads_hash` is a khash_t(barcode)* used to track unique CB+UB combinations
	  khash_t(bc_info_hash_t) *cb_h = kh_init(bc_info_hash_t); // Hash table for counting reads for CB with unique UMI
	  khash_t(bc_hash_t) *cb_ub_h = kh_init(bc_hash_t); // Hash table for counting unique CB+UB combos as STRING think to convert to integer as well

	  bam1_t *b = bam_init1(); // Initialize a container for a BAM record.

	  while (sam_itr_next(sam_fp[jj], iter, b) >= 0) {

	    // Check the mapping quality
	    if (b->core.qual < min_qual_score) continue; // Skip if the alignment quality is below 10. TODO: make it a parameter. 

	    // Extract CB tag (cell barcode)
	    uint8_t *cb_ptr = bam_aux_get(b, "CB");
	    if (cb_ptr==NULL) continue; // Skip if no CB tag
	    char* cb = bam_aux2Z(cb_ptr); // Convert to string
	    khint64_t k;

	    // Follow the CIGAR and get index to the basepair in the read
	    int read_index = find_read_index_for_ref_pos(b, rec[jj]->pos);
	    if (read_index < 0) continue; // Skip if the position is not covered by the read


	    uint8_t *quals = bam_get_qual(b);
	    if (quals[0] == 0xff) continue;  
	    if (quals[read_index] < min_base_qual) continue;  //Minimum base quality to consider (lower BQ will be skipped)

	    unsigned long int kmer=packKmer(cb,16);
	    // Check if the CB tag is in the list of valid barcodes
	    khint64_t kcb = kh_get(bc_hash_t, hbc, kmer);
	    if (kcb == kh_end(hbc)) continue; // Skip if CB not found among valid barcodes 
  
	    // Extract UB tag (unique molecular identifier)  
	    uint8_t *ub_ptr = bam_aux_get(b, "UB");  // Maybe parameter for this. 
	    // Consider how to handle missing UB tag. if ub_ptr==NULL // we could use the read first bases for ATAC. 
	    char ub[100]; // Assuming UB won't exceed this length
	    if (ub_ptr) strcpy(ub, bam_aux2Z(ub_ptr)); // Convert to string if present
	    else strcpy(ub, ""); // Use an empty string if no UB tag
    
	    // Construct a unique key for CB+UB combination.
	    // This can be handled more efficiently by pushing CB with << k and just converting UB to kmer and adding up. 
	    char cb_ub_key[200]; // Ensure this is large enough to hold CB+UB+delimiter  
	    sprintf(cb_ub_key, "%s-%s", cb, ub);

	    // Check if we have already processed this CB
	    int ret=0; 
	    
	    k = kh_put(bc_info_hash_t, cb_h, kmer , &ret); // Add it to the hash table
	    assert(ret>=0);
	    if(ret>0){ 
	      assert(kmer==kh_val(hbc,kcb).kmer); // otherwise hash maybe corrupted for the key.  
	      kh_val(cb_h,k).bc_index=kh_val(hbc,kcb).bc_index;
	      kh_val(cb_h,k).ref_count=0;
	      kh_val(cb_h,k).alt_count=0;
	      kh_val(cb_h,k).kmer=kmer;
	    }

	    // Check if we have already processed this CB+UB combination... May need to pass paramters for length. 
	    int k2 = kh_put(bc_hash_t, cb_ub_h, packKmer(cb_ub_key,strlen(cb_ub_key)), &ret);
	    if (ret>0) { // This is a new, unique CB+UB combination
	      //  I could do majority voting for multiple cominations, but I will pick the first one for now.
	      uint8_t *seq = bam_get_seq(b); // Get the sequence
	      // Get base at the read's position corresponding to the SNP
	      char base = seq_nt16_str[bam_seqi(seq, read_index)]; 
	      // Assuming ref_allele and alt_allele are char variables holding the reference and alternate alleles respectively
	      if (base == ref_allele) { 
		ref_count++;
		kh_val(cb_h,k).ref_count++; 
	      }
	      else if (base == alt_allele) { 
		alt_count++;
		kh_val(cb_h,k).alt_count++; 
	      }
	    }// if ret
	  }//while(bam)
	  bam_destroy1(b); // Clean up BAM record container


	  //////////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////////////////
	  // Traversing the barcode hash table with the allele counts for the current SNP //
	  //////////////////////////////////////////////////////////////////////////////////
	  //////////////////////////////////////////////////////////////////////////////////
	  //Traversing the barcode hash table with the allele counts for the current SNP. 
	  for (khint64_t k = kh_begin(cb_h); k != kh_end(cb_h); ++k)  // traverse
	    if (kh_exist(cb_h, k))            // test if a bucket contains data
	      if ((kh_val(cb_h,k).ref_count > 0) || (kh_val(cb_h,k).alt_count > 0)) {
		// If we want to write a pileup. 
		// You may want to write this information to a file instead of printing it. This will allow for ASE calculations down the line. 
		// But it maybe better to just do it on heterozygote genotypes as a second pass once we determine the (bc,individual) pairs. 
		//    printf("--%s\t%d\t%s\t%f\t%d\t%d\t%d\t%d\t%d\n", bcf_hdr_id2name(vcf_hdr, rec->rid), pos, rec->d.id, alf, ref_count, alt_count, 
		// 	   kh_val(cb_h,k).bc_index, kh_val(cb_h,k).ref_count,kh_val(cb_h,k).alt_count);

		//assert(kh_val(cb_h,k).bc_index<nbcs); // may not be necessary any more. as it should work. 

		// Here we do the fastdemux dlda calculations using dosage and store results in a matrix.		
#pragma omp critical
		bc_count_snps[kh_val(cb_h,k).bc_index]++; // Increase number of SNPs seen for this barcode. 
		int N=kh_val(cb_h,k).ref_count + kh_val(cb_h,k).alt_count;
#pragma omp critical
		bc_count_umis[kh_val(cb_h,k).bc_index]+= N; // Increase number of UMIs in SNPs seen for this barcode. 		
		//float val = kh_val(cb_h,k).alt_count - alf * (float)N;
		float val = kh_val(cb_h,k).alt_count/(float)N - alf; // This may work better. 
		//Loop across individuals for a barcode. 
		for(int ii = 0; ii < nsmpl; ++ii)
		  //The weight could be pre-calculated based on alf and dosage for each sample. dosage becomes the wight
		  // weights: dosage[i] = ((dosages[i] - 2 * alf) / (SQRT2 * alf * (1-alf)))
#pragma omp critical
		  bc_mat_dlda[ (kh_val(cb_h,k).bc_index * nsmpl) + ii ] += weights[ii] * val;
	      }// If reads. 
      // Clear the reads hash table for the next SNP
	  kh_destroy(bc_hash_t, cb_ub_h);
	  kh_destroy(bc_info_hash_t, cb_h);
	  free(weights);
	}//iter if


	if ((ref_count > 0 || alt_count > 0) && (alf > 0.01)) {
	  //	  int kk=0;
#pragma omp critical
	  {
	    for(int ii=0; ii < nsmpl; ii++)
	      for(int jj=0; jj <= ii; jj++)
		cormat[ii][jj] += ((dosages[ii] - 2.0 *alf) * (dosages[jj] - 2.0 *alf)) / (2 * alf * (1-alf));
	      //  corpairs[kk++] += ((dosages[ii] - 2.0 *alf) * (dosages[jj] - 2.0 *alf)) / (2 * alf * (1-alf));
	    snpcntused++;
	  }//omp critical
	}

	//Cleaning up. 
	sam_itr_destroy(iter);
	// When you're done with dosages, remember to free it to avoid memory leaks
	free(dosages);

	//Not needed.. 
       	// Output to file if there are reads covering the SNP
	//	if (ref_count > 0 || alt_count > 0) {
	  //printf("++%s\t%d\t%s\t%f\t%d\t%d\n", bcf_hdr_id2name(vcf_hdr, rec->rid), pos, rec->d.id, alf, ref_count, alt_count);
	  //        fprintf(stderr,"%s\t%d\t%s\t%d\t%d\n", bcf_hdr_id2name(vcf_hdr, rec->rid), pos, rec->d.id, ref_count, alt_count);
	  // You may want to write this information to a file instead of printing it
	//	}
      }// for omp rec vcf
    }//while(1) bcf

    fprintf(stderr,"Processed %d SNPs\n",snpcnt);
    fprintf(stderr,"Processed %d SNPs with >0 umis\n",snpcntused);

    strcpy(filename, output_prefix); // Copy the prefix into filename
    strcat(filename, ".corr.txt.gz"); // Append the extension to filename 
    fprintf(stderr," Writing ouput file:  %s\n",filename);
    gzFile fp = gzopen(filename, "w");
    gzprintf(fp,"Sample1\tSample2\tCorrelation\n");
    for(int ii=0;ii<nsmpl;ii++)
      for(int jj=0;jj<=ii;jj++){
	cormat[ii][jj] /= (float)snpcntused;
	gzprintf(fp,"%s\t%s\t%f\n",vcf_hdr->samples[ii],vcf_hdr->samples[jj],cormat[ii][jj]);
	if(cormat[ii][jj]>0.15 && (ii!=jj))
	  fprintf(stderr,"WARNING: Correlation: %f, %s, %s, %d, %d\n",cormat[ii][jj],vcf_hdr->samples[ii],vcf_hdr->samples[jj],ii,jj); 
      }
    gzclose(fp);
    




    int *bc_smpl_count = (int *) malloc(nsmpl * sizeof(int));
    for(i=0 ; i<nsmpl; i++)
      bc_smpl_count[i]=0;

    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    /////    Print BC/SNP/UMI information and results     //////
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    strcpy(filename, output_prefix); // Copy the prefix into filename
    strcat(filename, ".info.txt.gz"); // Append the extension to filename 
    fprintf(stderr," Writing ouput file:  %s\n",filename);
    fp = gzopen(filename, "w");
    //Make header
    //AAAAAAAAAATCGACC-1      339     202     431     1       0.690469        0.690469        AL-001  0.146772        AL-190  0.041283        AL-064 
    //    gzprintf(fp,"\t%d\t%d\t%f\t%s\t%f\t%s\t%f\t%s\n",i,bc_count_snps[i],bc_count_umis[i]
    gzprintf(fp,"BARCODE\tbcnum\tNsnp\tNumi\tDropType\tKletScore\tBestScore\tBestSample\tSecondBestScore\tSecondBestSample\tThirdBestScore\tThirdBestSample\n");
    for(i = 0; i<nbcs; i++){
      pair_t order[nsmpl];
      float kletval,kletvalmax; 
      int kletindex; 
      char kmerstring[32]; //long enough to hold the CB barcode. 
      
      // Not reporting droplets with 0 UMIs. 
      if(bc_count_umis[i]<1)
	continue;

      for (j = 0; j < nsmpl; j++) {
        order[j].index = j ;  // Store the index
        order[j].value = bc_mat_dlda[ i * nsmpl + j ] / bc_count_snps[i];  // Store the corresponding value
      }

      // Sort using 
      ks_mergesort(pair, nsmpl, order, 0);
      
      // Count droplet assigned to individual.
      if(bc_count_umis[i]>=100 && bc_count_snps[i]>=100)
	bc_smpl_count[order[0].index]++;

      // Correct for relatedness....
      int ii=0, jj=0;
      float rhomult=1;
      if(order[0].index < order[1].index){
	ii = order[1].index;
	jj = order[0].index;
      }else{
	ii = order[0].index;
	jj = order[1].index;	
      }
      if(cormat[ii][jj]>0.1)
	rhomult=1/sqrt(1+cormat[ii][jj]);      
      // Determine if singlet, doublet or M-let. 
      kletindex=2; 
      kletvalmax=0; //order[0].value; // Maybe start w/ 0
      kletval=order[0].value; 
      for (j = 1; j < nsmpl; j++) {
	kletval += order[j].value; 
	float aux=rhomult/sqrt(j+1)*kletval;  // add correlation value, but I have to be able to index it. 
	if(kletvalmax < aux){
	  kletvalmax = aux;
	  kletindex = j;
	}
	if(kletindex>15){
	  kletindex=16;
	  break;
	}
	if(order[j].value < 0)
	  break;
      }
      if(kletvalmax<order[0].value) //Determine if best model is singlet. 
	kletindex=0;

      // Consdier unpacking the kmer unPackKmer(unsigned long int kmer,int k,char *s) but I need the kmer stored in an array....
      // bc_kmer[i] has the kmer packed. 
      unPackKmer(bc_kmer[i],16,kmerstring);
      
      gzprintf(fp,"%s-1\t%d\t%d\t%d\t%d\t%f\t%f\t%s\t%f\t%s\t%f\t%s\n",kmerstring,i,bc_count_snps[i],bc_count_umis[i]
	       ,kletindex+1,kletvalmax
	       ,order[0].value,vcf_hdr->samples[order[0].index]
	       ,order[1].value,vcf_hdr->samples[order[1].index]
	       ,order[2].value,vcf_hdr->samples[order[2].index]
	       );
    }
    gzclose(fp);

    ///////////////////////////////////////////////////
    //// summarize individuals with most droplets  ////
    ///////////////////////////////////////////////////
    float *bc_smpl_pct = (float *) malloc(nsmpl * sizeof(float));
    int totaldroplets=0;
    for(i=0 ; i<nsmpl; i++)
      totaldroplets += bc_smpl_count[i];
    fprintf(stderr,"Invividuals with more >10 droplets that have more than 100 SNPs/UMIs\n");
    for(i=0 ; i<nsmpl; i++){
      bc_smpl_pct[i] = bc_smpl_count[i] / (float)totaldroplets;
      if(bc_smpl_count[i]>10)
	fprintf(stderr,"%s,%d,%f\n",vcf_hdr->samples[i],bc_smpl_count[i],bc_smpl_pct[i]*100.0);
    }
    // Second pass after prior update, would only consider a subset of samples to be more likely. 


    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    /////    2nd PASS Print BC/SNP/UMI results            //////
    ////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////
    strcpy(filename, output_prefix); // Copy the prefix into filename
    strcat(filename, ".info.2nd.txt.gz"); // Append the extension to filename 
    fprintf(stderr," Writing ouput file:  %s\n",filename);
    fp = gzopen(filename, "w");
    //Make header
    //AAAAAAAAAATCGACC-1      339     202     431     1       0.690469        0.690469        AL-001  0.146772        AL-190  0.041283        AL-064 
    //    gzprintf(fp,"\t%d\t%d\t%f\t%s\t%f\t%s\t%f\t%s\n",i,bc_count_snps[i],bc_count_umis[i]
    gzprintf(fp,"BARCODE\tbcnum\tNsnp\tNumi\tDropType\tKletScore\tBestScore\tBestSample\tSecondBestScore\tSecondBestSample\tThirdBestScore\tThirdBestSample\n");
    for(i = 0; i<nbcs; i++){
      pair_t order[nsmpl];
      float kletval,kletvalmax; 
      int kletindex; 
      char kmerstring[32]; //long enough to hold the CB barcode. 
      
      // Not reporting droplets with <100 UMIs or <100SNPs. 
      if(bc_count_umis[i]<100 || bc_count_snps[i]<100)
	continue;
      
      for (j = 0; j < nsmpl; j++) {
        order[j].index = j ;  // Store the index
	if(bc_smpl_pct[j] >= 0.01)
	  order[j].value = bc_mat_dlda[ i * nsmpl + j ] / bc_count_snps[i];  // Store the corresponding value
	else
	  order[j].value = -0.01;
      }
      // Sort using 
      ks_mergesort(pair, nsmpl, order, 0);
 
      // Count droplet assigned to individual.
      //      if(bc_count_umis[i]>=100 && bc_count_snps[i]>=100)
      //	bc_smpl_count[order[0].index]++;
      // Correct for relatedness....
      int ii=0, jj=0;
      float rhomult=1;
      if(order[0].index < order[1].index){
	ii = order[1].index;
	jj = order[0].index;
      }else{
	ii = order[0].index;
	jj = order[1].index;	
      }
      if(cormat[ii][jj]>0.1)
	rhomult=1/sqrt(1+cormat[ii][jj]);      
      // Determine if singlet, doublet or M-let. 
      kletindex=2; 
      kletvalmax=0; //order[0].value; // Maybe start w/ 0
      kletval=order[0].value; 
      for (j = 1; j < nsmpl; j++) {
	kletval += order[j].value; 
	float aux=rhomult/sqrt(j+1)*kletval;  // add correlation value, but I have to be able to index it. 
	if(kletvalmax < aux){
	  kletvalmax = aux;
	  kletindex = j;
	}
	if(kletindex>15){
	  kletindex=16;
	  break;
	}
	if(order[j].value < 0)
	  break;
      }
      if(kletvalmax<order[0].value) //Determine if best model is singlet. 
	kletindex=0;

      // Consdier unpacking the kmer unPackKmer(unsigned long int kmer,int k,char *s) but I need the kmer stored in an array....
      // bc_kmer[i] has the kmer packed. 
      unPackKmer(bc_kmer[i],16,kmerstring);
      
      gzprintf(fp,"%s-1\t%d\t%d\t%d\t%d\t%f\t%f\t%s\t%f\t%s\t%f\t%s\n",kmerstring,i,bc_count_snps[i],bc_count_umis[i]
	       ,kletindex+1,kletvalmax
	       ,order[0].value,vcf_hdr->samples[order[0].index]
	       ,order[1].value,vcf_hdr->samples[order[1].index]
	       ,order[2].value,vcf_hdr->samples[order[2].index]
	       );
    }
    gzclose(fp);


    /////////////////////////////////
    //// Print full DLDA matrix  ////
    /////////////////////////////////
    // maybe not always necessary, maybe use an option to see if we want to print. 
    strcpy(filename, output_prefix); // Copy the prefix into filename
    strcat(filename, ".dlda.txt.gz"); // Append the extension to filename 
    fprintf(stderr," Writing ouput file:  %s\n",filename);
    fp = gzopen(filename, "w");
    // Header, next 4 lines.
    gzprintf(fp,"BARCODE");
    for(j = 0;j<nsmpl; j++)
      gzprintf(fp,"\t%s",vcf_hdr->samples[j]);
    gzprintf(fp,"\n");
    for(i = 0; i<nbcs; i++){
      char kmerstring[32]; //long enough to hold the CB barcode. 
      unPackKmer(bc_kmer[i],16,kmerstring);
      gzprintf(fp,"%s-1",kmerstring);
      for(j = 0;j<nsmpl; j++)
        gzprintf(fp,"\t%f",bc_mat_dlda[i*nsmpl + j]);
     gzprintf(fp,"\n");
    }
    gzclose(fp);
    

    // Final Clean up
    fprintf(stderr, "Closing everything \n");

    free(bc_count_snps);
    free(bc_count_umis); 
    free(bc_mat_dlda); 
    free(bc_kmer);
    free(corpairs);
    free(cormat);
    free(bc_smpl_count);
    free(bc_smpl_pct);

    bcf_close(vcf_fp);
    kh_destroy(bc_hash_t, hbc);


    for(j=0;j<BAM_BUFF_SIZE*num_threads; j++){    
      bcf_destroy(rec[j]); //convert to loop
      hts_idx_destroy(bam_idx[j]);
      sam_close(sam_fp[j]);
    }
    bam_hdr_destroy(bam_hdr);

    return 0;
}

