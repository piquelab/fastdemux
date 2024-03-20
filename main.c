#include <stdio.h>
#include <stdlib.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include "khash.h"

// Define the hash table for storing barcodes (CB tags)
KHASH_SET_INIT_STR(barcode)

// Function to load valid barcodes from a file into a hash table
khash_t(barcode) *load_barcodes(char *filename) {
    khash_t(barcode) *h = kh_init(barcode);
    FILE *fp = fopen(filename, "r");
    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        int ret;
        khint_t k;
        line[strcspn(line, "\n")] = 0; // Remove newline character
        k = kh_put(barcode, h, strdup(line), &ret);
	//        if (!ret) kh_del(barcode, h, k); // If the barcode already exists, remove the duplicate
    }
    fclose(fp);
    return h;
}


// Function to find the read sequence index corresponding to a reference position
// based on the CIGAR string of the alignment. Returns -1 if the position is not covered.
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

    // Load barcodes
    fprintf(stderr, "Reading the barcode file %s\n", barcode_filename);
    khash_t(barcode) *barcodes = load_barcodes(barcode_filename);
    fprintf(stderr, "Number of barcodes is: %d\n", kh_size(barcodes));

    // Open BAM file
    fprintf(stderr, "Opening BAM file %s\n", bam_filename);
    samFile *sam_fp = sam_open(bam_filename, "r");
    if (sam_fp == NULL) {
        fprintf(stderr, "Failed to open BAM file %s\n", bam_filename);
        return 1;
    }
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
      //      khint_t k = kh_put(str, sample_labels, strdup(sample_name), &ret); // Insert sample name into hash table
      //      if (!ret) {
      //	free((char *)kh_key(sample_labels, k)); // If the sample name already exists, free the strdup'd name
      //	kh_del(str, sample_labels, k);
      //}
    }
    fprintf(stderr,"\n");

    khash_t(barcode) *reads_hash = kh_init(barcode); // Hash table for counting unique CB+UB
    int ret;

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

      fprintf(stderr, "Processing %s %d %c %c\n",bcf_hdr_id2name(vcf_hdr, rec->rid), rec->pos, ref_allele, alt_allele);
      

      //Extract dosages 
      float *dosages=NULL;
      int nds_arr = 0;

      if (bcf_get_format_float(vcf_hdr, rec, "DS", &dosages, &nds_arr) < 0) continue;

      /*
      for (int i = 0; i < nsmpl; i++) {
	fprintf(stderr,"%f ",dosages[i]);
      }
      fprintf(stderr,"\n");
      */


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

      // When you're done with dosages, remember to free it to avoid memory leaks
      free(dosages);



      // Setup pileup
      hts_itr_t *iter = sam_itr_queryi(bam_idx, rec->rid, rec->pos, rec->pos + 1);
      if (iter == NULL) continue;

      //      fprintf(stderr,"iter not null\n");

      // Assuming `barcodes` is a khash_t(barcode)* containing valid CB tags
      // and `reads_hash` is a khash_t(barcode)* used to track unique CB+UB combinations

      bam1_t *b = bam_init1(); // Initialize a container for a BAM record.

      while (sam_itr_next(sam_fp, iter, b) >= 0) {
	// Extract CB tag (cell barcode)
	uint8_t *cb_ptr = bam_aux_get(b, "CB");
	if (cb_ptr==NULL) continue; // Skip if no CB tag
	char* cb = bam_aux2Z(cb_ptr); // Convert to string

	//        fprintf(stderr,"%s\n",cb);

    
	// Check if the CB tag is in the list of valid barcodes
	khint_t k = kh_get(barcode, barcodes, cb);
	if (k == kh_end(barcodes)) continue; // Skip if CB not found among valid barcodes
    
	// Extract UB tag (unique molecular identifier)
	uint8_t *ub_ptr = bam_aux_get(b, "UB");
	char ub[100]; // Assuming UB won't exceed this length
	if (ub_ptr) strcpy(ub, bam_aux2Z(ub_ptr)); // Convert to string if present
	else strcpy(ub, ""); // Use an empty string if no UB tag
    
	// Construct a unique key for CB+UB combination
	char cb_ub_key[200]; // Ensure this is large enough to hold CB+UB+delimiter
	sprintf(cb_ub_key, "%s-%s", cb, ub);

	//	fprintf(stderr,"%s \n",cb_ub_key);
    
	int read_index = find_read_index_for_ref_pos(b, rec->pos);
	if (read_index < 0) continue; // Skip if the position is not covered by the read


	// Check if we have already processed this CB+UB combination
	k = kh_get(barcode, reads_hash, cb_ub_key);
	if (k == kh_end(reads_hash)) { // This is a new, unique CB+UB combination
	  int ret;
	  k = kh_put(barcode, reads_hash, strdup(cb_ub_key), &ret); // Add it to the hash table
        
	  // Determine if the read matches the reference or alternate allele
	  //	  uint8_t *seq = bam_get_seq(b); // Get the sequence
	  // Potential issue here, if index out of bounds of the read length need to check CIGAR????? #2
	  ///	  int allele_index = bam_seqi(seq, b->core.pos); // Get base at read's starting position
	  //	  char base = seq_nt16_str[allele_index >> 1 >> (allele_index & 1) * 4 & 0xf]; // Convert to base
	  ///      char base = seq_nt16_str[bam_seqi(seq, allele_index)];

	  uint8_t *seq = bam_get_seq(b); // Get the sequence
	  // Get base at the read's position corresponding to the SNP
	  char base = seq_nt16_str[bam_seqi(seq, read_index)]; 

	  // Assuming ref_allele and alt_allele are char variables holding the reference and alternate alleles respectively
	  if (base == ref_allele) ref_count++;
	  else if (base == alt_allele) alt_count++;
	}
      }

      //      fprintf(stderr,"after while \n");

      //      fprintf(stderr,"%s\t%d\t%s\t%d\t%d\n", bcf_hdr_id2name(vcf_hdr, rec->rid), pos, rec->d.id, ref_count, alt_count);


      // After processing, ref_count and alt_count hold the counts of reads matching the reference and alternate alleles,
      // considering only unique CB+UB combinations for the current SNP position.

      bam_destroy1(b); // Clean up BAM record container

      //      bam_plp_destroy(plp);
      sam_itr_destroy(iter);

      // Output to file if there are reads covering the SNP
      if (ref_count > 0 || alt_count > 0) {
        printf("%s\t%d\t%s\t%d\t%d\n", bcf_hdr_id2name(vcf_hdr, rec->rid), pos, rec->d.id, ref_count, alt_count);
	//        fprintf(stderr,"%s\t%d\t%s\t%d\t%d\n", bcf_hdr_id2name(vcf_hdr, rec->rid), pos, rec->d.id, ref_count, alt_count);
        // You may want to write this information to a file instead of printing it
      }


      // Clear the reads hash table for the next SNP
      kh_destroy(barcode, reads_hash);
      reads_hash = kh_init(barcode); // Re-initialize after clearing

    }

    // Clean up
    kh_destroy(barcode, reads_hash);

    // Clean up
    fprintf(stderr, "Closing everything \n");

    bcf_destroy(rec);
    bcf_hdr_destroy(vcf_hdr);
    bcf_close(vcf_fp);
    hts_idx_destroy(bam_idx);
    bam_hdr_destroy(bam_hdr);
    sam_close(sam_fp);
    kh_destroy(barcode, barcodes);

    return 0;
}

