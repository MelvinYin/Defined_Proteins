/**
 *    Search PSSM procedure
 *    Coded in C programming language (C99 standard) and MPI
 *    Copyright (c) 2009 Alexandr Goncearenco
 *    Affiliation: CBU, BCCS, UNIFOB AS, UiB, Berezovsky Group.
 */

/*
INPUT:
    extracted.matrix: matrices in frequency format. the needles
    subject: the haystack. SCOP for example
    search parameters (hardcoded): p-value (or E value recalculated on subject): E = p*N

OUTPUT:
    M id \tab score.float \tab full fasta description of the match

WARNING:
    matrices are numbered according to position in extracted.matrix file

///////////////////////////

    read: matrices into memory from frequency format (2) -> F
    read: subject fasta sequences into memory
    calculate composition from subject sequences
    for each F:
	M = calcPSSM(F, K, composition)
	runP(M, permutation_vector, p-value) return S(p-value)
	runS(M, S(p)):
	    report each match to result file

 DONE
///////////////////////////
*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <errno.h>
//#include <error.h>

#include "PSSM.h"
#include "searchPSSM.h"

#define MASTER 0
#define DATA_TAG 1
#define BREAK_TAG 2

// parameters
//
// with N = gcc add20 M: p= 5^10-6 => E=10. rho = 10. epsilon=3
//
const int delta = 10;		// cutting step for initial segments, in residues
const double omega = 1.0;	// pseudocount scaling factor in PSSM
//const double omega = 2.0;	// pseudocount scaling factor in PSSM
//const double p = 0.00000072;	// initial p-value, lower threshold
//const double p = 0.000005;	// initial p-value, lower threshold
const double rho = 10;		// threshold ratio shuffled/natural number of matches where score > score(p)
const double mu = 0.05;		// upper threshold for p-value
const double eta = 0.01;	// increment of p-value in clustering
const double epsilon = 3;	// max distance(M1, M2) to consider M2 as converged
const int theta = 30;		// max # iterations

//const char *extracted_matrix_filename = "extracted.matrix";
//const char *subject_filename = "subject.fasta";
//const char *results_filename = "search_matches.tab";

void master_thread(Sequence *subject) {
	int size = 0;
	FILE *fd = NULL;
	int message[5];
	int *profile_id;
	int *subject_id;
	int *position;
	int *score; // multiplied by 10000
	int *evalue;
	int break_counter = 0;

	profile_id = &message[0];
	subject_id = &message[1];
	position = &message[2];
	score = &message[3];
	evalue = &message[4];
	
//	fd = fopen(results_filename, "w");
//	assert(fd != NULL);
//
//	fclose(fd);
}

int main (int argc, char *argv[]) {

	//in R: sample(0:49)	
	int random_index[11][50] = 
	{{46,20,41,29,21,26,37,13,17,36,47,8,25,43,9,48,39,24,1,5,34,45,32,18,14,40,38,28,7,6,11,27,0,44,33,19,16,35,23,4,42,49,2,15,30,12,3,10,22,31},
	{26,41,10,42,36,32,44,16,21,5,20,34,37,33,35,0,46,45,19,18,39,9,27,1,7,13,8,14,28,15,4,43,48,12,29,25,24,40,31,30,38,47,23,11,3,22,2,49,17,6},
        {17,10,21,31,40,18,44,39,30,29,33,41,38,34,15,7,5,23,36,35,47,42,25,19,48,11,3,37,20,43,9,46,1,26,28,22,12,0,4,24,27,14,32,8,16,6,45,49,2,13},
	{48,1,42,34,15,25,30,35,38,2,46,43,45,47,33,0,20,11,14,7,21,32,49,23,37,5,13,39,18,36,28,26,9,22,31,41,19,10,24,8,6,17,44,40,16,29,4,3,27,12},
	{35,48,29,36,47,45,40,11,33,24,9,8,32,15,43,4,34,5,21,6,14,30,13,28,2,26,39,19,31,46,37,41,1,18,38,49,10,17,44,23,22,3,27,20,25,16,12,42,7,0},
	{37,2,8,0,19,36,39,20,13,47,22,31,7,48,32,26,11,42,17,41,44,10,5,49,45,40,1,30,38,29,24,18,43,6,15,3,14,35,34,46,4,23,33,12,28,21,16,27,25,9},
	{31,37,39,41,44,12,18,26,47,24,36,16,8,15,43,6,30,49,38,35,14,19,22,9,28,34,48,5,45,32,27,23,25,46,21,42,40,11,17,0,20,29,13,10,7,1,2,4,3,33},
	{2,35,10,39,47,28,31,44,43,29,37,6,46,45,30,1,16,32,5,25,24,18,11,34,26,21,38,14,22,7,20,23,17,42,33,48,0,19,41,49,12,9,40,4,13,8,27,3,36,15},
	{28,30,33,14,11,48,26,46,21,1,35,49,15,43,37,42,24,16,10,13,34,25,44,23,41,12,9,7,4,20,8,5,29,45,39,17,0,3,31,27,19,36,6,22,40,47,2,32,18,38},
	{42,33,40,26,0,4,36,22,37,19,8,17,30,46,31,48,44,29,45,28,15,41,13,47,27,16,10,5,24,3,35,38,21,49,34,43,1,25,39,20,12,2,9,7,14,18,6,32,23,11},
        {27,30,16,49,38,33,11,36,0,32,42,7,26,41,40,44,45,4,21,47,35,12,31,13,28,15,2,19,39,22,3,18,5,34,10,48,43,6,29,24,46,14,37,17,1,25,20,23,8,9}};

	int random_index_30[1][30] = 
	{{28,1,29,24,22,19,26,9,0,13,23,18,10,14,20,25,27,3,15,2,7,12,6,11,17,21,8,16,5,4}};
	
	int i=0, j=0, k=0, l=0, m=0, q=0, f=0;
	double composition[20]; for (j=20;j--;) composition[j] = 0.0;
	unsigned long aa_count[26]; for (j=26;j--;) aa_count[j] = 0;
	double PSSM[50][26];	for (i=50;i--;) for (j=26;j--;) PSSM[i][j] = 0.0;
	double Information[50];	for (i=50;i--;) Information[i] = 0.0;
	double S = 0.0;
		
	char ch = '\0';
	char *protein = NULL;
	
	Matrix *matrix = NULL, *matrices = NULL, **pmatrices = &matrices;
	int matrices_count = 0;
	
	Sequence *subject = NULL, **psubject = &subject;
	int N_subject = 0;  	// # subject or haystack
	int N = 0; 		// # unmasked fragments in all subject sequences
	int N_masked = 0;    	// # masked fragments in all subject sequences
	
	FILE *fd_aln = NULL;
	FILE *fd_scores = NULL;
	FILE *fd_pos_scores = NULL;

	double E = 1.0;
	double p = 0.0;
//	bool output_distributions = false;
	bool output_alignment = true;
	
	double information_threshold = 1.0; // min 1 bit is required on position to do the scoring
	bool scale_score_by_informative_positions = true;
	bool information_weighting = true;
//	bool print_positional_scores_table = false;
	
	int number_of_randomizations = 0;
	int a = 1;
	int file_format = 2;
  
  char *alignment_filename = malloc(10000);
  char *extracted_matrix_filename = malloc(10000);
  char *subject_filename = malloc(10000);
  
//  char *results_filename = malloc(10000);

  strcpy(extracted_matrix_filename, argv[1]);
  strcpy(subject_filename, argv[2]);
  strcpy(alignment_filename, argv[3]);
  int motif_len;
  sscanf(argv[4], "%d", &motif_len);
  srand(time(0));
  printf("Motif_len: %d\n", motif_len);
  int* random_index_arr = malloc(sizeof(int) * motif_len);
  int* taken = malloc(sizeof(int) * motif_len);
  int* backup_arr = malloc(sizeof(int) * motif_len);
  for (i=0;i<motif_len;i++){
    random_index_arr[i] = 0;
    taken[i] = -1;
    backup_arr[i] = -1;
  }

  int to_end = -1;
  for (i=0;i<motif_len;i++){
    int random_i;
    if (i == motif_len - 1){
      for (j=0;j<motif_len;j++){
        taken[j] = j;
      }
      for (j=0;j<motif_len;j++){
        taken[random_index_arr[j]] = -1;
      }
      for (j=0;j<motif_len;j++){
        if (taken[j] != -1){
          random_i = j;
          break;
        }
      }
      to_end = 1;
    }
    
    
    while (to_end == -1){
      random_i = rand() % motif_len;
      if (random_i == i) {
        continue;
      }
      int is_in = 0;
      for (j=0;j<motif_len;j++){
        if (taken[j] == -1){
          taken[j] = random_i;
          break;
        }
        if (taken[j] == random_i){
          is_in = 1;
          break;
        }
      }
      if (is_in == 1){
        continue;
      }
      break;
    }
    random_index_arr[i] = random_i;
  }


	N_subject = load_fasta_sequences(subject_filename, psubject, false);
	
	unsigned long size_of_subject = 0;
	unsigned long size_of_subject_aminoacid = 0;
	bool flag1 = false;
	for (k=0;k<N_subject;k++) {
		for (i=0;(ch = (subject + k)->sequence[i++]);) {
			aa_count[ch-'A']++;
			size_of_subject++;
		}
		for (i=0; i+motif_len <= strlen((subject + k)->sequence); i++) {
			// is masked?
			flag1 = false;
			for (j=motif_len;j--;) {
				q = (subject + k)->sequence[i+j];
				if (90==q || 74==q || 79==q || 85==q || 88==q || 66==q) {
					flag1 = true;
					break;
				}
			}
			if (flag1) {
				N_masked++;
			        continue;
			}
			N++;
		}
	}
	
	for (j=20; j--;) {
		q = amino_acids[j]-'A';
		size_of_subject_aminoacid += aa_count[q];
	}
	for (j=20; j--;) {
		q = amino_acids[j]-'A';
		composition[j] = (double)aa_count[q] / size_of_subject_aminoacid;
	}
//	for (j=20; j--;) {
//		#ifndef NDEBUG
////		printf("%c: %f\n", amino_acids[j], composition[j]);
//		#endif
//	}
		
//	if (file_format == 2 || file_format == 1) {
//		matrices_count = load_Matrices(extracted_matrix_filename, pmatrices);
//	}
  matrices_count = load_VariableMatricesCount(extracted_matrix_filename, pmatrices, composition);
//	if (file_format == 3) {
//		matrices_count = load_VariableMatricesCount(extracted_matrix_filename, pmatrices, composition);
//	}
//	if (file_format == 4) {
//		matrices_count = load_VariableMatricesFreq(extracted_matrix_filename, pmatrices, composition);
//	}
	
	p = E / N;


	#ifndef NDEBUG
	printf("searchPSSM: Matrices: %d, Subject sequences: %d, Masked segments: %d, unmasked: %d. p=%E E=%f\n"
	       "---------------------------------------------------------------\n",
	       matrices_count, N_subject, N_masked, N, p, E);
	#endif


	if (E < 1.0) E=1.0;

	int from = 0;
	int to = matrices_count;

	char scores_filename[50];
	int last_informative_position = 0;
	int informative_positions = 0;
	int last_position = 0;
	
	#ifndef NDEBUG
	printf("DISTRIBUTING WORK [] %d matrices, from=%d, to=%d\n", matrices_count, from, to);
	#endif
				
	for (m=from; m<to; m++) {
		matrix = matrices + m;
		
//		if (output_distributions) {
//			snprintf(scores_filename, 50, "search_scores_%d_shuffled.tab", matrix->id);
//			fd_scores = fopen(scores_filename, "w");
//		}

			for (i=motif_len;i--;) for (j=26;j--;) {
				PSSM[i][j] = 0.0;
			}
			for (i=motif_len;i--;) {
				Information[i] = 0.0;
			}
			
			double sum_aa = 0;
			double fragment_freq = 0.0;
			double pc = 0.0;
			double H_sum = 0.0; // attention: it is in nats, not in bits
			
			//print_PSSM(matrix->freq, 0.0);

			informative_positions = 0;
			last_informative_position = 0;
			last_position = 0;
			for (i=0;i<matrix->L;i++) {
				H_sum = 0.0;
				for (j=20;j--;) {
					q = amino_acids[j] - 'A';
					if (matrix->freq[i][q] > 0.0) {
						H_sum += matrix->freq[i][q] * (log(matrix->freq[i][q])/log(2));
					}
				}
				
				// if information weighting is disabled consider that each position is equal =~ 4.3 bits
				if (information_weighting) {
					Information[i] = log(20)/log(2) + H_sum;
				} else {
					Information[i] = log(20)/log(2);
				}

				// ultimate threshold
				if (Information[i] > 0.3) {
					last_position = i;
				}
				
				if (Information[i] < information_threshold) {
				    	Information[i] = 0.0;
				} else {
				    	last_informative_position = i;
					informative_positions++;
				}
//				printf("%d %f %d\n", i, Information[i], last_informative_position);
			}
			
			double composition_sum_squares = 0.0;
			for (j=20;j--;) {
				composition_sum_squares += composition[j]*composition[j];
			}

			for (j=20;j--;) {
				q = amino_acids[j] - 'A';
			    					
				for (i=matrix->L;i--;) {	
					PSSM[i][q] = Information[i] *
					    log( (
						(matrix->freq[i][q] * matrix->K + omega*composition[j] * composition[j]) /
						    (matrix->K + omega * composition_sum_squares)
					    ) / composition[j] );
					assert(isfinite(PSSM[i][q]));
			    	}
			}
			
			List *list = NULL, *item = NULL, *curr = NULL, *prev = NULL;
			int list_length = 0;
			double minS = 0.0;
			int top_E = (int)round(E);
			double EvalueScore[top_E + 1];
			int e = 0;

			for (e=1;e<=top_E;e++) {
				EvalueScore[e] = 0.0;
			}

			for (f=0; f<N_subject; f++) {
				protein = (subject + f)->sequence;
				
				int protein_length = strlen(protein);
				for (l=0; l+matrix->L <= protein_length; l++) {
					S = 0.0;
					for (i=0;i<matrix->L;i++) {
						// compare matrix->L residues in one fragment from proteome with reshuffled matrix position
						q = protein[l+i];
						if ('X' == q) {
							S += 0;
						} else {

						  S +=  PSSM[random_index_arr[i]][q-'A'];

						}
					}
					if (scale_score_by_informative_positions) {
						S /= informative_positions;
					}
					else {
						S /= matrix->L;
					}
					assert(isfinite(S));
					
//					if (output_distributions && S > -10000000) {
//						fprintf(fd_scores, "%f\n", S);
//					}
					
					if (list_length < top_E || S > minS) {
						curr = list;
						prev = NULL;
						while (curr != NULL && curr->score < S) {
							prev = curr;
							curr = curr->next;
						}
						
						//Dec'7 2010, Alex:
						//This removes duplicate scores from the background distribution of scores (reshuffled profile)
						//Note, this is incorrect to estimate the background with no duplicates
						//especially in case of profile with small number (<4) of score contributing positions
						//
						//if (curr != NULL && curr->score == S) {
						//    continue; // remove duplicate scores
						//}
						
						item = malloc(sizeof(List));
						item->score = S;
						item->next = curr;
						if (list == item->next) {
							list = item;
						}
						if (prev != NULL) {
							prev->next = item;
						}						
						list_length++;
						
						//printf("S=%f, minS=%f, list_length=%d\n", S, minS, list_length);
						
						if (list_length > top_E) {
							//truncate one from list head
							prev = list;
							list = list->next;
							free(prev);
							list_length--;
						}
					}
					minS = list->score;
				}			
			}
			while (list != NULL) {
				EvalueScore[list_length] = list->score;
				prev = list;
				list = list->next;
				free(prev);
				list_length--;
			}
			assert(list_length == 0);
			
//			#ifndef NDEBUG
//			write_EvalueScore(matrix->id, EvalueScore, top_E, (int)round(E));
//			#endif
			
			//minS += 5;
			
			// min score
			#ifndef NDEBUG
			printf("Matrix %d: |E(s>%f)| <= %f \n", matrix->id, EvalueScore[(int)E], E);
			#endif
			
			char match_header[256];
			char match_sequence[51];

			char pos_scores_filename[50];
			

    fd_aln = fopen(alignment_filename, "w");
    assert(fd_aln != NULL);

//			if (print_positional_scores_table) {
//	    			snprintf(pos_scores_filename, 50, "search_pos_scores_%d.txt", matrix->id);
//				fd_pos_scores = fopen(pos_scores_filename, "w");
//				assert(fd_pos_scores != NULL);
//			}
			
//			if (output_distributions) {
//				snprintf(scores_filename, 50, "search_scores_%d_natural.tab", matrix->id);
//				fd_scores = fopen(scores_filename, "w");
//			}
			
			// Use S(p) as threshold on natural matrix matches:
			int protein_length;
			int message[5];
			int real_E = 0;
			double pos_scores[50]; 
			
			for (f=0; f<N_subject; f++) {
				protein = (subject + f)->sequence;
				protein_length = strlen(protein);
				for (l=0; l+matrix->L <= protein_length; l++) {
					S = 0.0;
					for (i=0;i<matrix->L;i++) {
						q = protein[l+i];
						if ('X' == q) {
							S += 0;
							pos_scores[i] = S;
						} else {
							S += PSSM[i][q-'A'];
							pos_scores[i] = PSSM[i][q-'A'];
						}
					}
					if (scale_score_by_informative_positions) {
						S /= informative_positions;
					}
					else {
						S /= matrix->L;
					}
					
//					if (output_distributions && S > -10000000) {
//						fprintf(fd_scores, "%f\n", S);
//					}
					
//					if (print_positional_scores_table) {
//						snprintf(match_header, 120, "%s", (subject + f)->description);
//						fprintf(fd_pos_scores, ">%s\n", match_header);
//						snprintf(match_sequence, last_position+2, "%s", protein + l);
//						fprintf(fd_pos_scores, "%s\n", match_sequence);
//
//						fprintf(fd_pos_scores, "       ");
//						for (i=0;i<50;i++) {
//							fprintf(fd_pos_scores, "\t%2d(%c) ", i+1, protein[l+i]);
//						}
//						fprintf(fd_pos_scores, "\n");
//						fprintf(fd_pos_scores, "%2.4f", S);
//						for (i=0;i<50;i++) {
//							fprintf(fd_pos_scores, "\t%2.3f", pos_scores[i]);
//							pos_scores[i] = 0.0;
//						}
//						fprintf(fd_pos_scores, "\n");
//					}
					// ALGORITHM MODIFICATION HERE (Dec'7 2010, Alex):
					// it was a strict S<EvalueScore[E] condition before
					// which did not work in case of a large tail consisting of the same scores,
					// which can happen if profile only contains a small number of positions (2-3)
					// contributing to the score
					if (S >= EvalueScore[(int)E]) {
					
						for (real_E = (int)E; real_E > 1;) {
						    if (S <= EvalueScore[real_E]) {
							break;
						    } else {
							real_E--;
						    }
						}
						
						message[0] = matrix->id; 	  // profile
						message[1] = f;			  // sequence
						message[2] = l;			  // position
						message[3] = (int)round(S*10000); // score
						message[4] = real_E; 		  // Evalue
						
						if (output_alignment) {
              if (l < (15 - motif_len/2)){
                continue;
              }
							snprintf(match_header, 120, "%s", (subject + f)->description);
							fprintf(fd_aln, ">%s\n", match_header);
							//snprintf(match_sequence, last_informative_position+2, "%s", protein + l);
							snprintf(match_sequence, 31, "%s", protein + l - (15 -
							motif_len/2));
              fprintf(fd_aln, "%d\n", l);
              fprintf(fd_aln, "%s\n", match_sequence);
						}
					}		
				}
			}
			
			if (output_alignment) {
				fclose(fd_aln);
			}
//			if (output_distributions) {
//				fclose(fd_scores);
//			}
//			if (print_positional_scores_table) {
//				fclose(fd_pos_scores);
//			}
	}
		
	return(0);
}


