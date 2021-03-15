#include <errno.h>


int load_fasta_sequences(const char *filename, Sequence **psequences, bool removex) {
	FILE *fd = NULL;
	char *seq = NULL;
	char desc[100];
	char buf[10000];
	char *line = NULL;
	int N = 0;
	int i = 0;
  printf("%s ", filename);
	fd = fopen(filename, "r");
	if (fd == NULL) {
		printf("error while opening file %s\n", filename);
	}
	assert(fd != NULL);

	for (;;) {
		line = fgets(buf, 10000, fd);
		//if (line[0] == 0) continue;
		
		if (line != NULL) {
			if (line[strlen(line)] == 0x0A || line[strlen(line)] == 0x0D)
				line[strlen(line)] = 0x00;
			if (strlen(line) > 1 && (line[strlen(line)-1] == 0x0A || line[strlen(line)-1] == 0x0D))
				line[strlen(line) - 1] = 0x00;
		}
		
		if (line == NULL || line[0] == '>') {
			if (seq != NULL) {
				for (i=0; i<strlen(seq); i++) {
					seq[i] = (char)toupper((int)seq[i]);
					
					if (90==seq[i] || 74==seq[i] || 79==seq[i] || 85==seq[i] || 88==seq[i] || 66==seq[i]) {
						if (removex) {
							/* remove all proteins with occurences of non-natural amino acids:
					    		A CDEFGHI KLMN PQRST VW Y
					    		 B       J    O     U  X Z
					    		 66      74   79    85 8890
							*/
							free(seq);
							seq = NULL;
							break;
						} else {
							seq[i] = 'X';
						}
					}
				}
				if (seq != NULL) {
					N++;
					*psequences = realloc(*psequences, N * sizeof(Sequence));
					strncpy((*psequences + N - 1)->description, desc, 100);
					(*psequences + N - 1)->description[99] = '\0';
					(*psequences + N - 1)->sequence = seq;
					seq = NULL;
				}
			}
			if (line == NULL) {
				break;
			}
			strncpy(desc, line + 1, 100);
			desc[99] = '\0';
			continue;
		}
		
		if (strlen(line) == 0) {
			continue;
		}
						
		if (seq != NULL) {
			seq = realloc(seq, strlen(seq) + strlen(line) + 1);
			strcpy(seq + strlen(seq), line);
		} else {
			seq = calloc(strlen(line) + 1, sizeof(char));
			strcpy(seq, line);
		}
	}
	fclose(fd);
	return N;
}

void print_PSSM(double M[50][26], double threshold) {
	int i=0, j=0, a=0;
	printf("PSSM");
	
	for (j=0; j<20; j++) {
		printf(" %c    ", amino_acids[j]);
	}
	printf("\n");
	for (i=0; i<50; i++) {
		printf("%-2d ", i);
		for (j=0; j<20; j++) {
			a = amino_acids[j] - 'A';
			if (M[i][a] > threshold) {
				printf("% 1.2f ", M[i][a]);
			} else {
				printf("  .   ");
			}
		}
		printf("\n");
	}
	printf("\n");
}


// returns number of loaded matrices
int load_Matrices(const char *filename, Matrix **pmatrices) {
	//char segment[52] = "";
	char *buf = calloc(256, sizeof(char));
	int m=0, j=0;
	int pos=0, profile_format;
	int c[20];
	double b[20];
	int matrix_number = 0;
	FILE *fd = NULL;
	Matrix M;

	if (*pmatrices != NULL) {
		printf("Input error. load_Matrices expects an unallocated pointer to matrices (NULL)\n");
		return 0;
	}

	fd = fopen(filename, "r");
	if (fd == NULL){
	      fprintf (stderr, "Couldn't open file %s; %s\n", filename, strerror (errno));
	}
	assert(fd != NULL);

	while (fgets(buf, 256, fd)) {
		if (1 == sscanf(buf, "PROTOTYPE %d", &profile_format)) {
			if (profile_format != 1 && profile_format != 2) {
				printf("Format error. Wrong profile format version! Expecting version 1 or 2.\n");
				return 0;
			}
		}
		if (strstr(buf, "BEGIN")) {
			// no conditions
		}
		if (strstr(buf, "END")) {
			//printf("Sizeof Matrix = %d\n", m*sizeof(Matrix));
			m++;
			*pmatrices = realloc(*pmatrices, m*sizeof(Matrix));
			memcpy(*pmatrices + m - 1, &M, sizeof(Matrix));
		}
		if (profile_format == 1) {
			if (1 == sscanf(buf, "MATRIX K=%d", &M.K)) {
				M.id = matrix_number++;
				M.L = 50;
			}
			if (21 == sscanf(buf,
			    "%2d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d", 
			     &pos,&c[0],&c[1],&c[2],&c[3],&c[4],&c[5],&c[6],&c[7],&c[8],&c[9],&c[10],&c[11],&c[12],&c[13],&c[14],&c[15],&c[16],&c[17],&c[18],&c[19])) {
		    	
				if (pos<0 || pos>49) {
					printf("Format error: invalid position %d\n", pos);
					continue;
				}
				for (j=20; j--;) {
					if (c[j]<0) {
						printf("Format error: invalid position %d,%d = %d\n", pos, j, c[j]);
					}
					M.freq[pos][amino_acids[j] - 'A'] = (double)(c[j]) / M.K;
				}
			}
		}
		if (profile_format == 2) {
			if (2 == sscanf(buf, "MATRIX ID=%d K=%d", &M.id, &M.K)) {
				// no conditions
			}
			if (21 == sscanf(buf,
			    "%2d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
			     &pos,&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7],&b[8],&b[9],&b[10],&b[11],&b[12],&b[13],&b[14],&b[15],&b[16],&b[17],&b[18],&b[19])) {
		    	
				if (pos<0 || pos>49) {
					printf("Format error: invalid position %d\n", pos);
					continue;
				}
				for (j=20; j--;) {
					if (b[j]<0) {
						printf("Format error: invalid position %d,%d = %lf\n", pos, j, b[j]);
					}
					M.freq[pos][amino_acids[j] - 'A'] = b[j];
				}
			}
		}
	}
	//printf("Finished parsing file with matrices\n");
	fclose(fd);
	free(buf);
	return m;
}

// returns number of loaded matrices
// matrix containing gaps, and counts
int load_VariableMatricesCount(const char *filename, Matrix **pmatrices, double composition[20]) {
	const int profile_format_version = 3;
	//char segment[52] = "";
	char *buf = calloc(256, sizeof(char));
	int m=0, j=0;
	double b[21];
	FILE *fd = NULL;
	Matrix M;

	if (*pmatrices != NULL) {
		printf("Input error. load_Matrices expects an unallocated pointer to matrices (NULL)\n");
		return 0;
	}
	int num_residues = -1;

	fd = fopen(filename, "r");
	assert(fd != NULL);
  int row_no = 0;
	while (fgets(buf, 256, fd)) {
		if (strstr(buf, "BEGIN")) {
			for (int pos=0; pos<50; pos++) {
				for (j=20;j--;) {
					M.freq[pos][amino_acids[j] - 'A'] = composition[j];
				}
			}
		}
		if (strstr(buf, "END")) {
			//printf("Sizeof Matrix = %d\n", m*sizeof(Matrix));
				m++;
				*pmatrices = realloc(*pmatrices, m*sizeof(Matrix));
				memcpy(*pmatrices + m - 1, &M, sizeof(Matrix));

		}
		if (2 == sscanf(buf, "MATRIX K=%d L=%d", &M.K, &M.L)) {
			// no conditions
			M.id = m;
			//printf("Matrix %d %d %d\n", m, M.K, L);
		}
		if (20 == sscanf(buf,
		    "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",
		     &b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7],&b[8],&b[9],&b[10],&b[11],&b[12],&b[13],&b[14],&b[15],&b[16],&b[17],&b[18],&b[19])) {
		  
			if (row_no<0 || row_no>49 || row_no>M.L-1) {
				printf("Format error: invalid position %d (hard borders 0-49) (matrix borders 0-%d) \n", row_no, M.L);
				continue;
			}
//			if (num_residues == -1){
//        num_residues = 0;
//        for (int i=0;i<20;i++){
//          num_residues += (int) b[i];
//        }
//        M.K = num_residues;
//			}
      M.K = 450;
			for (j=0; j<20;j++) {
				if (b[j]<0) {
					printf("Format error: invalid position %d,%d = %lf\n", row_no, j, b[j]);
				}
        M.freq[row_no][amino_acids[j] - 'A'] = b[j];
//				M.freq[row_no][amino_acids[j] - 'A'] = b[j] / (double)M.K;
			}
      row_no++;
		}
	}
	//printf("Finished parsing file with matrices\n");
	fclose(fd);
	free(buf);
	return m;
}
// returns number of loaded matrices
int load_VariableMatricesFreq(const char *filename, Matrix **pmatrices, double composition[20]) {
	const int profile_format_version = 4;
	//char segment[52] = "";
	char *buf = calloc(256, sizeof(char));
	int m=0, j=0;
	int pos=0, tmp;
	double b[20];
	FILE *fd = NULL;
	Matrix M;

	for (pos=0; pos<50; pos++) {
		for (j=26;j--;) {
			M.freq[pos][j] = 0.0;
		}
	}

	if (*pmatrices != NULL) {
		printf("Input error. load_Matrices expects an unallocated pointer to matrices (NULL)\n");
		return 0;
	}

	fd = fopen(filename, "r");
	assert(fd != NULL);

	int check = 0;

	while (fgets(buf, 256, fd)) {
		if (1 == sscanf(buf, "PROFILE %d", &tmp)) {
			if (tmp != profile_format_version) {
				printf("Format error. Wrong profile format version! Expecting version %d.\n", profile_format_version);
				return 0;
			}
			check++;
		}
		if (strstr(buf, "BEGIN")) {
			for (pos=0; pos<50; pos++) {
				for (j=20;j--;) {
					M.freq[pos][amino_acids[j] - 'A'] = composition[j];
				}
			}
			check++;
		}
		if (strstr(buf, "END")) {
			if (check < 3) {
				printf("PROFILE 4: load_VariableMatricesFreq detected an invalid profile in your input! only %d checks passed\n", check);
				return 0;
			}
			check = 0;
			
			//printf("Sizeof Matrix = %d\n", m*sizeof(Matrix));
			*pmatrices = realloc(*pmatrices, (m+1)*sizeof(Matrix));
			memcpy(*pmatrices + m, &M, sizeof(Matrix));
			m++;
		}
		if (3 == sscanf(buf, "MATRIX ID=%d K=%d L=%d", &M.id, &M.K, &M.L)) {
		    	printf("PROFILE FORMAT 4. Read matrix %d, ID=%d K=%d L=%d\n", m, M.id, M.K, M.L);
			// no conditions
			check++;
		}
		if (21 == sscanf(buf,
		    "%4d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		     &pos,&b[0],&b[1],&b[2],&b[3],&b[4],&b[5],&b[6],&b[7],&b[8],&b[9],&b[10],&b[11],&b[12],&b[13],&b[14],&b[15],&b[16],&b[17],&b[18],&b[19])) {
		     
			if (pos<0 || pos>49 || pos>M.L-1) {
				printf("Format error: invalid position %d (hard borders 0-49) (matrix borders 0-%d) \n", pos, M.L);
				continue;
			}
			for (j=20; j--;) {
				if (b[j]<0 || b[j]>1) {
					printf("Format error: invalid position %d,%d = %lf\n", pos, j, b[j]);
				}
				M.freq[pos][amino_acids[j] - 'A'] = b[j];
			}
		}
	}
	//printf("Finished parsing file with matrices\n");
	fclose(fd);
	free(buf);
	return m;
}

void write_Matrix(const char *filename, Matrix *matrices, int n, int id) {
	FILE *fd = NULL;
	const int proto_format_version = 1;
	int m=0, i=0, j=0;
	Matrix *M=NULL;
	char full_filename[256] = "";

	if (n <= 0){
		return;
	}
	

	
	for (m=0; m<n; m++) {
		M = matrices + m;
		
		if (id != M->id)
		    continue;
		
		snprintf(full_filename, 256, "%s.%d", filename, id);
		fd = fopen(full_filename, "w");
		assert(fd != NULL);

		fprintf(fd, "PROTOTYPE %d\n", proto_format_version);
		fprintf(fd, "BEGIN\n");
		fprintf(fd, "SEGMENT %s\n", M->initial_segment);
		fprintf(fd, "MATRIX K=%d N=%d P=%1.8lf S=%f W=%f\n", M->K, M->N, M->p, M->score, M->omega);
	
		fprintf(fd, "50 ");
		for (j=0; j<strlen(amino_acids); j++) {
			fprintf(fd, "    %c ", amino_acids[j]);
		}
		fprintf(fd, "\n");
		for (i=0; i<50; i++) {
			fprintf(fd, "%2d ", i);
			for (j=0; j<strlen(amino_acids); j++) {
				assert(M->freq[i][amino_acids[j] - 'A'] < 100000);
				fprintf(fd, "%5d ", (int)round(M->freq[i][amino_acids[j] - 'A'] * M->K));
			}
			fprintf(fd, "\n");
		}
		fprintf(fd, "END\n");
		fclose(fd);
	}
}

void write_EvalueScore(int matrix_id, double EvalueScore[], int top_E, int E)
{
    FILE *fd = NULL;
    char full_filename[256] = "";
    int i = 0;
    
    snprintf(full_filename, 256, "evalues.%d.tab", matrix_id);
    fd = fopen(full_filename, "w");
    assert(fd != NULL);
    fprintf(fd, "profile\tevalue\tscore\t\n");
    for (i=1; i <= top_E; i++) {
	fprintf(fd, "%d\t%d\t%1.8lf\n", matrix_id, i, EvalueScore[i]);
    }
    fclose(fd);
}

