#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "abra_NativeSemiGlobalAligner.h"

using namespace std;

#define DIR_UP   1
#define DIR_DIAG 2
#define DIR_LEFT 3

#define I 0
#define M 1
#define D 2

#define MAX_REF_LEN 5000
#define MAX_CONTIG_LEN 2000

__thread int matrix[MAX_CONTIG_LEN][MAX_REF_LEN][3];
__thread char bt[MAX_CONTIG_LEN][MAX_REF_LEN][3];

__thread int match = 8;
__thread int mismatch_pen = -32;
__thread int gap_open = -48;
__thread int gap_extend = -1;

void populate(const char* seq1, const char* seq2, int seq1_len, int seq2_len) {

	for (int r=1; r<=seq1_len; r++) {
		matrix[r][0][I] = gap_open + (r*gap_extend);
		matrix[r][0][M] = gap_open + (r*gap_extend);
		matrix[r][0][D] = gap_open + (r*gap_extend);
	}

	for (int c=0; c<=seq2_len; c++) {
		matrix[0][c][I] = gap_open + (c*gap_extend);
		matrix[0][c][M] = 0;
		matrix[0][c][D] = gap_open + (c*gap_extend);
	}

	for (int r=1; r<=seq1_len; r++) {
		for (int c=1; c<=seq2_len; c++) {

			//
			// Insertion (lower) matrix
			int insertExt = matrix[r-1][c][I] + gap_extend;
			int insertOpen = matrix[r-1][c][M] + gap_open;

			if (insertExt >= insertOpen) {
				matrix[r][c][I] = insertExt;
				bt[r][c][I] = DIR_UP;
			} else {
				matrix[r][c][I] = insertOpen;
				bt[r][c][I] = DIR_DIAG;
			}

			// Deletion (upper) matrix
			int deleteExt = matrix[r][c-1][D] + gap_extend;
			int deleteOpen = matrix[r][c-1][M] + gap_open;

			if (deleteExt >= deleteOpen) {
				matrix[r][c][D] = deleteExt;
				bt[r][c][D] = DIR_LEFT;
			} else {
				matrix[r][c][D] = deleteOpen;
				bt[r][c][D] = DIR_DIAG;
			}

			//
			// Match/mismatch (middle) matrix
			int insertClose = matrix[r][c][I];
			int baseMatch = seq1[r-1] == seq2[c-1] ? (matrix[r-1][c-1][M] + match) : (matrix[r-1][c-1][M] + mismatch_pen);
			int deleteClose = matrix[r][c][D];

			if (baseMatch >=insertClose && baseMatch >= deleteClose) {
				matrix[r][c][M] = baseMatch;
				bt[r][c][M] = DIR_DIAG;
			} else if (insertClose >= deleteClose) {
				matrix[r][c][M] = insertClose;
				bt[r][c][M] = DIR_UP;
			} else {
				matrix[r][c][M] = deleteClose;
				bt[r][c][M] = DIR_LEFT;
			}
		}
	}
}

struct cigar_elem {
	char op;
	int len;
};

void update_curr_elem(char op, vector <cigar_elem> & elems) {

	if (elems.size() == 0 || elems.back().op != op) {
		cigar_elem new_elem;
		new_elem.op = op;
		new_elem.len = 1;
		elems.push_back(new_elem);
	} else {
		elems.back().len += 1;
	}
}


void backtrack(const char* seq1, const char* seq2, int seq1_len, int seq2_len, char* result) {
	int best_idx = -1;
	int best_score = -300000000;
	int second_best_score = -300000000;
	int row = seq1_len;

	for (int c=1; c<=seq2_len; c++) {
		if (matrix[row][c][M] > best_score) {
			best_idx = c;
			best_score = matrix[row][c][M];
		} else if (matrix[row][c][M] > second_best_score) {
			second_best_score = matrix[row][c][M];
		}
	}

	int r = seq1_len;
	int c = best_idx;
	int ref_end_idx = c;

	vector<cigar_elem> elems;

	int level = M;

	while (r > 0 && c > 0) {
		char curr_bt = bt[r][c][level];

		if (curr_bt == DIR_DIAG) {
			if (level == M) {
				r -= 1;
				c -= 1;
			} else if (level == I) {
				r -= 1;
			} else if (level == D) {
				c -= 1;
			}

			if (level == M) {
				// If moving back to M level from I or D, skip update.
				update_curr_elem('M', elems);
			}

			level = M;

		} else if (curr_bt == DIR_LEFT) {
			if (level == D) {
				c -= 1;
			} else if (level == M) {
				// noop
			}

			level = D;

			update_curr_elem('D', elems);
		} else if (curr_bt == DIR_UP) {
			if (level == I) {
				r -= 1;
			} else if (level == M) {
				// noop
			}

			level = I;

			update_curr_elem('I', elems);
		} else {
			break;
		}
	}

	int ref_idx = c;

	char cigar[2056];
	int idx = 0;
	for (int i=elems.size()-1; i>=0; i--) {
		snprintf(cigar+idx, 2056-idx, "%d%c", elems[i].len, elems[i].op);
		idx = strlen(cigar);
	}

	sprintf(result, "%d:%d:%d:%d:%s", best_score, second_best_score, ref_idx, ref_end_idx, cigar);
}

void align(const char* seq1, const char* seq2, char* result) {

	int seq1_len = strlen(seq1);
	int seq2_len = strlen(seq2);

	populate(seq1, seq2, seq1_len, seq2_len);
	backtrack(seq1, seq2, seq1_len, seq2_len, result);
}

extern "C"
 JNIEXPORT jstring JNICALL Java_abra_NativeSemiGlobalAligner_align
   (JNIEnv *env, jobject obj, jstring j_seq1, jstring j_seq2, jint j_match, jint j_mismatch,
		   jint j_gap_open, jint j_gap_extend) {

	match = j_match;
	mismatch_pen = j_mismatch;
	gap_open = j_gap_open;
	gap_extend = j_gap_extend;

	const char* seq1  = env->GetStringUTFChars(j_seq1, 0);
	const char* seq2  = env->GetStringUTFChars(j_seq2, 0);

	char result[4098];
	align(seq1, seq2, result);

//	fprintf(stderr, "SGA result: %s\n", result);

	jstring ret = env->NewStringUTF(result);

	env->ReleaseStringUTFChars(j_seq1, seq1);
	env->ReleaseStringUTFChars(j_seq2, seq2);

	return ret;
}

/*
int main(int argc, char* argv[]) {

	const char* ref = "CCAGATCAGCCTAGGCAACATGGTGAAACCCCGTCTCTACCAAAAATAAAAAACTTAGCTGAGCGTGGTGGTGCACGCCTGTAGCCCCAGCTGCTGAGGAGCCTGAGCCCAGGGGGTGGAGGCTGCAGTGAGCCATGATCACACTACTGTACTCCAGCCTAGGTGACAGAGTGAGACCCTGTCTCAAAAAAATAAAAGAAAATAAAAATAAACAAAGAGAGAAGTGGAAGAAGAGGTGGAGTTTTGTATTTATGACTTGAATTTTGTATTCATGACTGGGTTGACACCCCAATCCACTCCATTTTTAGCCTTGAAACATGGCAAACAGTAACCATTAAAAGGATGGAAAAGAGAAGAAGGCATGGGTGGGAAACTGTGCCTCCCATTTTTGTGCATCTTTGTTGCTGTCCTTCCACTATACTGTACCTTTCAGCATTTTGACGGCAACCTGGATTGAGACTCCTGTTTTGCTAATTCCATAAGCTGTTGCGTTCATCACTTTTCCAAAAGCACCTGATCCTAGTACCTTCCCTGCAAAGACAAATGGTGAGTACGTGCATTTTAAAGATTTTCCAATGGAAAAGAAATGCTGCAGAAACATTTGGCACATTCCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCTGTACCATCTGTAGCTGGCTTTCATACCTAAATTGCTTCAGAGATGAAATGATGAGTCAGTTAGGAATAGGCAGTTCTGCAGATAGAGGAAAGAATAATGAATTTTTACCTTTGCTTTTACCTTTTTGTACTTGTGACAAATTAGCAGGGTTAAAACGACAATGAAGAGGAGACAAACACCAATTGTTGCATAGAATGAGATGTTGTCTTGGATGAAAGGGAAGGGGCCTGCAACAAAAGAGTGTCACTCAGCGATGAAACAGAATTCCTGTGTGACATTATAAATAGTGGACAACTCATTATAATCTCTCACATCCTGTTTCAGTAATAATCATTTTCAGTCCTAACAACCACTCTACATATACTCTACTCCCCACAGACAATCAGGCAATGTCCCTGTAAAGGATACATTTCCTCCCTAGAAAATTGCGGATTATTCTCAATCCATTCTTTAAAACCATTTACTAGGGTAAATTTACAAGAATTACATCTGGTCCAGGCACGATGGCTCACGCCTGTAGTCCCAGCACTTTGGGAGGCCAAGATGGGAGGATCACTTGAGTCCAAGAATTAGACACCAGCCCAGGCAACACAGTGAAATCCCGTCTCTAAAAAAATTCAAAAATTAGCTGGGCGTGGTGGCAGGTGCCTGTAATCCCAGCTGCTCGGGAGGCTGAGGCAGGAG";
	const char* seq = "CCATTCTTACCAAACTCTAAATTTTCTCTTGGAAACTCCCATTTGAGATCATATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGCACATTCCATTCTTGCCAAACTCTAGATTTTCTCTTGGAAACTCCCATTTGAGATCACATTCATATTCTCTGAAATCAACGTAGAAGTACTCATTATCTGAGGAGCCGGTCACCCGTACCATCTGTAGC";


	for (int i=0; i<1000; i++) {
		char result[4098];
		align(seq, ref, result);

		printf("result: %s\n", result);
	}
}
*/
