/* Copyright 2013 University of North Carolina at Chapel Hill.  All rights reserved. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <stack>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/sparse_hash_set>
#include <stdexcept>
#include "abra_NativeAssembler.h"

using namespace std;
using google::sparse_hash_map;
using google::sparse_hash_set;

//#define READ_LENGTH 100
//#define KMER 63
//#define MIN_CONTIG_LENGTH 101
//#define MIN_NODE_FREQUENCY 3
//#define MIN_NODE_FREQUENCY 2
#define MAX_CONTIG_SIZE 50000
#define MAX_READ_LENGTH 1001
//#define MIN_BASE_QUALITY 20
#define INCREASE_MIN_NODE_FREQ_THRESHOLD 25000

#define OK 0
#define TOO_MANY_PATHS_FROM_ROOT -1
#define TOO_MANY_CONTIGS -2
#define STOPPED_ON_REPEAT -3

#define MAX_FREQUENCY 32766

int read_length;
int min_contig_length;
int kmer_size;
int min_node_freq;
int min_base_quality;

struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return (s1 == s2) || (s1 && s2 && strncmp(s1, s2, kmer_size) == 0);
  }
};

struct my_hash
{
	unsigned long operator()(const char* kmer) const
	{
		unsigned long hash = 0;
		int c;

//		while((c = *kmer++))
		for (int i=0; i<kmer_size; i++)
		{
			c = kmer[i];
			/* hash = hash * 33 ^ c */
			hash = ((hash << 5) + hash) ^ c;
		}

		return hash;
	}
};

struct struct_pool {
	struct node_pool* node_pool;
	struct read_pool* read_pool;
};

#define NODES_PER_BLOCK 10000
#define MAX_NODE_BLOCKS 500000
#define READS_PER_BLOCK 10000
#define MAX_READ_BLOCKS 100000

struct node_pool {
	struct node** nodes;
	int block_idx;
	int node_idx;
};

struct read_pool {
	char** reads;
	int block_idx;
	int read_idx;
};

struct node {

	//TODO: Collapse from 8 to 2 bits.  Only store as key.
	char* kmer;
	//TODO: Convert to stl?
	struct linked_node* toNodes;
	struct linked_node* fromNodes;
	char* contributingRead;
	unsigned short frequency;
	char hasMultipleUniqueReads;
	char contributing_strand;
};

struct linked_node {
	struct node* node;
	struct linked_node* next;
};

int compare(const char* s1, const char* s2) {
	return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
}

int compare_kmer(const char* s1, const char* s2) {
	return (s1 == s2) || (s1 && s2 && strncmp(s1, s2, kmer_size) == 0);
}

struct struct_pool* init_pool() {
	struct_pool* pool = (struct struct_pool*) malloc(sizeof(struct_pool));
	pool->node_pool = (struct node_pool*) malloc(sizeof(node_pool));
	// Allocate array of arrays
	pool->node_pool->nodes = (struct node**) malloc(sizeof(struct node*) * MAX_NODE_BLOCKS);
	// Allocate first array of nodes
	pool->node_pool->nodes[0] = (struct node*) malloc(sizeof(struct node) * NODES_PER_BLOCK);
	pool->node_pool->block_idx = 0;
	pool->node_pool->node_idx = 0;

	pool->read_pool = (struct read_pool*) malloc(sizeof(read_pool));
	pool->read_pool->reads = (char**) malloc(sizeof(char*) * MAX_READ_BLOCKS);
	pool->read_pool->reads[0] = (char*) malloc(sizeof(char) * (read_length+1) * READS_PER_BLOCK);
	pool->read_pool->block_idx = 0;
	pool->read_pool->read_idx = 0;

	return pool;
}

char* allocate_read(struct_pool* pool) {
	if (pool->read_pool->block_idx > MAX_READ_BLOCKS) {
		printf("READ BLOCK INDEX TOO BIG!!!!\n");
		exit(-1);
	}

	if (pool->read_pool->read_idx >= READS_PER_BLOCK) {
		pool->read_pool->block_idx++;
		pool->read_pool->read_idx = 0;
		pool->read_pool->reads[pool->read_pool->block_idx] = (char*) malloc(sizeof(char) * (read_length+1) * READS_PER_BLOCK);
	}

	return &pool->read_pool->reads[pool->read_pool->block_idx][pool->read_pool->read_idx++ * (read_length+1)];
}

struct node* allocate_node(struct_pool* pool) {
	if (pool->node_pool->block_idx >= MAX_NODE_BLOCKS) {
		printf("NODE BLOCK INDEX TOO BIG!!!!\n");
		exit(-1);
	}

	if (pool->node_pool->node_idx >= NODES_PER_BLOCK) {
		pool->node_pool->block_idx++;
		pool->node_pool->node_idx = 0;
		pool->node_pool->nodes[pool->node_pool->block_idx] = (struct node*) malloc(sizeof(struct node) * NODES_PER_BLOCK);
	}

	return &pool->node_pool->nodes[pool->node_pool->block_idx][pool->node_pool->node_idx++];
}

struct node* new_node(char* seq, char* contributingRead, struct_pool* pool, int strand) {

//	node* my_node = (node*) malloc(sizeof(node));
	node* my_node = allocate_node(pool);
	memset(my_node, 0, sizeof(node));
	my_node->kmer = seq;
//	strcpy(my_node->contributingRead, contributingRead);
	my_node->contributingRead = contributingRead;
	my_node->frequency = 1;
	my_node->hasMultipleUniqueReads = 0;
	my_node->contributing_strand = (char) strand;
	return my_node;
}

char* get_kmer(int idx, char* sequence, struct struct_pool* pool) {
	return &sequence[idx];
}

int is_node_in_list(struct node* node, struct linked_node* list) {
	struct linked_node* ptr = list;

	while (ptr != NULL) {
		if (compare_kmer(ptr->node->kmer, node->kmer)) {
			return 1;
		}
		ptr = ptr->next;
	}

	return 0;
}

void link_nodes(struct node* from_node, struct node* to_node) {
	if (!is_node_in_list(to_node, from_node->toNodes)) {
		struct linked_node* to_link = (linked_node*) malloc(sizeof(linked_node));
		to_link->node = to_node;
		to_link->next = from_node->toNodes;
		from_node->toNodes = to_link;
	}

	if (!is_node_in_list(from_node, to_node->fromNodes)) {
		struct linked_node* from_link = (linked_node*) malloc(sizeof(linked_node));
		from_link->node = from_node;
		from_link->next = to_node->fromNodes;
		to_node->fromNodes = from_link;
	}
}

void increment_node_freq(struct node* node, char* read_seq, int strand) {
	if (node->frequency < MAX_FREQUENCY-1) {
		node->frequency++;
	}

	if (!(node->hasMultipleUniqueReads) &&
		(!compare(node->contributingRead, read_seq) || node->contributing_strand != (char) strand)) {
		node->hasMultipleUniqueReads = 1;
	}
}

int include_kmer(char* sequence, char*qual, int idx) {

	int include = 1;

	for (int i=idx; i<idx+kmer_size; i++) {
		// Discard kmers with ambiguous bases
		if (sequence[i] == 'N') {
			include = 0;
			break;
		}

		// Discard kmers with low base qualities
		if (qual[i] - '!' < min_base_quality) {
			include = 0;
			break;
		}
	}

	return include;
}

void add_to_graph(char* sequence, sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, struct_pool* pool, char* qual, char* read_info) {

	int strand = atoi(read_info);

	struct node* prev = 0;

	for (int i=0; i<=read_length-kmer_size; i++) {

		if (include_kmer(sequence, qual, i)) {
			char* kmer = get_kmer(i, sequence, pool);
	//		printf("\tkmer: %s\n", kmer);

			struct node* curr = (*nodes)[kmer];

			if (curr == NULL) {
				curr = new_node(kmer, sequence, pool, strand);

				if (curr == NULL) {
					printf("Null node for kmer: %s\n", kmer);
					exit(-1);
				}

				(*nodes)[kmer] = curr;
			} else {
				increment_node_freq(curr, sequence, strand);
			}

			if (prev != NULL) {
				link_nodes(prev, curr);
			}

			prev = curr;
		} else {
			prev = NULL;
		}
	}
}

void build_graph(const char* read_file, sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, struct_pool* pool) {
	FILE *fp = fopen(read_file, "r");
	char read[MAX_READ_LENGTH];
	memset(read, 0, MAX_READ_LENGTH);

	char qual[MAX_READ_LENGTH];
	memset(qual, 0, MAX_READ_LENGTH);

	char read_info[128];

	int line = 0;
	while (fscanf(fp, "%s", read_info) != EOF) {
		fscanf(fp, "%s", read);
		fscanf(fp, "%s", qual);
//		printf("read: %d : %s\n", line++, read);
		if (strcmp(read, "") != 0) {
			char* read_ptr = allocate_read(pool);
			memcpy(read_ptr, read, read_length+1);
			add_to_graph(read_ptr, nodes, pool, qual, read_info);
			line++;

			if ((line % 100000) == 0) {
				printf("Processed %d reads.\n", line);
			}
		}
	}

	printf("Num reads: %d\n", line);
	printf("Num nodes: %d\n", nodes->size());

	fclose(fp);
}

struct linked_node* remove_node_from_list(struct node* node, struct linked_node* list) {
	struct linked_node* node_ptr = list;
	struct linked_node* prev_ptr = NULL;

	char is_found = false;
	while ((node_ptr != NULL) && (!is_found)) {
		if (strncmp(node_ptr->node->kmer, node->kmer, kmer_size) == 0) {
			if (prev_ptr == NULL) {
				// Set head of list to next elem
				list = list->next;
			} else {
				// Remove node from list
				prev_ptr->next = node_ptr->next;
			}

			// Free linked_node
			free(node_ptr);
			is_found = true;
		}

		prev_ptr = node_ptr;
		node_ptr = node_ptr->next;
	}

	return list;
}


void cleanup(struct linked_node* linked_nodes) {
	struct linked_node* ptr = linked_nodes;
	while (ptr != NULL) {
		struct linked_node* next = ptr->next;
		free(ptr);
		ptr = next;
	}
}

void prune_graph(sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, char isUnalignedRegion) {

	int freq = min_node_freq;

	if (!isUnalignedRegion) {
		int increase_freq = nodes->size() / INCREASE_MIN_NODE_FREQ_THRESHOLD;

		if (increase_freq > 0) {
			freq = freq + increase_freq;
			printf("Increased mnf to: %d for nodes size: %d\n", freq, nodes->size());
		}
	}

	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
		         it != nodes->end(); ++it) {

		const char* key = it->first;
		struct node* node = it->second;

		if ((node != NULL) && ((node->frequency < freq) || (!(node->hasMultipleUniqueReads)))) {

			// Remove node from "from" lists
			struct linked_node* to_node = node->toNodes;
			while (to_node != NULL) {
				to_node->node->fromNodes = remove_node_from_list(node, to_node->node->fromNodes);
				to_node = to_node->next;
			}

			// Remove node from "to" lists
			struct linked_node* from_node = node->fromNodes;
			while (from_node != NULL) {
				from_node->node->toNodes = remove_node_from_list(node, from_node->node->toNodes);
				from_node = from_node->next;
			}

			// Remove node from map
			nodes->erase(key);
			cleanup(node->toNodes);
			node->toNodes = NULL;
			cleanup(node->fromNodes);
			node->fromNodes = NULL;

			// Free memory
//			free(node);
		}
	}
}

struct linked_node* identify_root_nodes(sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes) {

	struct linked_node* root_nodes = NULL;
	int count = 0;

	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
	         it != nodes->end(); ++it) {
		struct node* node = it->second;
		if ((node != NULL) && (node->fromNodes == NULL)) {
			struct linked_node* next = root_nodes;
			root_nodes = (linked_node*) malloc(sizeof(linked_node));
			root_nodes->node = node;
			root_nodes->next = next;

//			printf("Root: %s\n", node->seq);
			count++;
		}
	}

	printf("num root nodes: %d\n", count);

	return root_nodes;
}

struct contig {
	char seq[MAX_CONTIG_SIZE];
	int size;
	char is_repeat;
	struct node* curr_node;
	sparse_hash_set<const char*, my_hash, eqstr>* visited_nodes;
};

struct contig* new_contig() {
	struct contig* curr_contig;
	curr_contig = (contig*) malloc(sizeof(contig));
//	printf("seq size: %d\n", sizeof(curr_contig->seq));
	memset(curr_contig->seq, 0, sizeof(curr_contig->seq));
	curr_contig->size = 0;
	curr_contig->is_repeat = 0;
	curr_contig->visited_nodes = new sparse_hash_set<const char*, my_hash, eqstr>();
	return curr_contig;
}

struct contig* copy_contig(struct contig* orig) {
	struct contig* copy = (contig*) malloc(sizeof(contig));
	memcpy(copy->seq, orig->seq, MAX_CONTIG_SIZE);
	strncpy(copy->seq, orig->seq, MAX_CONTIG_SIZE);
	copy->size = orig->size;
	copy->is_repeat = orig->is_repeat;
	copy->visited_nodes = new sparse_hash_set<const char*, my_hash, eqstr>(*orig->visited_nodes);
	return copy;
}

void free_contig(struct contig* contig) {
	delete contig->visited_nodes;
	memset(contig, 0, sizeof(contig));
	free(contig);
}

char is_node_visited(struct contig* contig, struct node* node) {
	 sparse_hash_set<const char*, my_hash, eqstr>::const_iterator it = contig->visited_nodes->find(node->kmer);
	 return it != contig->visited_nodes->end();
}

void output_contig(struct contig* contig, int& contig_count, FILE* fp, const char* prefix) {
	if (strlen(contig->seq) >= min_contig_length) {
		if (contig->is_repeat) {
			fprintf(fp, ">%s_%d_repeat\n%s\n", prefix, contig_count++, contig->seq);
		} else {
			fprintf(fp, ">%s_%d\n%s\n", prefix, contig_count++, contig->seq);
		}
	}
}

//#define OK 0
//#define TOO_MANY_PATHS_FROM_ROOT -1
//#define TOO_MANY_CONTIGS -2
//#define STOPPED_ON_REPEAT -3

int build_contigs(
		struct node* root,
		int& contig_count,
		FILE* fp,
		const char* prefix,
		int max_paths_from_root,
		int max_contigs,
		char stop_on_repeat,
		char shadow_mode) {

	int status = OK;
	stack<contig*> contigs;
	struct contig* root_contig = new_contig();
	root_contig->curr_node = root;
	contigs.push(root_contig);

	int paths_from_root = 1;

	while ((contigs.size() > 0) && (status == OK)) {
		// Get contig from stack
		struct contig* contig = contigs.top();

		if (is_node_visited(contig, contig->curr_node)) {
			// We've encountered a repeat
			contig->is_repeat = 1;
			if ((!shadow_mode) && (!stop_on_repeat)) {
				output_contig(contig, contig_count, fp, prefix);
			}
			contigs.pop();
			free_contig(contig);
			if (stop_on_repeat) {
				status = STOPPED_ON_REPEAT;
			}
		}
		else if (contig->curr_node->toNodes == NULL) {
			// We've reached the end of the contig.
			// Append entire current node.
			memcpy(&(contig->seq[contig->size]), contig->curr_node->kmer, kmer_size);

			// Now, write the contig
			if (!shadow_mode) {
				output_contig(contig, contig_count, fp, prefix);
			}
			contigs.pop();
			free_contig(contig);
		}
		else {
			// Append first base from current node
			contig->seq[contig->size++] = contig->curr_node->kmer[0];
			if (contig->size >= MAX_CONTIG_SIZE) {
				char kmer[1024];
				memset(kmer, 0, 1024);
				strncpy(kmer, contig->curr_node->kmer, kmer_size);
				printf("Max contig size exceeded at node: %s\n", kmer);
				exit(-1);
			}

			contig->visited_nodes->insert(contig->curr_node->kmer);

			// Move current contig to next "to" node.
			struct linked_node* to_linked_node = contig->curr_node->toNodes;
			contig->curr_node = to_linked_node->node;
			paths_from_root++;

			// If there are multiple "to" nodes, branch the contig and push on stack
			to_linked_node = to_linked_node->next;
			while (to_linked_node != NULL) {
				//TODO: Do not clone contig for first node.
				struct contig* contig_branch = copy_contig(contig);
//				printf("orig size: %d, copy size: %d\n", contig->visited_nodes->size(), contig_branch->visited_nodes->size());
				contig_branch->curr_node = to_linked_node->node;
				contigs.push(contig_branch);
				to_linked_node = to_linked_node->next;
				paths_from_root++;
			}
		}

		if (contig_count >= max_contigs) {
			status = TOO_MANY_CONTIGS;
		}

		if (paths_from_root >= max_paths_from_root) {
			status = TOO_MANY_PATHS_FROM_ROOT;
		}
	}

	// Cleanup stranded contigs in case processing stopped.
	while (contigs.size() > 0) {
		struct contig* contig = contigs.top();
		contigs.pop();
		free_contig(contig);
	}

	return status;
}

void cleanup(sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, struct struct_pool* pool) {

	// Free linked lists
	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
	         it != nodes->end(); ++it) {
		struct node* node = it->second;

		if (node != NULL) {
			cleanup(node->toNodes);
			cleanup(node->fromNodes);
		}
	}

	// Free nodes and keys
//	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
//	         it != nodes->end(); ++it) {
//		char* key = (char*) it->first;
//		struct node* node = it->second;
//
////		if (node != NULL) {
////			free(node);
////		}
//
//		if (key != NULL) {
//			free(key);
//		}
//	}
//

	for (int i=0; i<=pool->node_pool->block_idx; i++) {
		free(pool->node_pool->nodes[i]);
	}

	free(pool->node_pool->nodes);
	free(pool->node_pool);

	for (int i=0; i<=pool->read_pool->block_idx; i++) {
		free(pool->read_pool->reads[i]);
	}

	free(pool->read_pool->reads);
	free(pool->read_pool);

	free(pool);
}

int assemble(const char* input,
			  const char* output,
			  const char* prefix,
			  int truncate_on_repeat,
			  int max_contigs,
			  int max_paths_from_root,
			  int input_read_length,
			  int input_kmer_size) {

	read_length = input_read_length;

	min_contig_length = read_length + 1;

	//TODO: Parameterize mcl - shorter for unaligned region?
/*
	if (truncate_on_repeat) {
		min_contig_length = read_length + 1;
	} else {
		min_contig_length = 150;
	}
*/

	kmer_size = input_kmer_size;

	struct struct_pool* pool = init_pool();
	sparse_hash_map<const char*, struct node*, my_hash, eqstr>* nodes = new sparse_hash_map<const char*, struct node*, my_hash, eqstr>();

	long startTime = time(NULL);
	printf("Assembling: %s -> %s\n", input, output);
	nodes->set_deleted_key(NULL);

	printf("Building graph\n");
	fflush(stdout);
	build_graph(input, nodes, pool);
	printf("Pruning graph\n");
	fflush(stdout);

	//TODO: Set this explicitly
	char isUnalignedRegion = !truncate_on_repeat;
	prune_graph(nodes, isUnalignedRegion);

	printf("Pruning graph done\n");
	fflush(stdout);

	struct linked_node* root_nodes = identify_root_nodes(nodes);

	int contig_count = 0;
	char truncate_output = 0;

	FILE *fp = fopen(output, "w");
	while (root_nodes != NULL) {

		int shadow_count = 0;

		// Run in shadow mode first
		int status = build_contigs(root_nodes->node, shadow_count, fp, prefix, max_paths_from_root, max_contigs, truncate_on_repeat, true);

		if (status == OK) {
			// Now output the contigs
			build_contigs(root_nodes->node, contig_count, fp, prefix, max_paths_from_root, max_contigs, truncate_on_repeat, false);
		}

		switch(status) {
			case TOO_MANY_CONTIGS:
				printf("TOO_MANY_CONTIGS: %s\n", prefix);
				contig_count = 0;
				break;
			case STOPPED_ON_REPEAT:
				printf("STOPPED_ON_REPEAT: %s\n", prefix);
				contig_count = 0;
				break;
			case TOO_MANY_PATHS_FROM_ROOT:
				char kmer[1024];
				memset(kmer, 0, 1024);
				strncpy(kmer, root_nodes->node->kmer, kmer_size);
				printf("TOO_MANY_PATHS_FROM_ROOT: %s - %s\n", prefix, kmer);
				break;
		}

		// If too many contigs or abort due to repeat, break out of loop and truncate output.
		if ((status == TOO_MANY_CONTIGS) || (status == STOPPED_ON_REPEAT)) {
			truncate_output = 1;
			break;
		}

		root_nodes = root_nodes->next;
	}
	fclose(fp);


	if (truncate_output) {
		FILE *fp = fopen(output, "w");
		fclose(fp);
	}

	cleanup(nodes, pool);

	delete nodes;

	long stopTime = time(NULL);
	printf("Done assembling(%ld): %s -> %s\n", (stopTime-startTime), input, output);

	return contig_count;
}

extern "C"
 JNIEXPORT jint JNICALL Java_abra_NativeAssembler_assemble
   (JNIEnv *env, jobject obj, jstring j_input, jstring j_output, jstring j_prefix,
    jint j_truncate_on_output, jint j_max_contigs, jint j_max_paths_from_root,
    jint j_read_length, jint j_kmer_size, jint j_min_node_freq, jint j_min_base_quality)
 {

     //Get the native string from javaString
     //const char *nativeString = env->GetStringUTFChars(javaString, 0);
	const char* input  = env->GetStringUTFChars(j_input, 0);
	const char* output = env->GetStringUTFChars(j_output, 0);
	const char* prefix = env->GetStringUTFChars(j_prefix, 0);
	int truncate_on_output = j_truncate_on_output;
	int max_contigs = j_max_contigs;
	int max_paths_from_root = j_max_paths_from_root;
	int read_length = j_read_length;
	int kmer_size = j_kmer_size;
	min_node_freq = j_min_node_freq;
	min_base_quality = j_min_base_quality;

	printf("Abra JNI entry point v0.58\n");

	printf("input: %s\n", input);
	printf("output: %s\n", output);
	printf("prefix: %s\n", prefix);
	printf("truncate_on_output: %d\n", truncate_on_output);
	printf("max_contigs: %d\n", max_contigs);
	printf("max_paths_from_root: %d\n", max_paths_from_root);
	printf("read_length: %d\n", read_length);
	printf("kmer_size: %d\n", kmer_size);
	printf("min node freq: %d\n", min_node_freq);
	printf("min base quality: %d\n", min_base_quality);

	int ret = assemble(input, output, prefix, truncate_on_output, max_contigs, max_paths_from_root, read_length, kmer_size);

     //DON'T FORGET THIS LINE!!!
    env->ReleaseStringUTFChars(j_input, input);
    env->ReleaseStringUTFChars(j_output, output);
    env->ReleaseStringUTFChars(j_prefix, prefix);

    fflush(stdout);
    return ret;
 }


/*
//extern "C"
JNIEXPORT void JNICALL Java_abra_NativeAssembler_assemble
  (JNIEnv * env, jobject obj, jstring j_input, jstring j_output, jstring j_prefix)
// JNIEXPORT void JNICALL Java_abra_NativeAssembler_assemble
//   (JNIEnv *env, jobject obj, jstring j_input, jstring j_output, jstring j_prefix)
 {
     //Get the native string from javaString
     //const char *nativeString = env->GetStringUTFChars(javaString, 0);
	const char* input  = env->GetStringUTFChars(j_input, 0);
	const char* output = env->GetStringUTFChars(j_output, 0);
	const char* prefix = env->GetStringUTFChars(j_prefix, 0);

	printf("input: %s\n", input);
	printf("output: %s\n", output);
	printf("prefix: %s\n", prefix);

     //DON'T FORGET THIS LINE!!!
    env->ReleaseStringUTFChars(j_input, input);
    env->ReleaseStringUTFChars(j_output, output);
    env->ReleaseStringUTFChars(j_prefix, prefix);
 }
 */

int main(int argc, char* argv[]) {


        min_node_freq = 2;
        min_base_quality = 5;

        assemble(
                "/datastore/nextgenout4/seqware-analysis/lmose/platinum/test/mtest.reads",
                "/datastore/nextgenout4/seqware-analysis/lmose/platinum/test/mtest.fa",
                "foo",
                false,
                500000,
                5000,
                101,
                53);


/*
	assemble(
		"/home/lmose/code/abra/src/main/c/sim83/unaligned.bam.reads",
		"/home/lmose/code/abra/src/main/c/unaligned_c.fa",
		"foo",
		false,
		50000,
		5000);
*/


	/*
	assemble(
		"/home/lmose/code/abra/src/main/c/sim83_reads_filtered.txt",
		"/home/lmose/code/abra/src/main/c/sim83_c.fasta",
		"foo",
		false,
		50000,
		5000);
		*/

	/*
	assemble(
		"/home/lmose/code/abra/src/main/c/sim83/sm/reads.txt",
		"/home/lmose/code/abra/src/main/c/sim83/sm/reads.fa",
		"foo",
		false,
		50000,
		5000);
		*/

/*
	printf("node size: %d\n", sizeof(node));
	assemble(
		"/home/lmose/code/abra/src/main/c/dump/unaligned.bam.reads",
		"/home/lmose/code/abra/src/main/c/dump/unaligned.fa",
		"foo",
		false,
		500000,
		5000,
		100,
		63);
*/

//	printf("node size: %d\n", sizeof(node2));

	/*
	assemble(
		argv[1],
		argv[2],
		"foo",
		false,
		50000,
		5000);
*/
	/*
	assemble(
		argv[1],
		argv[2],
		"foo",
		false,
		50000,
		5000);
*/

	/*
	assemble(
		"/home/lmose/code/abra/src/main/c/sim83/250000.txt",
		"/home/lmose/code/abra/src/main/c/sim83/250000.fa",
		"foo",
		false,
		5000000,
		5000);
		*/

/*
	assemble(
		"/home/lmose/code/abra/src/main/c/sim83/500000.txt",
		"/home/lmose/code/abra/src/main/c/sim83/500000.fa",
		"foo",
		false,
		5000000,
		5000);
*/

	/*
	assemble(
		"/home/lmose/code/abra/src/main/c/sim83/filtered.txt",
		"/home/lmose/code/abra/src/main/c/sim83/filtered.fa",
		"foo",
		false,
		5000000,
		5000);
*/

/*
	for (int i=0; i<50; i++) {
		char file[200];
		memset(file, 0, sizeof(file));
		sprintf(file, "%s_%d", "/home/lmose/code/abra/src/main/c/sim83/250000_c_run", i);
		assemble(
			"/home/lmose/code/abra/src/main/c/sim83/250000.txt",
			file,
			"foo",
			false,
			1000000,
			5000);
	}
*/
}
