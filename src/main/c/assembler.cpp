#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <queue>
#include <iostream>
#include <stack>
#include <list>
#include <vector>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/sparse_hash_set>
#include <sparsehash/dense_hash_map>
#include <sparsehash/dense_hash_set>
#include <stdexcept>
#include "abra_NativeAssembler.h"

using namespace std;
using google::dense_hash_map;
using google::dense_hash_set;

//#define READ_LENGTH 100
//#define KMER 63
//#define MIN_CONTIG_LENGTH 101
//#define MIN_NODE_FREQUENCY 3
//#define MIN_NODE_FREQUENCY 2
#define MAX_CONTIG_SIZE 5000
#define MAX_READ_LENGTH 1001
//#define MIN_BASE_QUALITY 20
#define INCREASE_MIN_NODE_FREQ_THRESHOLD 1600

#define MAX_TOTAL_CONTIG_LEN 10000000

#define OK 0
#define TOO_MANY_PATHS_FROM_ROOT -1
#define TOO_MANY_CONTIGS -2
#define STOPPED_ON_REPEAT -3
#define TOO_MANY_NODES -4

#define MAX_FREQUENCY 32766
#define MAX_QUAL_SUM 255

// TODO: This is used to bound qual sum arrays.  Use a memory pool instead for this.
#define MAX_KMER_LEN 201

// TODO: Allocate dynamically
#define MAX_SAMPLES 8

// Kmers containing bases below this threshold are excluded from assembly.
#define MIN_BASE_QUALITY 8

// Minimum edge frequency as percent
// #define MIN_EDGE_FREQUENCY .01

//TODO: Better variable localization
__thread int read_length;
__thread int min_contig_length;
__thread int kmer_size;
__thread int min_node_freq;
__thread int min_base_quality;
__thread double min_edge_ratio;
__thread int debug;
__thread int max_nodes;
__thread int next_node_id = 1;

#define BIG_CONSTANT(x) (x##LLU)

uint64_t MurmurHash64A ( const void * key, int len, uint64_t seed )
{
  const uint64_t m = BIG_CONSTANT(0xc6a4a7935bd1e995);
  const int r = 47;

  uint64_t h = seed ^ (len * m);

  const uint64_t * data = (const uint64_t *)key;
  const uint64_t * end = data + (len/8);

  while(data != end)
  {
    uint64_t k = *data++;

    k *= m;
    k ^= k >> r;
    k *= m;

    h ^= k;
    h *= m;
  }

  const unsigned char * data2 = (const unsigned char*)data;

  switch(len & 7)
  {
  case 7: h ^= uint64_t(data2[6]) << 48;
  case 6: h ^= uint64_t(data2[5]) << 40;
  case 5: h ^= uint64_t(data2[4]) << 32;
  case 4: h ^= uint64_t(data2[3]) << 24;
  case 3: h ^= uint64_t(data2[2]) << 16;
  case 2: h ^= uint64_t(data2[1]) << 8;
  case 1: h ^= uint64_t(data2[0]);
          h *= m;
  };

  h ^= h >> r;
  h *= m;
  h ^= h >> r;

  return h;
}


struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return (s1 == s2) || (s1 && s2 && strncmp(s1, s2, kmer_size) == 0);
  }
};

struct my_hash
{
	uint64_t operator()(const char* kmer) const
	{
		return MurmurHash64A(kmer, kmer_size, 97);
		//return chunk;
	}
};

struct eqint
{
  bool operator()(int i1, int i2) const
  {
    return i1 == i2;
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
	char* seq;
	char kmer_seq[2];
	//TODO: Convert to stl?
	struct linked_node* toNodes;
	struct linked_node* fromNodes;
	char* contributingRead;
	unsigned char qual_sums[MAX_KMER_LEN];
	unsigned short sample_frequency[MAX_SAMPLES];
	int id;
	unsigned short frequency;
	char hasMultipleUniqueReads;
	char contributing_strand;
	char is_condensed;
	char is_filtered;
	char is_root;
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
	struct_pool* pool = (struct struct_pool*) calloc(1, sizeof(struct_pool));
	pool->node_pool = (struct node_pool*) calloc(1, sizeof(node_pool));
	// Allocate array of arrays
	pool->node_pool->nodes = (struct node**) calloc(MAX_NODE_BLOCKS, sizeof(struct node*));
	// Allocate first array of nodes
	pool->node_pool->nodes[0] = (struct node*) calloc(NODES_PER_BLOCK, sizeof(struct node));
	pool->node_pool->block_idx = 0;
	pool->node_pool->node_idx = 0;

	pool->read_pool = (struct read_pool*) calloc(1, sizeof(read_pool));
	pool->read_pool->reads = (char**) calloc(MAX_READ_BLOCKS, sizeof(char*));
	pool->read_pool->reads[0] = (char*) calloc(READS_PER_BLOCK, sizeof(char) * (read_length+1));
	pool->read_pool->block_idx = 0;
	pool->read_pool->read_idx = 0;

	return pool;
}

char* allocate_read(struct_pool* pool) {
	if (pool->read_pool->block_idx > MAX_READ_BLOCKS) {
		fprintf(stderr,"READ BLOCK INDEX TOO BIG!!!!\n");
		exit(-1);
	}

	if (pool->read_pool->read_idx >= READS_PER_BLOCK) {
		pool->read_pool->block_idx++;
		pool->read_pool->read_idx = 0;
		pool->read_pool->reads[pool->read_pool->block_idx] = (char*) calloc(READS_PER_BLOCK, sizeof(char) * (read_length+1));
	}

	return &pool->read_pool->reads[pool->read_pool->block_idx][pool->read_pool->read_idx++ * (read_length+1)];
}

struct node* allocate_node(struct_pool* pool) {
	if (pool->node_pool->block_idx >= MAX_NODE_BLOCKS) {
		fprintf(stderr,"NODE BLOCK INDEX TOO BIG!!!!\n");
		exit(-1);
	}

	if (pool->node_pool->node_idx >= NODES_PER_BLOCK) {
		pool->node_pool->block_idx++;
		pool->node_pool->node_idx = 0;
		pool->node_pool->nodes[pool->node_pool->block_idx] = (struct node*) calloc(NODES_PER_BLOCK, sizeof(struct node));
	}

	return &pool->node_pool->nodes[pool->node_pool->block_idx][pool->node_pool->node_idx++];
}

unsigned char phred33(char ch) {
	return ch - '!';
}

struct node* new_node(char sample_id, char* seq, char* contributingRead, struct_pool* pool, int strand, char* quals) {

//	node* my_node = (node*) malloc(sizeof(node));
	node* my_node = allocate_node(pool);

	my_node->id = next_node_id++;

//	memset(my_node, 0, sizeof(node));
	my_node->kmer = seq;
//	strcpy(my_node->contributingRead, contributingRead);
	my_node->contributingRead = contributingRead;
	my_node->frequency = 1;
	my_node->sample_frequency[sample_id-1] = 1;
	my_node->hasMultipleUniqueReads = 0;
	my_node->contributing_strand = (char) strand;
	for (int i=0; i<kmer_size; i++) {
		my_node->qual_sums[i] = phred33(quals[i]);
	}
	my_node->kmer_seq[0] = my_node->kmer[0];
	return my_node;
}

char* get_kmer(int idx, char* sequence) {
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

void increment_node_freq(char sample_id, struct node* node, char* read_seq, int strand, char* kmer_qual) {
	if (node->frequency < MAX_FREQUENCY-1) {
		node->frequency++;
	}

	if (node->sample_frequency[sample_id-1] < MAX_FREQUENCY-1) {
		node->sample_frequency[sample_id-1] += 1;
	}

	if (!(node->hasMultipleUniqueReads) &&
		(!compare(node->contributingRead, read_seq) || node->contributing_strand != (char) strand)) {
		node->hasMultipleUniqueReads = 1;
	}

	for (int i=0; i<kmer_size; i++) {
		unsigned char phred33_qual = phred33(kmer_qual[i]);
		if ((node->qual_sums[i] + phred33_qual) < MAX_QUAL_SUM) {
			node->qual_sums[i] += phred33_qual;
		} else {
			node->qual_sums[i] = MAX_QUAL_SUM;
		}
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

//		if (qual[i] - '!' < min_base_quality) {
		if (phred33(qual[i]) < MIN_BASE_QUALITY) {
			include = 0;
			break;
		}
	}

	return include;
}

void add_to_graph(char sample_id, char* sequence, dense_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, struct_pool* pool, char* qual, int strand) {

	struct node* prev = 0;

	for (int i=0; i<=read_length-kmer_size; i++) {

		if (include_kmer(sequence, qual, i)) {
			char* kmer = get_kmer(i, sequence);
			char* kmer_qual = get_kmer(i, qual);

			struct node* curr = (*nodes)[kmer];

			if (curr == NULL) {
				curr = new_node(sample_id, kmer, sequence, pool, strand, kmer_qual);

				if (curr == NULL) {
					fprintf(stderr,"Null node for kmer: %s\n", kmer);
					exit(-1);
				}

				(*nodes)[kmer] = curr;
			} else {
				increment_node_freq(sample_id, curr, sequence, strand, kmer_qual);
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

void build_graph2(const char* input, dense_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, struct_pool* pool) {
	int input_len = strlen(input);
	int record_len = read_length*2 + 2;
	int num_records = input_len / record_len;
	int record = 0;
	const char* ptr = input;
	int num_reads = 0;

	while ((record < num_records) && (nodes->size() < max_nodes)) {
		ptr = &(input[record*record_len]);

		char sample_id = ptr[0];

		int strand = 0;

		if (ptr[1] == '0') {
			strand = 0;
		} else if (ptr[1] == '1') {
			strand = 1;
		} else {
			fprintf(stderr,"Initial char in input invalid: %c\n", ptr[1]);
			fprintf(stderr,"ERROR!  INVALID INPUT:\n===========================%s\n===========================\n", input);
			exit(-1);
		}

		// TODO: skip copying the input read.  Downstream code appears to depend
		// upon null terminator.
		char* read_ptr = allocate_read(pool);
//		memset(read_ptr, 0, read_length+1);
		memcpy(read_ptr, &(ptr[2]), read_length);

//		char* read_ptr = (char*) &(ptr[1]);

		char* qual_ptr = (char*) &(ptr[read_length+2]);
		add_to_graph(sample_id, read_ptr, nodes, pool, qual_ptr, strand);
		record++;
	}

//	fprintf(stderr,"Num reads: %d\n", record);
//	fprintf(stderr,"Num nodes: %d\n", nodes->size());
}
/*
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
		if (strcmp(read, "") != 0) {
			char* read_ptr = allocate_read(pool);
			memcpy(read_ptr, read, read_length+1);
			int strand = atoi(read_info);
			add_to_graph(read_ptr, nodes, pool, qual, strand);
			line++;

			if ((line % 100000) == 0) {
				fprintf(stderr,"Processed %d reads.\n", line);
			}
		}
	}

	fprintf(stderr,"Num reads: %d\n", line);
	fprintf(stderr,"Num nodes: %d\n", nodes->size());

	fclose(fp);
}
*/

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

int is_base_quality_good(struct node* node) {
	int is_good = 1;

	for (int i=0; i<kmer_size; i++) {
		if (node->qual_sums[i] < min_base_quality) {
			is_good = 0;
			break;
		}
	}

	return is_good;
}

void remove_node_and_cleanup(const char* key, struct node* node, dense_hash_map<const char*, struct node*, my_hash, eqstr>* nodes) {
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
}

char is_min_edge_ratio_reached(int per_sample_total_freq[], struct node* node) {
	char exceeds_min_ratio = 0;

	for (int i=0; i<MAX_SAMPLES; i++) {
//		fprintf(stderr,"sample: %d, freq: %d, total_freq: %d\n", i, node->sample_frequency[i], per_sample_total_freq[i]);

		if ((per_sample_total_freq[i] > 0) &&
			((double) node->sample_frequency[i] / (double) per_sample_total_freq[i] >= min_edge_ratio)) {

			exceeds_min_ratio = 1;
			break;
		}
	}

	return exceeds_min_ratio;
}

void prune_low_frequency_edges(dense_hash_map<const char*, struct node*, my_hash, eqstr>* nodes) {

	long removed_edge_count = 0;

	for (dense_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
				 it != nodes->end(); ++it) {

		const char* key = it->first;
		node* curr_node = it->second;

		if (curr_node != NULL) {
			////////////////////////////////////////////////
			// Check to node list for low frequency edges
			struct linked_node* to_node = curr_node->toNodes;

			// Calculate total outgoing "edge" frequency
			int to_node_total_freq = 0;

			int per_sample_total_freq[MAX_SAMPLES];
			memset(per_sample_total_freq, 0, sizeof(int)*MAX_SAMPLES);

			while (to_node != NULL) {
				// Using node frequency as proxy for edge frequency here...
				to_node_total_freq = to_node_total_freq + to_node->node->frequency;
				for (int i=0; i<MAX_SAMPLES; i++) {
					per_sample_total_freq[i] += to_node->node->sample_frequency[i];
				}

				to_node = to_node->next;
			}

			// Identify edges to prune
			to_node = curr_node->toNodes;
			vector<node*> to_nodes_to_remove;

			while (to_node != NULL) {
				char exceeds_min_ratio = is_min_edge_ratio_reached(per_sample_total_freq, to_node->node);

				if (!exceeds_min_ratio) {
					to_nodes_to_remove.push_back(to_node->node);
				}

//				if ( ((double) to_node->node->frequency / (double) to_node_total_freq) < min_edge_ratio ) {
//					to_nodes_to_remove.push_back(to_node->node);
//				}

				to_node = to_node->next;
			}

			// Remove edges
			for (vector<node*>::const_iterator iter = to_nodes_to_remove.begin(); iter != to_nodes_to_remove.end(); ++iter ) {
				// Remove edges in each direction
				node* node_to_remove = *(iter);
				node_to_remove->fromNodes = remove_node_from_list(curr_node, node_to_remove->fromNodes);
				curr_node->toNodes = remove_node_from_list(node_to_remove, curr_node->toNodes);
				removed_edge_count += 1;
			}

			////////////////////////////////////////////////
			// Check from node list for low frequency edges
			struct linked_node* from_node = curr_node->fromNodes;

			// Calculate total outgoing "edge" frequency
			int from_node_total_freq = 0;
			memset(per_sample_total_freq, 0, sizeof(int)*MAX_SAMPLES);

			while (from_node != NULL) {
				// Using node frequency as proxy for edge frequency here...
				from_node_total_freq = from_node_total_freq + from_node->node->frequency;

				for (int i=0; i<MAX_SAMPLES; i++) {
					per_sample_total_freq[i] += from_node->node->sample_frequency[i];
				}

				from_node = from_node->next;
			}

			// Identify edges to prune
			from_node = curr_node->fromNodes;
			vector<node*> from_nodes_to_remove;

			while (from_node != NULL) {
				char exceeds_min_ratio = is_min_edge_ratio_reached(per_sample_total_freq, from_node->node);

				if (!exceeds_min_ratio) {
					from_nodes_to_remove.push_back(from_node->node);
				}

//				if ( ((double) from_node->node->frequency / (double) from_node_total_freq) < min_edge_ratio ) {
//					from_nodes_to_remove.push_back(from_node->node);
//				}
				from_node = from_node->next;
			}

			// Remove edges
			for (vector<node*>::const_iterator iter = from_nodes_to_remove.begin(); iter != from_nodes_to_remove.end(); ++iter ) {
				// Remove edges in each direction
				node* node_to_remove = *(iter);
				node_to_remove->toNodes = remove_node_from_list(curr_node, node_to_remove->toNodes);
				curr_node->fromNodes = remove_node_from_list(node_to_remove, curr_node->fromNodes);
				removed_edge_count += 1;
			}
		}
	}

//	fprintf(stderr,"Pruned %ld edges\n", removed_edge_count);
}


void prune_graph(dense_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, char isUnalignedRegion) {

	// First prune kmers that do not reach base quality sum threshold
	for (dense_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
				 it != nodes->end(); ++it) {

		const char* key = it->first;
		struct node* node = it->second;

		if (node != NULL && !is_base_quality_good(node)) {
			remove_node_and_cleanup(key, node, nodes);
		}
	}

//	fprintf(stderr,"Remaining nodes after pruning step 1: %d\n", nodes->size());

	// Now go back through and ensure that each node reaches minimum frequency threshold.
	int freq = min_node_freq;

	/*
	if (!isUnalignedRegion) {
		int increase_freq = nodes->size() / INCREASE_MIN_NODE_FREQ_THRESHOLD;

		if (increase_freq > 0) {
			freq = freq + increase_freq;
//			fprintf(stderr,"Increased mnf to: %d for nodes size: %d\n", freq, nodes->size());
		}
	}
	*/

	if (freq > 1) {
		for (dense_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
					 it != nodes->end(); ++it) {

			const char* key = it->first;
			struct node* node = it->second;

			if ((node != NULL) && ((node->frequency < freq) || (!(node->hasMultipleUniqueReads)))) {
				remove_node_and_cleanup(key, node, nodes);
			}
		}
	}

//	fprintf(stderr,"Remaining nodes after pruning step 2: %d\n", nodes->size());

	prune_low_frequency_edges(nodes);

	// Final pass through cleaning up nodes that are unreachable
	for (dense_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
				 it != nodes->end(); ++it) {

		const char* key = it->first;
		struct node* node = it->second;

		if (node != NULL && node->toNodes == NULL && node->fromNodes == NULL) {
			remove_node_and_cleanup(key, node, nodes);
		}
	}

//	fprintf(stderr,"Remaining nodes after edge pruning: %d\n", nodes->size());
}

void print_kmer(struct node* node) {
    for (int i=0; i<kmer_size; i++) {
            fprintf(stderr,"%c", node->kmer[i]);
    }
}

void print_node(struct node* node) {
        fprintf(stderr,"kmer: ");
        print_kmer(node);
        fprintf(stderr,"\tfrom: ");

        struct linked_node* from = node->fromNodes;
        while (from != NULL) {
        	print_kmer(from->node);
        	fprintf(stderr,",");
        	from = from->next;
        }

        fprintf(stderr,"\tto: ");
        struct linked_node* to = node->toNodes;
        while (to != NULL) {
        	print_kmer(to->node);
        	fprintf(stderr,",");
        	to = to->next;
        }
}

int is_root(struct node* node) {
	int is_root = 0;

	if (node != NULL) {
		if (node->fromNodes == NULL) {
			// No from nodes means this is a root node.
//			if (node->frequency > 1) {
//				// Don't allow singleton to be a root
//				is_root = 1;
//			}
			is_root = 1;
		} else {
			// Identify nodes that point to themselves with no other incoming edges.
			// This will be cleaned up during contig building.
			struct linked_node* from = node->fromNodes;
			if (from->next == NULL && (strncmp(node->kmer, from->node->kmer, kmer_size) == 0)) {
				is_root = 1;
			}
		}
	}

	return is_root;
}

struct linked_node* identify_root_nodes(dense_hash_map<const char*, struct node*, my_hash, eqstr>* nodes) {

	struct linked_node* root_nodes = NULL;
	int count = 0;

	for (dense_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
	         it != nodes->end(); ++it) {
		struct node* node = it->second;

//		if (node != NULL) {
//			print_node(node);
//		}

		if (is_root(node)) {
			struct linked_node* next = root_nodes;
			root_nodes = (linked_node*) malloc(sizeof(linked_node));
			root_nodes->node = node;
			root_nodes->next = next;
			node->is_root = 1;

//			fprintf(stderr,"\tROOT");

			count++;
		}

//		fprintf(stderr,"\n");
	}

//	fprintf(stderr,"num root nodes: %d\n", count);

	return root_nodes;
}

struct contig {
	vector<char*>* fragments;
//	char seq[MAX_CONTIG_SIZE];
	struct node* curr_node;
	dense_hash_set<int, std::tr1::hash<int>, eqint>* visited_nodes;
	double score;
//	int size;
	int real_size;
	char is_repeat;
};

struct contig* new_contig() {
	struct contig* curr_contig;
	curr_contig = (contig*) calloc(1, sizeof(contig));
//	memset(curr_contig->seq, 0, sizeof(curr_contig->seq));
	//curr_contig->size = 0;
	curr_contig->real_size = 0;
	curr_contig->is_repeat = 0;
	curr_contig->visited_nodes = new dense_hash_set<int, std::tr1::hash<int>, eqint>();
	curr_contig->visited_nodes->set_empty_key(0);
	curr_contig->score = 0;
	curr_contig->fragments = new vector<char*>();

	return curr_contig;
}

struct contig* copy_contig(struct contig* orig) {
	struct contig* copy = (contig*) calloc(1, sizeof(contig));

//	memcpy(copy->seq, orig->seq, MAX_CONTIG_SIZE);
//	strncpy(copy->seq, orig->seq, MAX_CONTIG_SIZE);

	// Copy original fragments to new contig
	copy->fragments = new vector<char*>(*(orig->fragments));

//	copy->size = orig->size;
	copy->real_size = orig->real_size;
	copy->is_repeat = orig->is_repeat;
	copy->visited_nodes = new dense_hash_set<int, std::tr1::hash<int>, eqint>(*orig->visited_nodes);
	copy->score = orig->score;
	return copy;
}

void free_contig(struct contig* contig) {
	delete contig->visited_nodes;
	delete contig->fragments;
	free(contig);
}

char is_node_visited(struct contig* contig, struct node* node) {
	dense_hash_set<int, std::tr1::hash<int>, eqint>::const_iterator it = contig->visited_nodes->find(node->id);
	return it != contig->visited_nodes->end();
}

void output_contig(struct contig* contig, int& contig_count, const char* prefix, char* contigs) {


//	if (strlen(contigs) + strlen(contig->seq) > MAX_TOTAL_CONTIG_LEN) {
//		fprintf(stderr,"contig string too long: %s\n", prefix);
//		exit(-1);
//	}

//	if (strlen(contig->seq) >= min_contig_length) {
	if (contig->real_size >= min_contig_length) {
		if (!contig->is_repeat) {

			char buf[MAX_CONTIG_SIZE*2+1];
			buf[0] = '\0';

			for (vector<char*>::iterator it = contig->fragments->begin(); it != contig->fragments->end(); ++it) {
				int to_cat = MAX_CONTIG_SIZE - strlen(buf);
				if (to_cat <= 0) {
					break;
				}
				strncat(buf, *it, to_cat);
			}

			char id_line[1024];
			sprintf(id_line, ">%s_%d_%f\n", prefix, contig_count++, contig->score);
			strcat(contigs, id_line);
			strcat(contigs, buf);
			strcat(contigs, "\n");
		}
	}
}

//#define OK 0
//#define TOO_MANY_PATHS_FROM_ROOT -1
//#define TOO_MANY_CONTIGS -2
//#define STOPPED_ON_REPEAT -3

// Return true if minimum contig score exceeded (where min = lowest of 128 top scoring contigs)
// Update priority queue as needed
char is_contig_score_ok(std::priority_queue<double, std::vector<double>, std::greater<double> > & contig_scores, double score) {
	char is_score_ok = 0;

//	fprintf(stderr, "contig_scores size: %ld\n", contig_scores.size());

	// 128 = max contigs
	// TODO: parameterize
	if (contig_scores.size() == 128) {
		double min_score = contig_scores.top();
		if (score >= min_score) {
			is_score_ok = 1;
		} else {
//			fprintf(stderr, "contig_score filtered: %lf, min: %lf\n", score, min_score);
		}
	} else if (contig_scores.size() < 128) {
		is_score_ok = 1;
	} else {
		fprintf(stderr, "ERROR: Invalid contig score size: %ld\n", contig_scores.size());
		fflush(stderr);
//		exit(-1);
	}

	return is_score_ok;
}

void update_contig_scores(std::priority_queue<double, std::vector<double>, std::greater<double> > & contig_scores, double score) {
	if (contig_scores.size() == 128 && score >= contig_scores.top()) {
		contig_scores.pop();
		contig_scores.push(score);
	}
	else if (contig_scores.size() < 128) {
		contig_scores.push(score);
	}
}

// Append curr node to contig
void append_to_contig(struct contig* contig, vector<char*>& all_contig_fragments, char entire_kmer) {

	if (contig->curr_node->is_condensed) {
		// Add condensed node sequence to fragment vector
		contig->fragments->push_back(contig->curr_node->seq);
		contig->real_size += strlen(contig->curr_node->seq);
	} else {

		if (!entire_kmer) {
			contig->fragments->push_back(contig->curr_node->kmer_seq);
			contig->real_size += 1;
		} else {
			char* fragment = (char*) calloc(kmer_size+1, sizeof(char));
			strncpy(fragment, contig->curr_node->kmer, kmer_size);
			contig->real_size += kmer_size;
			all_contig_fragments.push_back(fragment);
		}
	}
}

int build_contigs(
		struct node* root,
		int& contig_count,
		const char* prefix,
		int max_paths_from_root,
		int max_contigs,
		char stop_on_repeat,
		char shadow_mode,
		char* contig_str,
		std::priority_queue<double, std::vector<double>, std::greater<double> > & contig_scores,
		vector<char*> & all_contig_fragments) {

	int status = OK;
	stack<contig*> contigs;
	stack<contig*> contigs_to_output;
	stack<contig*> popped_contigs;
	struct contig* root_contig = new_contig();
	root_contig->curr_node = root;
	contigs.push(root_contig);

	// Track all contig fragments
	// Initialize to reasonably large number to avoid reallocations
	// TODO: These are just for last kmer in contig.  Use smaller number here?
	int INIT_FRAGMENTS_PER_THREAD = 10000;
	all_contig_fragments.clear();
	all_contig_fragments.reserve(INIT_FRAGMENTS_PER_THREAD);


	int paths_from_root = 1;

	while ((contigs.size() > 0) && (status == OK)) {
		// Get contig from stack
		struct contig* contig = contigs.top();

		if (is_node_visited(contig, contig->curr_node)) {
			// We've encountered a repeat
			contig->is_repeat = 1;
			popped_contigs.push(contig);
			contigs.pop();

			if (stop_on_repeat) {
				status = STOPPED_ON_REPEAT;
			}
		}
		else if (contig->curr_node->toNodes == NULL || contig->real_size >= (MAX_CONTIG_SIZE-1)) {
			// We've reached the end of the contig.
			// Append entire current node.

			append_to_contig(contig, all_contig_fragments, 1);

//			memcpy(&(contig->seq[contig->size]), contig->curr_node->kmer, kmer_size);

			// Now, write the contig
			contigs_to_output.push(contig);
			update_contig_scores(contig_scores, contig->score);

			contigs.pop();
		}
		else {
			// Append first base from current node
			append_to_contig(contig, all_contig_fragments, 0);

			if (contig->real_size >= MAX_CONTIG_SIZE) {
				char kmer[MAX_KMER_LEN];
				memset(kmer, 0, MAX_KMER_LEN);
				strncpy(kmer, contig->curr_node->kmer, kmer_size);
				fprintf(stderr,"Max contig size exceeded at node: %s\n", kmer);

				//TODO: Provide different status
				status = TOO_MANY_CONTIGS;
				break;
			}

			contig->visited_nodes->insert(contig->curr_node->id);

			// Count total edges
			int total_edge_count = 0;
			struct linked_node* to = contig->curr_node->toNodes;

			while (to != NULL) {
				total_edge_count = total_edge_count + to->node->frequency;
				to = to->next;
			}

			// Move current contig to next "to" node.
			struct linked_node* to_linked_node = contig->curr_node->toNodes;
			contig->curr_node = to_linked_node->node;
			paths_from_root++;

			double prev_contig_score = contig->score;

			double log10_total_edge_count = 0;

			// Update contig score if there is fork here.
			// Otherwise, we would just multiply by 1.
			if (to_linked_node->next != NULL) {
				log10_total_edge_count = log10(total_edge_count);
				contig->score = contig->score + log10(contig->curr_node->frequency) - log10_total_edge_count;
			}

			if (!is_contig_score_ok(contig_scores, contig->score)) {
				popped_contigs.push(contig);
				contigs.pop();
			}

			// If there are multiple "to" nodes, branch the contig and push on stack
			to_linked_node = to_linked_node->next;
			while (to_linked_node != NULL) {

				double contig_branch_score = prev_contig_score + log10(to_linked_node->node->frequency) - log10_total_edge_count;

				if (is_contig_score_ok(contig_scores, contig_branch_score)) {
					struct contig* contig_branch = copy_contig(contig);
					contig_branch->curr_node = to_linked_node->node;
					contig_branch->score = contig_branch_score;
					contigs.push(contig_branch);
				}

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

	if (status == OK) {

		while (contigs_to_output.size() > 0) {
			struct contig* contig = contigs_to_output.top();

			if (is_contig_score_ok(contig_scores, contig->score)) {
				output_contig(contig, contig_count, prefix, contig_str);
			}

			contigs_to_output.pop();
			free_contig(contig);
		}
	}

	// Cleanup stranded contigs in case processing stopped.
	while (contigs_to_output.size() > 0) {
		struct contig* contig = contigs_to_output.top();

		contigs_to_output.pop();
		free_contig(contig);
	}

	while (contigs.size() > 0) {
		struct contig* contig = contigs.top();
		contigs.pop();
		free_contig(contig);
	}

	while (popped_contigs.size() > 0) {
		struct contig* contig = popped_contigs.top();
		popped_contigs.pop();
		free_contig(contig);
	}

	for (vector<char*>::iterator it=all_contig_fragments.begin(); it != all_contig_fragments.end(); ++it) {
		free(*it);
	}

	all_contig_fragments.clear();

	return status;
}

void cleanup(dense_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, struct struct_pool* pool) {

	// Free linked lists
	for (dense_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
	         it != nodes->end(); ++it) {
		struct node* node = it->second;

		if (node != NULL) {
			cleanup(node->toNodes);
			cleanup(node->fromNodes);
		}
	}

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

char has_one_incoming_edge(struct node* node) {
	return node->fromNodes != NULL && node->fromNodes->next == NULL;
}

char has_one_outgoing_edge(struct node* node) {
	return node->toNodes != NULL && node->toNodes->next == NULL;
}

char prev_has_multiple_outgoing_edges(struct node* node) {
	char prev_bifurcates = 0;
	if (has_one_incoming_edge(node)) {
		struct node* prev = node->fromNodes->node;
		if (prev->toNodes != NULL && prev->toNodes->next != NULL) {
			prev_bifurcates = 1;
		}
	}

	return prev_bifurcates;
}

#define CONDENSED_SEQ_SIZE 100000
__thread int condensed_seq_idx = 0;
__thread int condensed_seq_cnt = 0;
__thread char* condensed_seq;

char* get_condensed_seq_buf() {

	if (condensed_seq_cnt == 0) {
		condensed_seq_cnt  = 1;
		condensed_seq = (char*) malloc(CONDENSED_SEQ_SIZE * sizeof(char));
	}

	if (condensed_seq_idx + MAX_CONTIG_SIZE+1 >= CONDENSED_SEQ_SIZE*condensed_seq_cnt) {
		condensed_seq_cnt += 1;
		condensed_seq = (char*) realloc(condensed_seq, CONDENSED_SEQ_SIZE*condensed_seq_cnt * sizeof(char));
	}

	return condensed_seq + condensed_seq_idx;
}

// NOTE: From nodes are invalid after this step!!!
void condense_graph(dense_hash_map<const char*, struct node*, my_hash, eqstr>* nodes) {
	for (dense_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
	         it != nodes->end(); ++it) {
		struct node* node = it->second;

		// Starting point 0 or >1 incoming edges or previous node with multiple outgoing edges and curr node has 1 outgoing edge
		if ((!has_one_incoming_edge(node) || prev_has_multiple_outgoing_edges(node)) && has_one_outgoing_edge(node)) {
			struct node* next = node->toNodes->node;

			if (has_one_incoming_edge(next)) {
				struct linked_node* last = next->toNodes;

				int idx = 0;
				char* seq = get_condensed_seq_buf();
				seq[idx++] = node->kmer[0];

				int nodes_condensed = 1;

				while (next != NULL && has_one_incoming_edge(next) && nodes_condensed < MAX_CONTIG_SIZE) {
					last = next->toNodes;

					if (next->toNodes != NULL) {
						seq[idx++] = next->kmer[0];
					} else {
						// End of path, copy entire kmer
						strncpy(&(seq[idx]), next->kmer, kmer_size);
						idx += kmer_size;
					}

					struct node* temp = NULL;

					if (has_one_outgoing_edge(next)) {
						temp = next->toNodes->node;
					} else {
						temp = NULL;
					}

					next->is_filtered = 1;

					next = temp;

					nodes_condensed += 1;
				}

				seq[idx] = '\0';

				// Advance condensed seq buffer idx
				condensed_seq_idx += (strlen(seq) + 1);

				// Update node
				node->seq = seq;
				node->is_condensed = 1;

				// Free original toNodes list
				cleanup(node->toNodes);
				// copy last toNodes to this node

				if (last == NULL) {
					node->toNodes = NULL;
				} else {
					struct linked_node* to_link = (linked_node*) malloc(sizeof(linked_node));
					node->toNodes = to_link;

					while (last != NULL) {
						to_link->node = last->node;
						to_link->next = last->next;
						last = last->next;
						if (last != NULL) {
							to_link->next = (linked_node*) malloc(sizeof(linked_node));
							to_link = to_link->next;
						}
					}
				}
			}
		}
	}
}

void dump_graph(dense_hash_map<const char*, struct node*, my_hash, eqstr>* nodes, const char* filename) {

	fprintf(stderr, "Filename: %s\n", filename);
	FILE* fp = fopen(filename, "w");

	// Output edges
	fprintf(fp, "digraph vdjer {\n//\tEdges\n");
	for (dense_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
				 it != nodes->end(); ++it) {

		const char* key = it->first;
		node* curr_node = it->second;

		if (!curr_node->is_filtered) {
			struct linked_node* to_node = curr_node->toNodes;

			while (to_node != NULL) {
				fprintf(fp, "\tv_%d -> v_%d\n", curr_node->id, to_node->node->id);
				to_node = to_node->next;
			}
		}
	}

	int num_vertices = 0;
	int num_condensed = 0;

	// Output vertices
	fprintf(fp, "//\tVertices\n");
	for (dense_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes->begin();
				 it != nodes->end(); ++it) {

		const char* key = it->first;
		node* curr_node = it->second;

		// Skip orphans
		if (!curr_node->is_filtered) {
			if (curr_node->is_condensed) {
				if (curr_node->is_root) {
					fprintf(fp, "\tv_%d [label=\"%s\",shape=box,color=green]\n", curr_node->id, curr_node->seq);
				} else {
					fprintf(fp, "\tv_%d [label=\"%s\",shape=box,color=blue]\n", curr_node->id, curr_node->seq);
				}
				num_condensed += 1;
			} else {
				char buf[50];
				strncpy(buf, curr_node->kmer, kmer_size);
				buf[kmer_size] = '\0';
				if (curr_node->is_root) {
//					fprintf(fp, "\tv_%d [label=\"%s\",shape=box,color=red]\n", curr_node->id, buf);
					fprintf(fp, "\tv_%d [label=\"%c\",shape=box,color=red]\n", curr_node->id, curr_node->kmer[0]);
				} else {
//					fprintf(fp, "\tv_%d [label=\"%s\",shape=box]\n", curr_node->id, buf);
					fprintf(fp, "\tv_%d [label=\"%c\",shape=box]\n", curr_node->id, curr_node->kmer[0]);
				}
			}

			num_vertices += 1;
		}
	}

	fprintf(fp, "}\n");

	fclose(fp);

	fprintf(stderr, "Num traversable vertices: %d\n", num_vertices);
	fprintf(stderr, "Num condensed vertices: %d\n", num_condensed);
}


char* assemble(const char* input,
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
	dense_hash_map<const char*, struct node*, my_hash, eqstr>* nodes = new dense_hash_map<const char*, struct node*, my_hash, eqstr>();
	nodes->set_empty_key(NULL);
	char* deleted_key = (char*) calloc(kmer_size, sizeof(char));
	nodes->set_deleted_key(deleted_key);

	long startTime = time(NULL);
	if (debug) {
		fprintf(stderr,"Assembling: -> %s\n", output);
	}

	build_graph2(input, nodes, pool);

	int status = -1;

	if (nodes->size() >= max_nodes) {
		status = TOO_MANY_NODES;
		if (debug) {
			fprintf(stderr,"Graph too complex for region: %s\n", prefix);
		}
	}

	//TODO: Set this explicitly
	char isUnalignedRegion = !truncate_on_repeat;
	prune_graph(nodes, isUnalignedRegion);

	struct linked_node* root_nodes = NULL;

	if (status != TOO_MANY_NODES) {
		root_nodes = identify_root_nodes(nodes);
	}

	condense_graph(nodes);

//	char graph_dump[1024];
//	sprintf(graph_dump, "%s.dot", prefix);
//	dump_graph(nodes, graph_dump);

	int contig_count = 0;
	char truncate_output = 0;

	char* contig_str = (char*) calloc(MAX_TOTAL_CONTIG_LEN, sizeof(char));
//	memset(contig_str, 0, MAX_TOTAL_CONTIG_LEN);

	std::priority_queue<double, std::vector<double>, std::greater<double> > contig_scores;
	vector<char*> all_contig_fragments;


	while (root_nodes != NULL) {

		int shadow_count = 0;

		status = build_contigs(root_nodes->node, contig_count, prefix, max_paths_from_root, max_contigs,
				truncate_on_repeat, false, contig_str, contig_scores, all_contig_fragments);

		switch(status) {
			case TOO_MANY_CONTIGS:
				fprintf(stderr,"TOO_MANY_CONTIGS: %s\n", prefix);
				contig_count = 0;
				break;
			case STOPPED_ON_REPEAT:
				if (debug) {
					fprintf(stderr,"STOPPED_ON_REPEAT: %s\n", prefix);
				}
				contig_count = 0;
				break;
			case TOO_MANY_PATHS_FROM_ROOT:
				char kmer[MAX_KMER_LEN];
				memset(kmer, 0, MAX_KMER_LEN);
				strncpy(kmer, root_nodes->node->kmer, kmer_size);
				fprintf(stderr,"TOO_MANY_PATHS_FROM_ROOT: %s - %s\n", prefix, kmer);
				break;
		}

		// If too many contigs or abort due to repeat, break out of loop and truncate output.
		if ((status == TOO_MANY_CONTIGS) || (status == STOPPED_ON_REPEAT)) {
			truncate_output = 1;
			break;
		}

		root_nodes = root_nodes->next;
	}

	cleanup(nodes, pool);

	delete nodes;

	free(deleted_key);

	// Cleanup condensed seq buffer
	if (condensed_seq_cnt > 0) {
		free(condensed_seq);
	}
	condensed_seq_idx = 0;
	condensed_seq_cnt = 0;

	long stopTime = time(NULL);

	if (kmer_size != input_kmer_size) {
		fprintf(stderr,"What!!?? %d : %d\n", kmer_size, input_kmer_size);
	}
	assert(kmer_size == input_kmer_size);
	if (debug) {
		fprintf(stderr,"Done assembling(%ld): %s, %d\n", (stopTime-startTime), output, contig_count);
	}

	if (status == OK || status == TOO_MANY_PATHS_FROM_ROOT) {
		return contig_str;
	} else if (status == STOPPED_ON_REPEAT) {
		strcpy(contig_str, "<REPEAT>");
		return contig_str;
	} else {
		strcpy(contig_str, "<ERROR>");
		return contig_str;
	}
}

extern "C"
 JNIEXPORT jstring JNICALL Java_abra_NativeAssembler_assemble
   (JNIEnv *env, jobject obj, jstring j_input, jstring j_output, jstring j_prefix,
    jint j_truncate_on_output, jint j_max_contigs, jint j_max_paths_from_root,
    jint j_read_length, jint j_kmer_size, jint j_min_node_freq, jint j_min_base_quality,
    jdouble j_min_edge_ratio, jint j_debug, jint j_max_nodes)
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
	min_edge_ratio = j_min_edge_ratio;
	debug = j_debug;
	max_nodes = j_max_nodes;

	if (debug) {
		fprintf(stderr,"Abra JNI entry point, prefix: %s, read_length: %d, kmer_size: %d, min_node_freq: %d, min_base_qual: %d, min_edge_ratio %f, debug: %d, max_nodes: %d\n",
				prefix, read_length, kmer_size, min_node_freq, min_base_quality, min_edge_ratio, debug, max_nodes);
	}

//	printf("input len: %s : %d\n", prefix, strlen(input));
//	printf("output: %s\n", output);
//	printf("prefix: %s\n", prefix);
//	printf("truncate_on_output: %d\n", truncate_on_output);
//	printf("max_contigs: %d\n", max_contigs);
//	printf("max_paths_from_root: %d\n", max_paths_from_root);
//	printf("read_length: %d\n", read_length);
//	printf("kmer_size: %d\n", kmer_size);
//	printf("min node freq: %d\n", min_node_freq);
//	printf("min base quality: %d\n", min_base_quality);
//	printf("min edge ratio: %f\n", min_edge_ratio);

	char* contig_str = assemble(input, output, prefix, truncate_on_output, max_contigs, max_paths_from_root, read_length, kmer_size);
	jstring ret = env->NewStringUTF(contig_str);

     //DON'T FORGET THIS LINE!!!
    env->ReleaseStringUTFChars(j_input, input);
    env->ReleaseStringUTFChars(j_output, output);
    env->ReleaseStringUTFChars(j_prefix, prefix);
    free(contig_str);

    fflush(stdout);

    return ret;
 }

int main(int argc, char* argv[]) {


        min_node_freq = 2;
        min_base_quality = 5;
        min_edge_ratio = .05;

        assemble(
//                "/datastore/nextgenout4/seqware-analysis/lmose/platinum/test/mtest.reads",
//                "/datastore/nextgenout4/seqware-analysis/lmose/platinum/test/mtest.fa",
//                "/datastore/nextgenout4/seqware-analysis/lmose/platinum/test/ftest.reads",
//                "/datastore/nextgenout4/seqware-analysis/lmose/platinum/test/ftest.fa",
		"/datastore/nextgenout4/seqware-analysis/lmose/platinum/long_d/test1.reads",
		"/datastore/nextgenout4/seqware-analysis/lmose/platinum/long_d/test1.fa",
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
