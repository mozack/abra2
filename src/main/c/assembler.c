#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <stack>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/sparse_hash_set>
#include <stdexcept>
using namespace std;
using google::sparse_hash_map;
using google::sparse_hash_set;
//using ext::hash;
//using std::tr1::hash;

#define READ_LENGTH 100
#define KMER 63
#define MIN_CONTIG_LENGTH 101
#define MIN_NODE_FREQUENCY 3

#define MAX_CONTIG_SIZE 10000

struct eqstr
{
  bool operator()(const char* s1, const char* s2) const
  {
    return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
  }
};

struct my_hash
{
	unsigned long operator()(const char* s1) const
	{
		unsigned long hash = 0;
		int c;

		while((c = *s1++))
		{
			/* hash = hash * 33 ^ c */
			hash = ((hash << 5) + hash) ^ c;
		}

		return hash;
	}
};

struct node {

	//TODO: Collapse from 8 to 2 bits.  Only store as key.
	char* seq;
	int frequency;
	struct linked_node* toNodes;
	struct linked_node* fromNodes;
//	sparse_hash_map<const char*, struct node*, hash<const char*>, eqstr> to_nodes;
//	sparse_hash_map<const char*, struct node*, hash<const char*>, eqstr> from_nodes;
	char contributingRead[READ_LENGTH+1];
	char hasMultipleUniqueReads;
};

struct linked_node {
	struct node* node;
	struct linked_node* next;
};

int compare(const char* s1, const char* s2) {
	return (s1 == s2) || (s1 && s2 && strcmp(s1, s2) == 0);
}

sparse_hash_map<const char*, struct node*, my_hash, eqstr> nodes;
//sparse_hash_map<const char*, struct node*, hash<const char*>, eqstr> nodes;

struct node* new_node(char* seq, char* contributingRead) {
	node* my_node = (node*) malloc(sizeof(node));
	memset(my_node, 0, sizeof(node));
	my_node->seq = seq;
//	my_node->contributingRead = contributingRead;
	strcpy(my_node->contributingRead, contributingRead);
	my_node->frequency = 1;
	my_node->hasMultipleUniqueReads = 0;
	return my_node;
}

char* get_kmer(int idx, char* sequence) {
	char* kmer = (char*) malloc(sizeof(char) * KMER+1);
	memset(kmer, 0, KMER+1);

	memcpy(kmer, &sequence[idx], KMER);

	return kmer;
}

int is_node_in_list(struct node* node, struct linked_node* list) {
	struct linked_node* ptr = list;

	while (ptr != NULL) {
		if (compare(ptr->node->seq, node->seq)) {
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

void increment_node_freq(struct node* node, char* read_seq) {
	node->frequency++;

	if (!(node->hasMultipleUniqueReads) && !compare(node->contributingRead, read_seq)) {
		node->hasMultipleUniqueReads = 1;
	}
}

void add_to_graph(char* sequence) {

	struct node* prev = 0;

	for (int i=0; i<=READ_LENGTH-KMER; i++) {

		char* kmer = get_kmer(i, sequence);
//		printf("\tkmer: %s\n", kmer);

		struct node* curr = nodes[kmer];

		if (curr == NULL) {
			curr = new_node(kmer, sequence);

			if (curr == NULL) {
				printf("Null node for kmer: %s\n", kmer);
				exit(-1);
			}

			nodes[kmer] = curr;
		} else {
			increment_node_freq(curr, sequence);
		}

		if (prev != NULL) {
			link_nodes(prev, curr);
		}

		prev = curr;
	}
}

void build_graph(const char* read_file) {
	FILE *fp = fopen(read_file, "r");
	char read[101];
	memset(read, 0, 101);

	int line = 0;
	while (fscanf(fp, "%s", read) != EOF) {
//		printf("read: %d : %s\n", line++, read);
		add_to_graph(read);
		line++;
	}

	printf("Num reads: %d\n", line);

	fclose(fp);
}

struct linked_node* remove_node_from_list(struct node* node, struct linked_node* list) {
	struct linked_node* node_ptr = list;
	struct linked_node* prev_ptr = NULL;

	char is_found = false;
	while ((node_ptr != NULL) && (!is_found)) {
		if (strcmp(node_ptr->node->seq, node->seq) == 0) {
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

void prune_graph() {
	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes.begin();
		         it != nodes.end(); ++it) {

		const char* key = it->first;
		struct node* node = it->second;

		if ((node != NULL) && ((node->frequency < MIN_NODE_FREQUENCY) || (!(node->hasMultipleUniqueReads)))) {

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
			nodes.erase(key);

			// Free memory
			free(node);
		}
	}
}

struct linked_node* identify_root_nodes() {

	struct linked_node* root_nodes = NULL;
	int count = 0;

	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes.begin();
	         it != nodes.end(); ++it) {
		struct node* node = it->second;
		if ((node != NULL) && (node->fromNodes == NULL)) {
			struct linked_node* next = root_nodes;
			root_nodes = (linked_node*) malloc(sizeof(linked_node));
			root_nodes->node = node;
			root_nodes->next = next;

			printf("Root: %s\n", node->seq);
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
	memset(curr_contig->seq, 0, sizeof(curr_contig->seq));
	curr_contig->size = 0;
	curr_contig->is_repeat = 0;
	curr_contig->visited_nodes = new sparse_hash_set<const char*, my_hash, eqstr>();
	return curr_contig;
}

struct contig* copy_contig(struct contig* orig) {
	struct contig* copy = (contig*) malloc(sizeof(contig));
	memcpy(copy->seq, orig->seq, MAX_CONTIG_SIZE);
	copy->size = orig->size;
	copy->is_repeat = orig->is_repeat;
	copy->visited_nodes = new sparse_hash_set<const char*, my_hash, eqstr>(*orig->visited_nodes);
}

void free_contig(struct contig* contig) {
	delete contig->visited_nodes;
	memset(contig, 0, sizeof(contig));
	free(contig);
}

char is_node_visited(struct contig* contig, struct node* node) {
	 sparse_hash_set<const char*, my_hash, eqstr>::const_iterator it = contig->visited_nodes->find(node->seq);
	 return it != contig->visited_nodes->end();
}

void output_contig(struct contig* contig, int& contig_count, FILE* fp, const char* prefix) {
	if (strlen(contig->seq) >= MIN_CONTIG_LENGTH) {
		if (contig->is_repeat) {
			fprintf(fp, ">%s_%d_repeat\n%s\n", prefix, contig_count++, contig->seq);
		} else {
			fprintf(fp, ">%s_%d\n%s\n", prefix, contig_count++, contig->seq);
		}
	}
}

void build_contigs(struct node* root, int& contig_count, FILE* fp, const char* prefix) {

	stack<contig*> contigs;
	struct contig* root_contig = new_contig();
//	contig->seq[0] = root->seq[0];
//	contig->size = 1;
	root_contig->curr_node = root;
	contigs.push(root_contig);

	while (contigs.size() > 0) {
		// Get contig from stack
		struct contig* contig = contigs.top();

		if (is_node_visited(contig, contig->curr_node)) {
			// We've encountered a repeat
			contig->is_repeat = 1;
			output_contig(contig, contig_count, fp, prefix);
			contigs.pop();
			free_contig(contig);
		}
		else if (contig->curr_node->toNodes == NULL) {
			// We've reached the end of the contig.
			// Append entire current node.
			memcpy(&(contig->seq[contig->size]), contig->curr_node->seq, KMER);

			// Now, write the contig
			output_contig(contig, contig_count, fp, prefix);
			contigs.pop();
			free_contig(contig);
		}
		else {
			// Append first base from current node
			contig->seq[contig->size++] = contig->curr_node->seq[0];
			if (contig->size >= MAX_CONTIG_SIZE) {
				printf("Max contig size exceeded at node: %s\n", contig->curr_node->seq);
			}

			contig->visited_nodes->insert(contig->curr_node->seq);

			// Move current contig to next "to" node.
			struct linked_node* to_linked_node = contig->curr_node->toNodes;
			contig->curr_node = to_linked_node->node;
			to_linked_node = to_linked_node->next;

			// If there are multiple "to" nodes, branch the contig and push on stack
			while (to_linked_node != NULL) {
				struct contig* contig_branch = copy_contig(contig);
//				memcpy(contig_branch, contig, sizeof(contig));
				contig_branch->curr_node = to_linked_node->node;
				contigs.push(contig_branch);
				to_linked_node = to_linked_node->next;
			}
		}
	}
}

void cleanup() {
	for (sparse_hash_map<const char*, struct node*, my_hash, eqstr>::const_iterator it = nodes.begin();
	         it != nodes.end(); ++it) {
		char* key = (char*) it->first;
		struct node* node = it->second;

		if (node != NULL) {
			free(node);
		}

		if (key != NULL) {
			free(key);
		}
	}
}

void assemble(const char* input,
			  const char* output,
			  const char* prefix,
			  char truncate_on_repeat,
			  int max_contigs,
			  int max_paths_from_root) {

	printf("Assembling: %s -> %s", input, output);
	nodes.set_deleted_key(NULL);

	build_graph(input);
	prune_graph();

	struct linked_node* root_nodes = identify_root_nodes();

	int contig_count = 0;
	FILE *fp = fopen(output, "w");
	while (root_nodes != NULL) {
		build_contigs(root_nodes->node, contig_count, fp, prefix);
		root_nodes = root_nodes->next;
	}
	fclose(fp);

	cleanup();

	printf("Done assembling: %s -> %s", input, output);
}

int main(int argc, char* argv[]) {


	assemble(
		"/home/lmose/code/abra/src/main/c/sim83_reads_filtered.txt",
		"/home/lmose/code/abra/src/main/c/sim83_c.fasta",
		"foo",
		true,
		50000,
		5000);
}
