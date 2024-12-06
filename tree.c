#include <stdlib.h>
#include "tree.h"
#include "distsNmore.h"

typedef struct tnode *TreeNodePtr;

typedef struct tnode {
  int index;
  TreeNodePtr left;
  TreeNodePtr right;
} TreeNode;

typedef struct TreeContainer {
	TreeNodePtr root;
	int size;
} TreeContainer;

void initTree(TreePtr* tree) {
	*tree = safeMalloc(sizeof(TreeContainer));
	(*tree)->root = NULL;
	(*tree)->size = 0;
}

TreeNodePtr addNode(TreePtr tree, TreeNodePtr p, int index) {
	if (p == NULL) {
		p = malloc(sizeof(TreeNode));
		p->index = index;
		p->left = NULL;
		p->right = NULL;
		if (tree != NULL) tree->size++;
	} else if (index < p->index) {
		p->left = addNode(tree, p->left, index);
	} else if (index > p->index) {
		p->right = addNode(tree, p->right, index);
	}
	return p;
}

void treeInsert(TreePtr tree, int item) {
	tree->root = addNode(tree, tree->root, item);
}

int getTreeSize(TreePtr tree) {
	return tree->size;
}

void destroyNodes(TreeNodePtr* p) {
	if (p == NULL) return;
	if ((*p)->left != NULL) destroyNodes(&((*p)->left));
	if ((*p)->right != NULL) destroyNodes(&((*p)->right));
	free(*p);
}

void destroyTree(TreePtr* p) {
	destroyNodes(&((*p)->root));
	free(*p);
}

void printNodes(TreePtr tree, TreeNodePtr p, FILE* output) {
	if (p == NULL) return;
	printNodes(tree, p->left, output);
	fprintf(output, "%d\t", p->index);
	printNodes(tree, p->right, output);
}

void treePrint(TreePtr tree, FILE* output) {
	printNodes(tree, tree->root, output);
}
