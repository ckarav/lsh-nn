#ifndef __TREE__
#define __TREE__

#include <stdio.h>

typedef struct TreeContainer* TreePtr;

void initTree(TreePtr* tree);
void treeInsert(TreePtr tree, int item);
int getTreeSize(TreePtr tree);
void destroyTree(TreePtr* p);
void treePrint(TreePtr tree, FILE* output);

#endif
