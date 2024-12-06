all: main.c multi.h enums.h parsers.c parsers.h ClusterElementsNum.c ClusterElementsNum.h cosine.c cosine.h hamming.c hamming.h distsNmore.c distsNmore.h rands.c rands.h vector.c vector.h matrix.c matrix.h assignment.c assignment.h evaluation.c evaluation.h init.c init.h update.c update.h recommend.c recommend.h testing.h tree.c tree.h
	gcc -O3 -o recommendation main.c parsers.c ClusterElementsNum.c cosine.c hamming.c distsNmore.c rands.c vector.c matrix.c assignment.c evaluation.c update.c init.c recommend.c tree.c -ldl -lm -lpthread -lgsl -lgslcblas -lcunit
clean:
	rm -rf recommendation
