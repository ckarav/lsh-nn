#ifndef _RANDZ_
#define _RANDZ_


/* returns a uniform random int in range [M,N] */
int uniform_randInt(int M, int N);

double uniform_randR(int N); //uniform random double in [0,N)

/*stores d random values(that follow the normal distribution ~N(0,1)) in the "v" vector*/
void gauss_rand(double *v, int d);
void gauss_rand2(double *v, int d);//same as above but it may have negative coordinates

double uniformRandDouble(double N); //uniform random double in [0,N)

#endif // _RANDZ_
