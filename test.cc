#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include <stdlib.h>
#include <stdio.h>

#include "blas.h"

#include "rand.h"
#include "mat.h"
#include "nmf.h"

#include <fstream>
#include <string>
#include <sstream>

int main(int argc, char **argv)
{

  int M;
  int N;
  int K = atoi(argv[2]);
  int max_iter = atoi(argv[3]);

  float *_W = {0, };// = (float *)malloc(sizeof(float)*M * 10);
  float *_H = {0, };// = (float *)malloc(sizeof(float) * 10 * N);
  float *_V = {0, };// = (float *)malloc(sizeof(float)*M*N);
  
  std::string line, entry;
  std::ifstream f(argv[1]);
  std::stringstream ss;

  int M_prev = 0;

  int idx;

  M = 0;
  N = 0;
  // Phase 1 : Get total count
  bool first = true;
  while (std::getline(f, line)) {
	  ss.clear();
	  ss << line;
	  M_prev = 0;
	  while (std::getline(ss, entry, '\t')) {
		  M_prev++;
	  }
	  if (!first && M != M_prev) {
		  printf("Incorrect size!\n");
		  exit(1);
	  }
	  else {
		  first = false;
		  M = M_prev;
	  }
	  N++;
  }

  // Phase 2 : Allocate spaces
  _W = (float *)malloc(sizeof(float)*M*K);
  _H = (float *)malloc(sizeof(float)*K*N);
  _V = (float *)malloc(sizeof(float)*M*N);

  // Phase 3 : Read actual values
  f.clear(); f.seekg(0, f.beg);
  idx = 0;
  while (std::getline(f, line)) {
	  ss.clear();
	  ss << line;
	  while (std::getline(ss, entry, '\t')) {
		  _V[idx++] = atof(entry.c_str());
	  }
  }
  f.close();

  init_rng();

  for (int j=0; j<M; j++)
    for (int i=0; i<K; i++)
      _W[j*K+i]=uniff();
  
  for (int j=0; j<K; j++)
    for (int i=0; i<N; i++)
      _H[j*N+i]=uniff();

  mmultf(_W,_H,_V,M,N,K);

  for (int i=0; i<M*N; i++)
    _V[i]+=uniff()*0.0;//5;

  float *WH = (float *)malloc(sizeof(float)*M*N);

  std::vector<double> ssds;

 // FILE *fp = fopen("nmf_ssd.txt","wt");
    printf("factorig for k=%d:\n",K); fflush(stdout);
    nmfssdf(_V,
	    NULL,
	    M,
	    N,
	    _W,
	    _H,
	    K,
	    WH,
	    1000,
	    max_iter,
	    true,
	    true,
	    NULL,
	    NULL,
	    &ssds);
    //fprintf(fp,"%d %0.5f\n",K,ssds.back());
    ssds.clear();
  //fclose(fp);

  FILE *fp = fopen("nmf_w.txt", "wt");
  for (int j = 0; j < M; j++) {
	  first = true;
	  for (int i = 0; i < K; i++) {
		  if (!first) fputc('\t', fp);
		  fprintf(fp, "%f", _W[j*K + i]);
		  first = false;
	  }
	  fputc('\n', fp);
  }
  fclose(fp);

  fp = fopen("nmf_h.txt", "wt");
  for (int j = 0; j < K; j++) {
	  first = true;
	  for (int i = 0; i < N; i++) {
		  if (!first) fputc('\t', fp);
		  fprintf(fp, "%f", _H[j*K + i]);
		  first = false;
	  }
	  fputc('\n', fp);
  }
  fclose(fp);

  free((void *)_V);
  free((void *)_H);
  free((void *)_W);
  
  return 0;
}
