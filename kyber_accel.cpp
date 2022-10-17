#include "kyber_accel.hpp"

void print_polyvec_hls(polyvec *p)
{
  unsigned int i,j;

  printf("\n");
  for(i=0 ; i < KYBER_K ; i++)
  {
	  printf("[%d]=",i);
	  for(j = 0 ; j < KYBER_N ; j++)
	  {
		  printf(" %d", p->vec[i].coeffs[j]);
	  }
	  printf("\n");
  }
}

void print_poly_hls(poly *p)
{
  unsigned int j;

  printf("\n");

  for(j = 0 ; j < KYBER_N ; j++)
  {
	  printf(" %d", p->coeffs[j]);
  }
  printf("\n");

}

int16_t montgomery_reduce_hls(int32_t a)
{
  int16_t t;

  t = (int16_t)a*QINV_hls;
  t = (a - (int32_t)t*KYBER_Q_hls) >> 16;
  return t;
}

int16_t fqmul_hls(int16_t a, int16_t b) {
  return montgomery_reduce_hls((int32_t)a*b);
}

void ntt_hls(int16_t r[256]) {
  unsigned int len, start, j, k;
  int16_t t, zeta;

  k = 1;
  for(len = 128; len >= 2; len >>= 1) {
//#pragma HLS PIPELINE
    for(start = 0; start < 256; start = j + len) {
//#pragma HLS PIPELINE
      zeta = zetas_hls[k++];
      for(j = start; j < start + len; j++) {
//#pragma HLS PIPELINE
        t = fqmul_hls(zeta, r[j + len]);
        r[j + len] = r[j] - t;
        r[j] = r[j] + t;
      }
    }
  }
}

int16_t barrett_reduce_hls(int16_t a) {
  int16_t t;
  const int16_t v = ((1<<26) + KYBER_Q_hls/2)/KYBER_Q_hls;

  t  = ((int32_t)v*a + (1<<25)) >> 26;
  t *= KYBER_Q_hls;

  return a - t;
}


void poly_reduce_hls(poly *r)
{
  unsigned int i;
  for(i=0;i<KYBER_N_hls;i++)
//#pragma HLS PIPELINE
    r->coeffs[i] = barrett_reduce_hls(r->coeffs[i]);
}


void poly_ntt_hls(poly *r)
{
  ntt_hls(r->coeffs);
  poly_reduce_hls(r);
}

void polyvec_ntt_hls(polyvec *r)
{

  unsigned int i;
  for(i=0;i<KYBER_K_hls;i++)
//#pragma HLS PIPELINE
    poly_ntt_hls(&r->vec[i]);
}

void poly_add_hls(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<KYBER_N;i++)
    r->coeffs[i] = a->coeffs[i] + b->coeffs[i];
}

void basemul_hls(int16_t r[2], const int16_t a[2], const int16_t b[2], int16_t zeta)
{
  r[0]  = fqmul_hls(a[1], b[1]);
  r[0]  = fqmul_hls(r[0], zeta);
  r[0] += fqmul_hls(a[0], b[0]);
  r[1]  = fqmul_hls(a[0], b[1]);
  r[1] += fqmul_hls(a[1], b[0]);
}


void poly_basemul_montgomery_hls(poly *r, const poly *a, const poly *b)
{
  unsigned int i;
  for(i=0;i<KYBER_N/4;i++) {
    basemul_hls(&r->coeffs[4*i], &a->coeffs[4*i], &b->coeffs[4*i], zetas_hls[64+i]);
    basemul_hls(&r->coeffs[4*i+2], &a->coeffs[4*i+2], &b->coeffs[4*i+2], -zetas_hls[64+i]);
  }
}


void kyber_accelerator(volatile poly *r, volatile const polyvec *a, volatile const polyvec *b)
{
	#pragma HLS INTERFACE mode=m_axi depth=1 port=r offset=slave
	#pragma HLS INTERFACE mode=m_axi depth=1 port=a offset=slave
	#pragma HLS INTERFACE mode=m_axi depth=1 port=b offset=slave
	#pragma HLS INTERFACE mode=s_axilite port=r bundle=CONTROL
	#pragma HLS INTERFACE mode=s_axilite port=a bundle=CONTROL
	#pragma HLS INTERFACE mode=s_axilite port=b bundle=CONTROL
	#pragma HLS INTERFACE mode=s_axilite port=return bundle=CONTROL

	unsigned int i,j;
	poly t;
	poly r_hls;
	polyvec a_hls;
	polyvec b_hls;

	//copy data into the accelerator
	for (j = 0 ; j < KYBER_N_hls ; j++)
	{
		#pragma HLS PIPELINE
		r_hls.coeffs[j] = r->coeffs[j];
	}

	//copy data into the accelerator
	for (i = 0 ; i < KYBER_K_hls ; i++)
	{
		for (j = 0 ; j < KYBER_N_hls ; j++)
		{
			#pragma HLS PIPELINE
			a_hls.vec[i].coeffs[j] = a->vec[i].coeffs[j];
		}
	}

	//copy data into the accelerator
	for (i = 0 ; i < KYBER_K_hls ; i++)
	{
		for (j = 0 ; j < KYBER_N_hls ; j++)
		{
			#pragma HLS PIPELINE
			b_hls.vec[i].coeffs[j] = b->vec[i].coeffs[j];
		}
	}

	/*********************************************************/
	/************** polynomials multiplication ***************/
	/*********************************************************/
	poly_basemul_montgomery_hls(&r_hls, &a_hls.vec[0], &b_hls.vec[0]);

	for(i=1;i<KYBER_K_hls;i++) {
		poly_basemul_montgomery_hls(&t, &a_hls.vec[i], &b_hls.vec[i]);
		poly_add_hls(&r_hls, &r_hls, &t);
	}

	poly_reduce_hls(&r_hls);
	/*********************************************************/
	/*********************************************************/
	/*********************************************************/

	//copy computed data back to main processor
	for (j = 0 ; j < KYBER_N_hls ; j++)
	{
		#pragma HLS PIPELINE
		r->coeffs[j] = r_hls.coeffs[j];
	}

}
