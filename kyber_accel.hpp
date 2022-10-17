//#include "ap_int.h"

//typedef ap_int<16> int16_t;
//typedef ap_int<32> integer32;

#include <stdio.h>
#include "polyvec.hpp"
#include "poly.hpp"


#define QINV_hls -3327 // q^-1 mod 2^16
#define KYBER_Q_hls 3329
#define KYBER_N_hls 256
#define KYBER_K_hls 3	// Change this for different security strengths

/*
typedef struct{
  integer16 coeffs[KYBER_N_hls];
} poly_hls;

typedef struct{
  poly_hls vec[KYBER_K_hls];
} polyvec_hls;
*/

//static integer16 polyvec_array[KYBER_K*KYBER_N];


static const int16_t zetas_hls[128] = {
  -1044,  -758,  -359, -1517,  1493,  1422,   287,   202,
   -171,   622,  1577,   182,   962, -1202, -1474,  1468,
    573, -1325,   264,   383,  -829,  1458, -1602,  -130,
   -681,  1017,   732,   608, -1542,   411,  -205, -1571,
   1223,   652,  -552,  1015, -1293,  1491,  -282, -1544,
    516,    -8,  -320,  -666, -1618, -1162,   126,  1469,
   -853,   -90,  -271,   830,   107, -1421,  -247,  -951,
   -398,   961, -1508,  -725,   448, -1065,   677, -1275,
  -1103,   430,   555,   843, -1251,   871,  1550,   105,
    422,   587,   177,  -235,  -291,  -460,  1574,  1653,
   -246,   778,  1159,  -147,  -777,  1483,  -602,  1119,
  -1590,   644,  -872,   349,   418,   329,  -156,   -75,
    817,  1097,   603,   610,  1322, -1285, -1465,   384,
  -1215,  -136,  1218, -1335,  -874,   220, -1187, -1659,
  -1185, -1530, -1278,   794, -1510,  -854,  -870,   478,
   -108,  -308,   996,   991,   958, -1460,  1522,  1628
};

void kyber_accelerator(volatile poly *r, volatile const polyvec *a, volatile const polyvec *b);

//void kyber_accelerator(volatile polyvec *pv);

