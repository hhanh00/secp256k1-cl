#define WINDOW_A 5
#define WINDOW_G 14

/** The number of entries a table with precomputed multiples needs to have. */
#define ECMULT_TABLE_SIZE(w) (1 << ((w)-2))

#define ECMULT_TABLE_GET(r,pre,n,w,neg) do { \
    if ((n) > 0) \
        *(r) = (pre)[((n)-1)/2]; \
    else \
        (neg)((r), &(pre)[(-(n)-1)/2]); \
} while(0)


#define ECMULT_TABLE_GET_GEJ(r,pre,n,w) ECMULT_TABLE_GET((r),(pre),(n),(w),secp256k1_gej_neg)
#define ECMULT_TABLE_GET_GE(r,pre,n,w)  ECMULT_TABLE_GET((r),(pre),(n),(w),secp256k1_ge_neg)

typedef struct {
    // X = sum(i=0..9, elem[i]*2^26) mod n
    unsigned int n[10];
} secp256k1_fe_t;

typedef struct {
    secp256k1_fe_t x;
    secp256k1_fe_t y;
    int infinity; // whether this represents the point at infinity
} secp256k1_ge_t;

    typedef struct {
    secp256k1_fe_t x; // actual X: x/z^2
    secp256k1_fe_t y; // actual Y: y/z^3
    secp256k1_fe_t z;
    int infinity; // whether this represents the point at infinity
} secp256k1_gej_t;

typedef struct {
	secp256k1_gej_t a;
	char wnaf_na[257];
	short bits_na;
	short wnaf_ng_1[129];
	short bits_ng_1;
	short wnaf_ng_128[129];
	short bits_ng_128;
} ecmult_params_device_t;

typedef struct {
    // For accelerating the computation of a*P + b*G:
    secp256k1_ge_t pre_g[ECMULT_TABLE_SIZE(WINDOW_G)];    // odd multiples of the generator
    secp256k1_ge_t pre_g_128[ECMULT_TABLE_SIZE(WINDOW_G)]; // odd multiples of 2^128*generator

    // For accelerating the computation of a*G:
    // To harden against timing attacks, use the following mechanism:
    // * Break up the multiplicand into groups of 4 bits, called n_0, n_1, n_2, ..., n_63.
    // * Compute sum((n_i + 1) * 16^i * G, i=0..63).
    // * Subtract sum(1 * 16^i * G, i=0..63).
    // For each i, and each of the 16 possible values of n_i, ((n_i + 1) * 16^i * G) is
    // precomputed (call it prec(i,n_i), as well as -sum(1 * 16^i * G) (called fin).
    // The formula now becomes sum(prec(i, n_i), i=0..63) + fin.
    // To make memory access uniform, the bytes of prec(i,n_i) are sliced per value of n_i.
    unsigned char prec[64][sizeof(secp256k1_ge_t)][16]; // prec[j][k][i] = k'th byte of (16^j * (i+1) * G)
    secp256k1_ge_t fin; // -(sum(prec[j][0], j=0..63))
} secp256k1_ecmult_consts_t;

    
void secp256k1_gej_set_infinity(secp256k1_gej_t *r) {
    r->infinity = 1;
}
void secp256k1_gej_neg(secp256k1_gej_t *r, __global const secp256k1_gej_t *a);
void secp256k1_ge_neg(secp256k1_ge_t *r, __constant const secp256k1_ge_t *a);


void secp256k1_ecmult_table_precomp_gej(__global secp256k1_gej_t *pre, __global secp256k1_gej_t *a, int w);
void secp256k1_gej_double(secp256k1_gej_t *r, secp256k1_gej_t *a);
void secp256k1_gej_add(__global secp256k1_gej_t *r, const secp256k1_gej_t *a, __global const secp256k1_gej_t *b);
void secp256k1_gej_add_local(secp256k1_gej_t *pr, const secp256k1_gej_t *pa, const secp256k1_gej_t *pb);
void secp256k1_gej_add_ge(secp256k1_gej_t *r, const secp256k1_gej_t *a, const secp256k1_ge_t *b);

__kernel void secp256k1_ecmult_table_precomp_gej(__global secp256k1_gej_t *pre, __global secp256k1_gej_t *pa, int w) {
    int ind = get_global_id(0);
    int sz = (1 << (w-2));
    secp256k1_gej_t a = pa[ind];
    pre[ind * sz] = a;
    secp256k1_gej_t d; secp256k1_gej_double(&d, &a);
    for(int i = 1; i < sz; i++) {
        secp256k1_gej_add(&pre[ind * sz + i], &d, &pre[ind * sz + i - 1]);
    }
}

__kernel void test_ecmult_table_precomp_gej(__global secp256k1_gej_t *pre, __global secp256k1_gej_t *pa, int w) {
    int ind = get_global_id(0);
    pre[ind] = pa[ind];
    // pre[ind].infinity = 2;
}

__kernel void secp256k1_ecmult(__global secp256k1_gej_t *pr, __global secp256k1_gej_t *ppre_a, __global ecmult_params_device_t *pParams, __constant const secp256k1_ecmult_consts_t *mc) {
    int ind = get_global_id(0);
    int sz = (1 << (WINDOW_A-2));

    ecmult_params_device_t params = pParams[ind];
    __global secp256k1_gej_t *pre_a = &ppre_a[ind * sz];
    
    secp256k1_gej_t r;
    secp256k1_gej_set_infinity(&r);
    secp256k1_gej_t tmpj;
    secp256k1_ge_t tmpa;

    short bits = max(params.bits_na, params.bits_ng_1);
    bits = max(bits, params.bits_ng_128);
    for (int i=bits-1; i>=0; i--) {
        secp256k1_gej_double(&r, &r);
        int n;
        if (i < params.bits_na && (n = params.wnaf_na[i])) {
            ECMULT_TABLE_GET_GEJ(&tmpj, pre_a, n, WINDOW_A);
            secp256k1_gej_add_local(&r, &r, &tmpj);
        }

        if (i < params.bits_ng_1 && (n = params.wnaf_ng_1[i])) {
            ECMULT_TABLE_GET_GE(&tmpa, mc->pre_g, n, WINDOW_G);
            secp256k1_gej_add_ge(&r, &r, &tmpa);
        }
        if (i < params.bits_ng_128 && (n = params.wnaf_ng_128[i])) {
            ECMULT_TABLE_GET_GE(&tmpa, mc->pre_g_128, n, WINDOW_G);
            secp256k1_gej_add_ge(&r, &r, &tmpa);
        }
    }
    
    pr[ind] = r;
}
                      
void secp256k1_fe_normalize(secp256k1_fe_t *r);    
void secp256k1_fe_mul(secp256k1_fe_t *r, const secp256k1_fe_t *a, const secp256k1_fe_t *b);
void secp256k1_fe_mul_int(secp256k1_fe_t *r, int a);
void secp256k1_fe_add(secp256k1_fe_t *r, const secp256k1_fe_t *a);
void secp256k1_fe_sqr(secp256k1_fe_t *r, const secp256k1_fe_t *a);
void secp256k1_fe_negate(secp256k1_fe_t *r, const secp256k1_fe_t *a, int m);
int secp256k1_fe_is_zero(const secp256k1_fe_t *a);

void secp256k1_gej_neg(secp256k1_gej_t *r, __global const secp256k1_gej_t *a) {
    r->infinity = a->infinity;
    r->x = a->x;
    r->y = a->y;
    r->z = a->z;
    secp256k1_fe_normalize(&r->y);
    secp256k1_fe_negate(&r->y, &r->y, 1);
}

void secp256k1_ge_neg(secp256k1_ge_t *r, __constant const secp256k1_ge_t *a) {
    r->infinity = a->infinity;
    r->x = a->x;
    r->y = a->y;
    secp256k1_fe_normalize(&r->y);
    secp256k1_fe_negate(&r->y, &r->y, 1);
}

typedef unsigned int uint32_t;
typedef unsigned long uint64_t;

void secp256k1_fe_normalize(secp256k1_fe_t *r) {
//    fog("normalize in: ", r);
    uint32_t c;
    c = r->n[0];
    uint32_t t0 = c & 0x3FFFFFFUL;
    c = (c >> 26) + r->n[1];
    uint32_t t1 = c & 0x3FFFFFFUL;
    c = (c >> 26) + r->n[2];
    uint32_t t2 = c & 0x3FFFFFFUL;
    c = (c >> 26) + r->n[3];
    uint32_t t3 = c & 0x3FFFFFFUL;
    c = (c >> 26) + r->n[4];
    uint32_t t4 = c & 0x3FFFFFFUL;
    c = (c >> 26) + r->n[5];
    uint32_t t5 = c & 0x3FFFFFFUL;
    c = (c >> 26) + r->n[6];
    uint32_t t6 = c & 0x3FFFFFFUL;
    c = (c >> 26) + r->n[7];
    uint32_t t7 = c & 0x3FFFFFFUL;
    c = (c >> 26) + r->n[8];
    uint32_t t8 = c & 0x3FFFFFFUL;
    c = (c >> 26) + r->n[9];
    uint32_t t9 = c & 0x03FFFFFUL;
    c >>= 22;
/*    r->n[0] = t0; r->n[1] = t1; r->n[2] = t2; r->n[3] = t3; r->n[4] = t4;
    r->n[5] = t5; r->n[6] = t6; r->n[7] = t7; r->n[8] = t8; r->n[9] = t9;
    fog("         tm1: ", r);
    fprintf(stderr, "out c= %08lx\n", (unsigned long)c);*/

    // The following code will not modify the t's if c is initially 0.
    uint32_t d = c * 0x3D1UL + t0;
    t0 = d & 0x3FFFFFFUL;
    d = (d >> 26) + t1 + c*0x40;
    t1 = d & 0x3FFFFFFUL;
    d = (d >> 26) + t2;
    t2 = d & 0x3FFFFFFUL;
    d = (d >> 26) + t3;
    t3 = d & 0x3FFFFFFUL;
    d = (d >> 26) + t4;
    t4 = d & 0x3FFFFFFUL;
    d = (d >> 26) + t5;
    t5 = d & 0x3FFFFFFUL;
    d = (d >> 26) + t6;
    t6 = d & 0x3FFFFFFUL;
    d = (d >> 26) + t7;
    t7 = d & 0x3FFFFFFUL;
    d = (d >> 26) + t8;
    t8 = d & 0x3FFFFFFUL;
    d = (d >> 26) + t9;
    t9 = d & 0x03FFFFFUL;
/*    r->n[0] = t0; r->n[1] = t1; r->n[2] = t2; r->n[3] = t3; r->n[4] = t4;
    r->n[5] = t5; r->n[6] = t6; r->n[7] = t7; r->n[8] = t8; r->n[9] = t9;
    fog("         tm2: ", r); */

    // Subtract p if result >= p
    uint64_t low = ((uint64_t)t1 << 26) | t0;
    uint64_t mask = -(long)((t9 < 0x03FFFFFUL) | (t8 < 0x3FFFFFFUL) | (t7 < 0x3FFFFFFUL) | (t6 < 0x3FFFFFFUL) | (t5 < 0x3FFFFFFUL) | (t4 < 0x3FFFFFFUL) | (t3 < 0x3FFFFFFUL) | (t2 < 0x3FFFFFFUL) | (low < 0xFFFFEFFFFFC2FUL));
    t9 &= mask;
    t8 &= mask;
    t7 &= mask;
    t6 &= mask;
    t5 &= mask;
    t4 &= mask;
    t3 &= mask;
    t2 &= mask;
    low -= (~mask & 0xFFFFEFFFFFC2FUL);

    // push internal variables back
    r->n[0] = low & 0x3FFFFFFUL; r->n[1] = (low >> 26) & 0x3FFFFFFUL; r->n[2] = t2; r->n[3] = t3; r->n[4] = t4;
    r->n[5] = t5; r->n[6] = t6; r->n[7] = t7; r->n[8] = t8; r->n[9] = t9;
/*    fog("         out: ", r);*/

#ifdef VERIFY
    r->magnitude = 1;
    r->normalized = 1;
#endif
}

void static inline secp256k1_fe_set_int(secp256k1_fe_t *r, int a) {
    r->n[0] = a;
    r->n[1] = r->n[2] = r->n[3] = r->n[4] = r->n[5] = r->n[6] = r->n[7] = r->n[8] = r->n[9] = 0;
#ifdef VERIFY
    r->magnitude = 1;
    r->normalized = 1;
#endif
}

// TODO: not constant time!
int secp256k1_fe_is_zero(const secp256k1_fe_t *a) {
#ifdef VERIFY
    assert(a->normalized);
#endif
    return (a->n[0] == 0 && a->n[1] == 0 && a->n[2] == 0 && a->n[3] == 0 && a->n[4] == 0 && a->n[5] == 0 && a->n[6] == 0 && a->n[7] == 0 && a->n[8] == 0 && a->n[9] == 0);
}

int static inline secp256k1_fe_is_odd(const secp256k1_fe_t *a) {
#ifdef VERIFY
    assert(a->normalized);
#endif
    return a->n[0] & 1;
}

// TODO: not constant time!
int secp256k1_fe_equal(const secp256k1_fe_t *a, const secp256k1_fe_t *b) {
#ifdef VERIFY
    assert(a->normalized);
    assert(b->normalized);
#endif
    return (a->n[0] == b->n[0] && a->n[1] == b->n[1] && a->n[2] == b->n[2] && a->n[3] == b->n[3] && a->n[4] == b->n[4] &&
            a->n[5] == b->n[5] && a->n[6] == b->n[6] && a->n[7] == b->n[7] && a->n[8] == b->n[8] && a->n[9] == b->n[9]);
}

void secp256k1_fe_set_b32(secp256k1_fe_t *r, const unsigned char *a) {
    r->n[0] = r->n[1] = r->n[2] = r->n[3] = r->n[4] = 0;
    r->n[5] = r->n[6] = r->n[7] = r->n[8] = r->n[9] = 0;
    for (int i=0; i<32; i++) {
        for (int j=0; j<4; j++) {
            int limb = (8*i+2*j)/26;
            int shift = (8*i+2*j)%26;
            r->n[limb] |= (uint32_t)((a[31-i] >> (2*j)) & 0x3) << shift;
        }
    }
#ifdef VERIFY
    r->magnitude = 1;
    r->normalized = 1;
#endif
}

/** Convert a field element to a 32-byte big endian value. Requires the input to be normalized */
void static secp256k1_fe_get_b32(unsigned char *r, const secp256k1_fe_t *a) {
#ifdef VERIFY
    assert(a->normalized);
#endif
    for (int i=0; i<32; i++) {
        int c = 0;
        for (int j=0; j<4; j++) {
            int limb = (8*i+2*j)/26;
            int shift = (8*i+2*j)%26;
            c |= ((a->n[limb] >> shift) & 0x3) << (2 * j);
        }
        r[31-i] = c;
    }
}

void secp256k1_fe_negate(secp256k1_fe_t *r, const secp256k1_fe_t *a, int m) {
#ifdef VERIFY
    assert(a->magnitude <= m);
    r->magnitude = m + 1;
    r->normalized = 0;
#endif
    r->n[0] = 0x3FFFC2FUL * (m + 1) - a->n[0];
    r->n[1] = 0x3FFFFBFUL * (m + 1) - a->n[1];
    r->n[2] = 0x3FFFFFFUL * (m + 1) - a->n[2];
    r->n[3] = 0x3FFFFFFUL * (m + 1) - a->n[3];
    r->n[4] = 0x3FFFFFFUL * (m + 1) - a->n[4];
    r->n[5] = 0x3FFFFFFUL * (m + 1) - a->n[5];
    r->n[6] = 0x3FFFFFFUL * (m + 1) - a->n[6];
    r->n[7] = 0x3FFFFFFUL * (m + 1) - a->n[7];
    r->n[8] = 0x3FFFFFFUL * (m + 1) - a->n[8];
    r->n[9] = 0x03FFFFFUL * (m + 1) - a->n[9];
}

void secp256k1_fe_mul_int(secp256k1_fe_t *r, int a) {
#ifdef VERIFY
    r->magnitude *= a;
    r->normalized = 0;
#endif
    r->n[0] *= a;
    r->n[1] *= a;
    r->n[2] *= a;
    r->n[3] *= a;
    r->n[4] *= a;
    r->n[5] *= a;
    r->n[6] *= a;
    r->n[7] *= a;
    r->n[8] *= a;
    r->n[9] *= a;
}

void secp256k1_fe_add(secp256k1_fe_t *r, const secp256k1_fe_t *a) {
#ifdef VERIFY
    r->magnitude += a->magnitude;
    r->normalized = 0;
#endif
    r->n[0] += a->n[0];
    r->n[1] += a->n[1];
    r->n[2] += a->n[2];
    r->n[3] += a->n[3];
    r->n[4] += a->n[4];
    r->n[5] += a->n[5];
    r->n[6] += a->n[6];
    r->n[7] += a->n[7];
    r->n[8] += a->n[8];
    r->n[9] += a->n[9];
}

void secp256k1_fe_mul_inner(const uint32_t *a, const uint32_t *b, uint32_t *r) {
    uint64_t c = (uint64_t)a[0] * b[0];
    uint32_t t0 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[0] * b[1] +
            (uint64_t)a[1] * b[0];
    uint32_t t1 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[0] * b[2] +
            (uint64_t)a[1] * b[1] +
            (uint64_t)a[2] * b[0];
    uint32_t t2 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[0] * b[3] +
            (uint64_t)a[1] * b[2] +
            (uint64_t)a[2] * b[1] +
            (uint64_t)a[3] * b[0];
    uint32_t t3 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[0] * b[4] +
            (uint64_t)a[1] * b[3] +
            (uint64_t)a[2] * b[2] +
            (uint64_t)a[3] * b[1] +
            (uint64_t)a[4] * b[0];
    uint32_t t4 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[0] * b[5] +
            (uint64_t)a[1] * b[4] +
            (uint64_t)a[2] * b[3] +
            (uint64_t)a[3] * b[2] +
            (uint64_t)a[4] * b[1] +
            (uint64_t)a[5] * b[0];
    uint32_t t5 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[0] * b[6] +
            (uint64_t)a[1] * b[5] +
            (uint64_t)a[2] * b[4] +
            (uint64_t)a[3] * b[3] +
            (uint64_t)a[4] * b[2] +
            (uint64_t)a[5] * b[1] +
            (uint64_t)a[6] * b[0];
    uint32_t t6 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[0] * b[7] +
            (uint64_t)a[1] * b[6] +
            (uint64_t)a[2] * b[5] +
            (uint64_t)a[3] * b[4] +
            (uint64_t)a[4] * b[3] +
            (uint64_t)a[5] * b[2] +
            (uint64_t)a[6] * b[1] +
            (uint64_t)a[7] * b[0];
    uint32_t t7 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[0] * b[8] +
            (uint64_t)a[1] * b[7] +
            (uint64_t)a[2] * b[6] +
            (uint64_t)a[3] * b[5] +
            (uint64_t)a[4] * b[4] +
            (uint64_t)a[5] * b[3] +
            (uint64_t)a[6] * b[2] +
            (uint64_t)a[7] * b[1] +
            (uint64_t)a[8] * b[0];
    uint32_t t8 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[0] * b[9] +
            (uint64_t)a[1] * b[8] +
            (uint64_t)a[2] * b[7] +
            (uint64_t)a[3] * b[6] +
            (uint64_t)a[4] * b[5] +
            (uint64_t)a[5] * b[4] +
            (uint64_t)a[6] * b[3] +
            (uint64_t)a[7] * b[2] +
            (uint64_t)a[8] * b[1] +
            (uint64_t)a[9] * b[0];
    uint32_t t9 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[1] * b[9] +
            (uint64_t)a[2] * b[8] +
            (uint64_t)a[3] * b[7] +
            (uint64_t)a[4] * b[6] +
            (uint64_t)a[5] * b[5] +
            (uint64_t)a[6] * b[4] +
            (uint64_t)a[7] * b[3] +
            (uint64_t)a[8] * b[2] +
            (uint64_t)a[9] * b[1];
    uint32_t t10 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[2] * b[9] +
            (uint64_t)a[3] * b[8] +
            (uint64_t)a[4] * b[7] +
            (uint64_t)a[5] * b[6] +
            (uint64_t)a[6] * b[5] +
            (uint64_t)a[7] * b[4] +
            (uint64_t)a[8] * b[3] +
            (uint64_t)a[9] * b[2];
    uint32_t t11 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[3] * b[9] +
            (uint64_t)a[4] * b[8] +
            (uint64_t)a[5] * b[7] +
            (uint64_t)a[6] * b[6] +
            (uint64_t)a[7] * b[5] +
            (uint64_t)a[8] * b[4] +
            (uint64_t)a[9] * b[3];
    uint32_t t12 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[4] * b[9] +
            (uint64_t)a[5] * b[8] +
            (uint64_t)a[6] * b[7] +
            (uint64_t)a[7] * b[6] +
            (uint64_t)a[8] * b[5] +
            (uint64_t)a[9] * b[4];
    uint32_t t13 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[5] * b[9] +
            (uint64_t)a[6] * b[8] +
            (uint64_t)a[7] * b[7] +
            (uint64_t)a[8] * b[6] +
            (uint64_t)a[9] * b[5];
    uint32_t t14 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[6] * b[9] +
            (uint64_t)a[7] * b[8] +
            (uint64_t)a[8] * b[7] +
            (uint64_t)a[9] * b[6];
    uint32_t t15 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[7] * b[9] +
            (uint64_t)a[8] * b[8] +
            (uint64_t)a[9] * b[7];
    uint32_t t16 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[8] * b[9] +
            (uint64_t)a[9] * b[8];
    uint32_t t17 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[9] * b[9];
    uint32_t t18 = c & 0x3FFFFFFUL; c = c >> 26;
    uint32_t t19 = c;

    c = t0 + (uint64_t)t10 * 0x3D10UL;
    t0 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t1 + (uint64_t)t10*0x400UL + (uint64_t)t11 * 0x3D10UL;
    t1 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t2 + (uint64_t)t11*0x400UL + (uint64_t)t12 * 0x3D10UL;
    t2 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t3 + (uint64_t)t12*0x400UL + (uint64_t)t13 * 0x3D10UL;
    r[3] = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t4 + (uint64_t)t13*0x400UL + (uint64_t)t14 * 0x3D10UL;
    r[4] = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t5 + (uint64_t)t14*0x400UL + (uint64_t)t15 * 0x3D10UL;
    r[5] = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t6 + (uint64_t)t15*0x400UL + (uint64_t)t16 * 0x3D10UL;
    r[6] = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t7 + (uint64_t)t16*0x400UL + (uint64_t)t17 * 0x3D10UL;
    r[7] = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t8 + (uint64_t)t17*0x400UL + (uint64_t)t18 * 0x3D10UL;
    r[8] = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t9 + (uint64_t)t18*0x400UL + (uint64_t)t19 * 0x1000003D10UL;
    r[9] = c & 0x03FFFFFUL; c = c >> 22;
    uint64_t d = t0 + c * 0x3D1UL;
    r[0] = d & 0x3FFFFFFUL; d = d >> 26;
    d = d + t1 + c*0x40;
    r[1] = d & 0x3FFFFFFUL; d = d >> 26;
    r[2] = t2 + d;
}

void secp256k1_fe_sqr_inner(const uint32_t *a, uint32_t *r) {
    uint64_t c = (uint64_t)a[0] * a[0];
    uint32_t t0 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[0]*2) * a[1];
    uint32_t t1 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[0]*2) * a[2] +
            (uint64_t)a[1] * a[1];
    uint32_t t2 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[0]*2) * a[3] +
            (uint64_t)(a[1]*2) * a[2];
    uint32_t t3 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[0]*2) * a[4] +
            (uint64_t)(a[1]*2) * a[3] +
            (uint64_t)a[2] * a[2];
    uint32_t t4 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[0]*2) * a[5] +
            (uint64_t)(a[1]*2) * a[4] +
            (uint64_t)(a[2]*2) * a[3];
    uint32_t t5 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[0]*2) * a[6] +
            (uint64_t)(a[1]*2) * a[5] +
            (uint64_t)(a[2]*2) * a[4] +
            (uint64_t)a[3] * a[3];
    uint32_t t6 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[0]*2) * a[7] +
            (uint64_t)(a[1]*2) * a[6] +
            (uint64_t)(a[2]*2) * a[5] +
            (uint64_t)(a[3]*2) * a[4];
    uint32_t t7 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[0]*2) * a[8] +
            (uint64_t)(a[1]*2) * a[7] +
            (uint64_t)(a[2]*2) * a[6] +
            (uint64_t)(a[3]*2) * a[5] +
            (uint64_t)a[4] * a[4];
    uint32_t t8 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[0]*2) * a[9] +
            (uint64_t)(a[1]*2) * a[8] +
            (uint64_t)(a[2]*2) * a[7] +
            (uint64_t)(a[3]*2) * a[6] +
            (uint64_t)(a[4]*2) * a[5];
    uint32_t t9 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[1]*2) * a[9] +
            (uint64_t)(a[2]*2) * a[8] +
            (uint64_t)(a[3]*2) * a[7] +
            (uint64_t)(a[4]*2) * a[6] +
            (uint64_t)a[5] * a[5];
    uint32_t t10 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[2]*2) * a[9] +
            (uint64_t)(a[3]*2) * a[8] +
            (uint64_t)(a[4]*2) * a[7] +
            (uint64_t)(a[5]*2) * a[6];
    uint32_t t11 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[3]*2) * a[9] +
            (uint64_t)(a[4]*2) * a[8] +
            (uint64_t)(a[5]*2) * a[7] +
            (uint64_t)a[6] * a[6];
    uint32_t t12 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[4]*2) * a[9] +
            (uint64_t)(a[5]*2) * a[8] +
            (uint64_t)(a[6]*2) * a[7];
    uint32_t t13 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[5]*2) * a[9] +
            (uint64_t)(a[6]*2) * a[8] +
            (uint64_t)a[7] * a[7];
    uint32_t t14 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[6]*2) * a[9] +
            (uint64_t)(a[7]*2) * a[8];
    uint32_t t15 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[7]*2) * a[9] +
            (uint64_t)a[8] * a[8];
    uint32_t t16 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)(a[8]*2) * a[9];
    uint32_t t17 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + (uint64_t)a[9] * a[9];
    uint32_t t18 = c & 0x3FFFFFFUL; c = c >> 26;
    uint32_t t19 = c;

    c = t0 + (uint64_t)t10 * 0x3D10UL;
    t0 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t1 + (uint64_t)t10*0x400UL + (uint64_t)t11 * 0x3D10UL;
    t1 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t2 + (uint64_t)t11*0x400UL + (uint64_t)t12 * 0x3D10UL;
    t2 = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t3 + (uint64_t)t12*0x400UL + (uint64_t)t13 * 0x3D10UL;
    r[3] = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t4 + (uint64_t)t13*0x400UL + (uint64_t)t14 * 0x3D10UL;
    r[4] = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t5 + (uint64_t)t14*0x400UL + (uint64_t)t15 * 0x3D10UL;
    r[5] = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t6 + (uint64_t)t15*0x400UL + (uint64_t)t16 * 0x3D10UL;
    r[6] = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t7 + (uint64_t)t16*0x400UL + (uint64_t)t17 * 0x3D10UL;
    r[7] = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t8 + (uint64_t)t17*0x400UL + (uint64_t)t18 * 0x3D10UL;
    r[8] = c & 0x3FFFFFFUL; c = c >> 26;
    c = c + t9 + (uint64_t)t18*0x400UL + (uint64_t)t19 * 0x1000003D10UL;
    r[9] = c & 0x03FFFFFUL; c = c >> 22;
    uint64_t d = t0 + c * 0x3D1UL;
    r[0] = d & 0x3FFFFFFUL; d = d >> 26;
    d = d + t1 + c*0x40;
    r[1] = d & 0x3FFFFFFUL; d = d >> 26;
    r[2] = t2 + d;
}


void secp256k1_fe_mul(secp256k1_fe_t *r, const secp256k1_fe_t *a, const secp256k1_fe_t *b) {
#ifdef VERIFY
    assert(a->magnitude <= 8);
    assert(b->magnitude <= 8);
    r->magnitude = 1;
    r->normalized = 0;
#endif
    secp256k1_fe_mul_inner(a->n, b->n, r->n);
}

void secp256k1_fe_sqr(secp256k1_fe_t *r, const secp256k1_fe_t *a) {
#ifdef VERIFY
    assert(a->magnitude <= 8);
    r->magnitude = 1;
    r->normalized = 0;
#endif
    secp256k1_fe_sqr_inner(a->n, r->n);
}

void secp256k1_gej_double(secp256k1_gej_t *r, secp256k1_gej_t *a) {
    secp256k1_fe_t t5 = a->y;
    secp256k1_fe_normalize(&t5);
    if (a->infinity || secp256k1_fe_is_zero(&t5)) {
        r->infinity = 1;
        return;
    }

    secp256k1_fe_t t1,t2,t3,t4;
    secp256k1_fe_mul(&r->z, &t5, &a->z);
    secp256k1_fe_mul_int(&r->z, 2);       // Z' = 2*Y*Z (2)
    secp256k1_fe_sqr(&t1, &a->x);
    secp256k1_fe_mul_int(&t1, 3);         // T1 = 3*X^2 (3)
    secp256k1_fe_sqr(&t2, &t1);           // T2 = 9*X^4 (1)
    secp256k1_fe_sqr(&t3, &t5);
    secp256k1_fe_mul_int(&t3, 2);         // T3 = 2*Y^2 (2)
    secp256k1_fe_sqr(&t4, &t3);
    secp256k1_fe_mul_int(&t4, 2);         // T4 = 8*Y^4 (2)
    secp256k1_fe_mul(&t3, &a->x, &t3);    // T3 = 2*X*Y^2 (1)
    r->x = t3;
    secp256k1_fe_mul_int(&r->x, 4);       // X' = 8*X*Y^2 (4)
    secp256k1_fe_negate(&r->x, &r->x, 4); // X' = -8*X*Y^2 (5)
    secp256k1_fe_add(&r->x, &t2);         // X' = 9*X^4 - 8*X*Y^2 (6)
    secp256k1_fe_negate(&t2, &t2, 1);     // T2 = -9*X^4 (2)
    secp256k1_fe_mul_int(&t3, 6);         // T3 = 12*X*Y^2 (6)
    secp256k1_fe_add(&t3, &t2);           // T3 = 12*X*Y^2 - 9*X^4 (8)
    secp256k1_fe_mul(&r->y, &t1, &t3);    // Y' = 36*X^3*Y^2 - 27*X^6 (1)
    secp256k1_fe_negate(&t2, &t4, 2);     // T2 = -8*Y^4 (3)
    secp256k1_fe_add(&r->y, &t2);         // Y' = 36*X^3*Y^2 - 27*X^6 - 8*Y^4 (4)
    r->infinity = 0;
}

void secp256k1_gej_add(__global secp256k1_gej_t *pr, const secp256k1_gej_t *pa, __global const secp256k1_gej_t *pb) {
    secp256k1_gej_t b, r;
    b = *pb;
    secp256k1_gej_add_local(&r, pa, &b);
    *pr = r;
}

void secp256k1_gej_add_local(secp256k1_gej_t *pr, const secp256k1_gej_t *pa, const secp256k1_gej_t *pb) {
    if (pa->infinity) {
        *pr = *pb;
        return;
    }
    if (pb->infinity) {
        *pr = *pa;
        return;
    }
    
    secp256k1_gej_t a = *pa;
    secp256k1_gej_t b = *pb;
    secp256k1_gej_t r;
    
    r.infinity = 0;
    secp256k1_fe_t z22; secp256k1_fe_sqr(&z22, &b.z);
    secp256k1_fe_t z12; secp256k1_fe_sqr(&z12, &a.z);
    secp256k1_fe_t u1; secp256k1_fe_mul(&u1, &a.x, &z22);
    secp256k1_fe_t u2; secp256k1_fe_mul(&u2, &b.x, &z12);
    secp256k1_fe_t s1; secp256k1_fe_mul(&s1, &a.y, &z22); secp256k1_fe_mul(&s1, &s1, &b.z);
    secp256k1_fe_t s2; secp256k1_fe_mul(&s2, &b.y, &z12); secp256k1_fe_mul(&s2, &s2, &a.z);
    secp256k1_fe_normalize(&u1);
    secp256k1_fe_normalize(&u2);
    if (secp256k1_fe_equal(&u1, &u2)) {
        secp256k1_fe_normalize(&s1);
        secp256k1_fe_normalize(&s2);
        if (secp256k1_fe_equal(&s1, &s2)) {
            secp256k1_gej_double(&r, &a);
        } else {
            r.infinity = 1;
        }
        return;
    }
    secp256k1_fe_t h; secp256k1_fe_negate(&h, &u1, 1); secp256k1_fe_add(&h, &u2);
    secp256k1_fe_t i; secp256k1_fe_negate(&i, &s1, 1); secp256k1_fe_add(&i, &s2);
    secp256k1_fe_t i2; secp256k1_fe_sqr(&i2, &i);
    secp256k1_fe_t h2; secp256k1_fe_sqr(&h2, &h);
    secp256k1_fe_t h3; secp256k1_fe_mul(&h3, &h, &h2);
    secp256k1_fe_mul(&r.z, &a.z, &b.z); secp256k1_fe_mul(&r.z, &r.z, &h);
    secp256k1_fe_t t; secp256k1_fe_mul(&t, &u1, &h2);
    r.x = t; secp256k1_fe_mul_int(&r.x, 2); secp256k1_fe_add(&r.x, &h3); secp256k1_fe_negate(&r.x, &r.x, 3); secp256k1_fe_add(&r.x, &i2);
    secp256k1_fe_negate(&r.y, &r.x, 5); secp256k1_fe_add(&r.y, &t); secp256k1_fe_mul(&r.y, &r.y, &i);
    secp256k1_fe_mul(&h3, &h3, &s1); secp256k1_fe_negate(&h3, &h3, 1);
    secp256k1_fe_add(&r.y, &h3);
    
    *pr = r;
    }


void secp256k1_gej_add_ge(secp256k1_gej_t *r, const secp256k1_gej_t *a, const secp256k1_ge_t *b) {
    if (a->infinity) {
        r->infinity = b->infinity;
        r->x = b->x;
        r->y = b->y;
        secp256k1_fe_set_int(&r->z, 1);
        return;
    }
    if (b->infinity) {
        *r = *a;
        return;
    }
    r->infinity = 0;
    secp256k1_fe_t z12; secp256k1_fe_sqr(&z12, &a->z);
    secp256k1_fe_t u1 = a->x; secp256k1_fe_normalize(&u1);
    secp256k1_fe_t u2; secp256k1_fe_mul(&u2, &b->x, &z12);
    secp256k1_fe_t s1 = a->y; secp256k1_fe_normalize(&s1);
    secp256k1_fe_t s2; secp256k1_fe_mul(&s2, &b->y, &z12); secp256k1_fe_mul(&s2, &s2, &a->z);
    secp256k1_fe_normalize(&u1);
    secp256k1_fe_normalize(&u2);
    if (secp256k1_fe_equal(&u1, &u2)) {
        secp256k1_fe_normalize(&s1);
        secp256k1_fe_normalize(&s2);
        if (secp256k1_fe_equal(&s1, &s2)) {
            secp256k1_gej_double(r, a);
        } else {
            r->infinity = 1;
        }
        return;
    }
    secp256k1_fe_t h; secp256k1_fe_negate(&h, &u1, 1); secp256k1_fe_add(&h, &u2);
    secp256k1_fe_t i; secp256k1_fe_negate(&i, &s1, 1); secp256k1_fe_add(&i, &s2);
    secp256k1_fe_t i2; secp256k1_fe_sqr(&i2, &i);
    secp256k1_fe_t h2; secp256k1_fe_sqr(&h2, &h);
    secp256k1_fe_t h3; secp256k1_fe_mul(&h3, &h, &h2);
    r->z = a->z; secp256k1_fe_mul(&r->z, &r->z, &h);
    secp256k1_fe_t t; secp256k1_fe_mul(&t, &u1, &h2);
    r->x = t; secp256k1_fe_mul_int(&r->x, 2); secp256k1_fe_add(&r->x, &h3); secp256k1_fe_negate(&r->x, &r->x, 3); secp256k1_fe_add(&r->x, &i2);
    secp256k1_fe_negate(&r->y, &r->x, 5); secp256k1_fe_add(&r->y, &t); secp256k1_fe_mul(&r->y, &r->y, &i);
    secp256k1_fe_mul(&h3, &h3, &s1); secp256k1_fe_negate(&h3, &h3, 1);
    secp256k1_fe_add(&r->y, &h3);
}
