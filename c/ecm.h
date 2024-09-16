#ifndef _ECM_H_
#define _ECM_H_
#include "gmp.h"
#include "stdint.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "curves.h"
#define DIM 82025
#define TRUE 1
#define FALSE 0
#define BLUE  "\x1B[34m"
#define RESET "\x1B[0m"

mpz_t lamb, a1, a2;
int MSB(const uintmax_t k){
    int i = 63;
    while(i>0){
        if((1ull<<i) & k)
            return i;
        i--;
    }
    return 0;
}
void erathostenes(uintmax_t primes[], const uintmax_t B){
    char numbers[B+1];
    for(int i = 0; i<B; i++) numbers[i]=1;
    uintmax_t jp = 0;
    for(uintmax_t i = 2; i<=B; i++){
        if(numbers[i]){
            primes[jp++]=i;
            for(uintmax_t j = i*i; j<=B; j+=i){
                numbers[j] = 0; 
            }
        }
    }
}

int add(const mpz_t a, Point *R, const Point *P, const Point *Q, const mpz_t N){

    if(P->z == 0){
        mpz_set(R->x, Q->x), mpz_set(R->y, Q->y), R->z=Q->z;
        return 0;
    }
    if(Q->z == 0){
        mpz_set(R->x, P->x), mpz_set(R->y, P->y), R->z=P->z;
        return 0;
    }

    if(!mpz_cmp(P->x,Q->x) && (mpz_cmp(P->y, Q->y) && !mpz_cmpabs(P->y,Q->y))){
        mpz_set_ui(R->x,0), mpz_set_ui(R->y,1), R->z=0;
        return 0;
    }


    //case xp!=xq
    if(mpz_cmp(P->x,Q->x)){
        mpz_sub(a1, P->x, Q->x);
        mpz_gcdext(a2, a1, a2, a1, N);
        if(mpz_cmp_ui(a2,1)){
            gmp_printf("%sp = %Zd%s\n",BLUE, a2,RESET);
            return 1;
        }
        mpz_sub(lamb, P->y, Q->y);
        mpz_mul(lamb, lamb, a1);
        mpz_mod(lamb, lamb, N);
    }
    else{
        // yr = yp*2
        mpz_mul_ui(a2, P->y,2); 
        mpz_gcdext(a2, a1, a2, a2, N);
        if(mpz_cmp_ui(a2,1)){
            gmp_printf("%sp = %Zd%s\n",BLUE,a2,RESET);
            return 1;
        }

        //lamb = ((3*xp^2+a)*(2yp)^-1) mod N
        mpz_pow_ui(lamb, P->x, 2); mpz_mul_ui(lamb, lamb, 3);
        //(mpz_sgn(a) > 0) ? mpz_add(lamb, lamb, a) : mpz_sub_ui(lamb, lamb, abs(a)); 
        mpz_add(lamb,lamb,a);
        mpz_mul(lamb, lamb, a1);
        mpz_mod(lamb, lamb, N);
    }

    //xr = (lamb*lamb-xp-xq) mod N
    mpz_mul(a1, lamb, lamb), mpz_sub(a1, a1, Q->x), mpz_sub(a1, a1, P->x); 
    mpz_mod(a1, a1, N);  

    //yr = (lamb * (xp-xr)-yp) mod N;
    mpz_sub(a2, P->x, a1); mpz_mul(a2,lamb,a2); mpz_sub(a2, a2, P->y);
    mpz_mod(a2, a2, N);
    R->z = 1;
    mpz_set(R->x,a1), mpz_set(R->y, a2);
    return 0;
}

int double_add(const mpz_t a, Point *R, const Point *P, const void* kp, const mpz_t N, uint8_t anomalous){
    mpz_set_ui(R->x,0); mpz_set_ui(R->y,1); R->z=0;
    int i;
    if(anomalous){ 
        i = mpz_sizeinbase(*((mpz_t*) kp), 2);
        for(; i>=0; i--){
            add(a,R,R,R,N);
            if( mpz_tstbit(*((mpz_t*) kp),i) ){
                if(add(a,R,R,P,N)) return 1;
            }
        }
    }
    else{
        i = MSB(*((uintmax_t* ) kp));
        for(; i>=0; i--){
            add(a,R,R,R,N);
            if( *((uintmax_t*) kp) & (1ull<<i) ){
                if(add(a,R,R,P,N)) return 1;
            }
        }
    }
    return 0;
}

int ECM_anomalous(Curve* C, const mpz_t N){
    size_t bits = mpz_sizeinbase(N, 2);
    mpz_inits(lamb,a1,a2,NULL);
    Point R; mpz_init(R.x), mpz_init_set_ui(R.y,1), R.z=0;
    Point* P = &(C->point);
    if(double_add(C->a, &R,P, N, N,TRUE)) return 1;
    mpz_clears(lamb,a1,a2,NULL);
    return -1;
}

void ECM(uintmax_t B, Curve* C, mpz_t N){
    mpz_inits(lamb,a1,a2,NULL);
    uintmax_t ks[DIM], primes[DIM];
    uintmax_t k;
    erathostenes(primes, B);
    for(uintmax_t i = 0; i<DIM; i++){
        uintmax_t ki = primes[i];
        k=primes[i];
        while(ki < B && ki>0){
            ki *= primes[i];
            if(ki<=B)
                k *= primes[i];
        } 
        ks[i] = k;
    } printf("%lu, %lu\n",ks[0],ks[1]);

    Point R; mpz_init(R.x), mpz_init_set_ui(R.y,1), R.z=0;
    Point* P = &(C->point);

    for(uintmax_t j = 0; j<DIM; j++){
        if(j%2==0) {
            if(double_add(C->a,&R,&(C->point),&ks[j],N,FALSE)) 
                return;
        }
        else {
            if(double_add(C->a,&(C->point),&R,&ks[j],N,FALSE)) 
                return;
        }
    }
    mpz_clears(lamb,a1,a2,NULL);
    return;
}




#endif