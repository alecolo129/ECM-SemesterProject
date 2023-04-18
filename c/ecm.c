#include "stdio.h"
#include "stdlib.h"
#include "gmp.h"
#include "time.h"
#define DIM 82025
struct point{
    mpz_t x;
    mpz_t y;
    int z;
};
typedef struct point Point;

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

int add(const int a, const int b, Point *R, const Point *P, const Point *Q, const mpz_t *N){

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
        mpz_gcdext(a2, a1, a2, a1, *N);
        if(mpz_cmp_ui(a2,1)){
            gmp_printf("%Zd\n",a2);
            return 1;
        }
        mpz_sub(lamb, P->y, Q->y);
        mpz_mul(lamb, lamb, a1);
        mpz_mod(lamb, lamb, *N);
    }
    else{
        // yr = yp*2
        mpz_mul_ui(a2, P->y,2); 
        mpz_gcdext(a2, a1, a2, a2, *N);
        if(mpz_cmp_ui(a2,1)){
            gmp_printf("%Zd\n",a2);
            return 1;
        }

        //lamb = ((3*xp^2+a)*(2yp)^-1) mod N
        mpz_pow_ui(lamb, P->x, 2); mpz_mul_ui(lamb, lamb, 3);
        (a>0) ? mpz_add_ui(lamb,  lamb, a) : mpz_sub_ui(lamb, lamb, abs(a)); 
        mpz_mul(lamb, lamb, a1);
        mpz_mod(lamb, lamb, *N);
    }

    //xr = (lamb*lamb-xp-xq) mod N
    mpz_mul(a1, lamb, lamb), mpz_sub(a1, a1, Q->x), mpz_sub(a1, a1, P->x); 
    mpz_mod(a1, a1, *N);  

    //yr = (lamb * (xp-xr)-yp) mod N;
    mpz_sub(a2, P->x, a1); mpz_mul(a2,lamb,a2); mpz_sub(a2, a2, P->y);
    mpz_mod(a2, a2, *N);
    R->z = 1;
    mpz_set(R->x,a1), mpz_set(R->y, a2);
    return 0;
}

int double_add(const int a, const int b, Point *R, const Point *P, const uintmax_t k, const mpz_t *N){
    mpz_set_ui(R->x,0); mpz_set_ui(R->y,1); R->z=0;
    for(int i = MSB(k); i>=0; i--){
        add(a,b,R,R,R,N);
        if(k & (1ull<<i)){
            if(add(a,b,R,R,P,N)) return 1;
        }
    }
    return 0;
}

void ECM(uintmax_t B, mpz_t N){
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
    for(int i = 7; i<=7; i++){
        int a = -152, b= -722;
        Point P; mpz_init_set_si(P.x,19), mpz_init_set_si(P.y,57), P.z=1;
        Point R; mpz_init(R.x), mpz_init_set_ui(R.y,1), R.z=0;
        mpz_t N; mpz_init_set_str(N,"1277363318506817942182519338423573481457994335556110607873995437465702672936563589384714346315827851651437147635758847500078874864149284979043565589635766201408749733227105731204907857990552053734202218060161764749808536025189091939157036622975191144016990167639754671557986776742905325413958199315014282356144176368119",10);
        int count = 0;
        for(uintmax_t j = 0; j<DIM; j++){
            printf("%d\n",count);
            count++;
            if(j%2==0) {if(double_add(a,b,&R,&P,ks[j],&N)) return;}
            else {if(double_add(a,b,&P,&R,ks[j],&N)) return;}
        }
    }
    mpz_clears(lamb,a1,a2,NULL);
    return;
}

int main(int argc, char** argv){
    Point R; mpz_init(R.x), mpz_init_set_ui(R.y,1), R.z=0;
    Point P; mpz_init_set_si(P.x,19), mpz_init_set_si(P.y,57), P.z=1;

     
    int a = -152, b= -722;
    clock_t t1,t2;
    mpz_t N; mpz_init_set_str(N,"1277363318506817942182519338423573481457994335556110607873995437465702672936563589384714346315827851651437147635758847500078874864149284979043565589635766201408749733227105731204907857990552053734202218060161764749808536025189091939157036622975191144016990167639754671557986776742905325413958199315014282356144176368119",10);
    int count=0;
    t1 = clock();
    uintmax_t B = 1048576;
    ECM(B,N);
    t2=clock();
    printf("%f\n",((double)t2-t1)/CLOCKS_PER_SEC);
    gmp_printf("(%Zd : %Zd : %d)\n",R.x, R.y, R.z);

    
    return 0;
}