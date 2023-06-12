from sage.rings.finite_rings.integer_mod import (lucas, lucas_q1)

def williams_factor(N, B1, rep = 5):
    Zn = Integers(N)
    p0s = [Integers(2**10).random_element() for i in range(rep)]
    ri = 1
    Ris = []
    
    tresh = 2**800
    count = 1
    
    for p in prime_range(B1):
        ri *= p
        tmp = p
        while(tmp < B1):
            tmp *= p
            if(tmp<B1):
                ri *= p

        #we limit the calls to lucas_q1 t0 improve efficiency
        if(ri>tresh):
            count+=1
            p0s[0] = lucas_q1(ri, Zn(p0s[0])) #lucas(rj, lucas(ri,p0)) = lucas(ri*rj, p0)
            d = gcd(p0s[0]-2,N)
            if(d!=0 and d!=1 and d!=N): 
                return d
            Ris.append(ri)
            ri=1

    #if the algo failed with the first p0 we try with the others (may have (p0/delta) = 1)
    for i in range(1, rep):
        for j in range(len(Ris)):
            p0s[i] = lucas_q1(Ris[j], Zn(p0s[i]))
            d = gcd(p0s[i]-2, N)
            if(d!=0 and d!=1 and d!=N):
                return d