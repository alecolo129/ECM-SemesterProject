from sage.all_cmdline import *
from sage.rings.finite_rings.integer_mod import (lucas, lucas_q1)
import timeit
import time
import parameters

def WilliamsFactor(N, B1, rep = 5)->(int or None):
    """
    Given a number N having a prime divisor p s.t. p+1 (or p-1) has only smooth factors, returns p using the Williams p+1 method

    Input:
    - N = p*q
    - B1 = bound on the divisors of p+1
    - rep = maximum number of attempts for getting (∆i/p)=-1
    
    Output:
    - p if we found a non trivial divisor
    - None otherwise
    """
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



Ns = parameters.Ns
count = 0
factors = 10
t = time.time()
for i in range(factors):
    B = Ns[i][1]
    N = Ns[i][0]
    d = WilliamsFactor(N,B,10)
    print(f"N{i}: found p={d}\n")
    if d:
        count+=1
t = time.time()-t
print(f"Percentage of factoring = {count/factors}, time:{t}")


"""
exec = timeit.repeat(lambda: FastWilliamsFactor(11080647170605001519502158144114954921368208695151749530900403255465715110639267780896952621957060195532806935465298224407255026318244226256558541691214067391855632936226487891032236023151690791942869930919103174056622014325980098920448833729898937122816152043557244249222118067756023700844761380689528189672358643704955606390078697723770270630662506388207121079748294251031125144819655633639230913121176977912367909438802702985598586326338712846946612428657752526335430463531540373677470601795584252531593165495739060424328876458405428719028848596073974060041115720780603068753414506260138812675247333883397880248717, 2**20, 4), number=1, repeat=1)
print("Total time fastWill: {}".format(sum(exec)))
"""