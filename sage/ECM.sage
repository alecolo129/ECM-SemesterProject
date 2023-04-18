from sage.all_cmdline import *
from sage.arith.misc import xgcd
import time
import parameters

def add(a, b, P, Q, n):
    """
    We assume that (xp,yp,zp), (xq,yq,zq) have already been casted to Integers(n)
    """
    xp, yp, zp = P[0], P[1], P[2]
    xq, yq, zq = Q[0], Q[1], Q[2]
    
    """
    if(4*(a**3)+27*(b**2)==0): 
        raise Exception("Invalid parameters: the curve has to be non-singular") 
    if((zp * power_mod(yp,2,n))%n != (zp*power_mod(xp,3,n) + power_mod(zp,2,n)*a*xp + power_mod(zp,3,n)*b)%n or (zq*power_mod(yq,2,n))%n != (zq*power_mod(xq,3,n) + power_mod(zq,2,n)*a*xq + power_mod(zq,3,n)*b)%n):    
        raise Exception("Invalid parameters: the two points need to stay on the curve") 
    """
    
    if((xp,yp,zp) == (0,1,0)):
        return Q
    if((xq,yq,zq) == (0,1,0)):
        return P
    if (xp==xq and yp==-yq):
        return (0,1,0)
    
    d, a1, lamb = 0,0,0
    if (xp!=xq):
        d,a1,a2 = xgcd(int(xp-xq), int(n))
        if (d != 1):
            raise Exception(d)
        lamb = (yp-yq)*a1
    else:
        d,a1,a2 = xgcd(int(2*yp), int(n))
        if (d != 1):
            raise Exception(d)
        lamb = (3*pow(xp,2)+a)*a1

    #print("lamb:", pow(lamb,2)-xq)   
    xr = pow(lamb,2)-xp-xq
    yr = lamb * (xp-xr) - yp
    return (xr,yr,1)

def double_add(a,b,P,k,n):
    Zn = Integers(n)
    R = (Zn(0),Zn(1),Zn(0))
    P = (Zn(P[0]),Zn(P[1]),Zn(P[2]))
    nb = k.nbits()-1
    for i in range(nb,-1,-1):
        R = add(a,b,R,R,n)
        if k & (1<<(i)):
            R = add(a,b,R,P,n)
    return R
    

def ECMfactor(curves, B, n, w=30):
    """
    Tries to factor the RSA modulus n using the ECM algorithm

    input:
    - curves: the list of elliptic curves
    - B: bound
    - n: RSA modulus 
    - w: maximum time to spend on each curve

    ouput:
    - A non trivial factor p or None
    """
    k = 1
    ks = []

    #populate list ks with all pi^alphai
    for p in prime_range(B):
        ki = p
        k *= p
        while ki<B:
            ki *= p
            if(ki<=B):
                k *= p
        ks.append(k)
        k=1
    
    #try to factor with all the input curves until either we get a non-trivial factor or we fail
    for i in range(len(curves)):
        P = curves[i][1]
        a, b = (curves[i][0][0]%n), (curves[i][0][1]%n)
        print(f"Trying curve {i} ...")
        c=0
        t = time.time()
        for ki in ks:
            c+=1
            try:
                P = double_add(a,b,P,ki,n)    
            except Exception as e:
                return e.args[0] 
            if(time.time()-t>w):
                break
    return None
    


Ns = parameters.Ns
curves = parameters.CURVES_SUPERSINGULAR

B = Ns[5][1]
N = Ns[5][0]

print(ECMfactor(curves, B, N))
