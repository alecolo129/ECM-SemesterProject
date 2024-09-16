from sage.arith.misc import xgcd
import time


def add(a, b, P, Q, n):

    # We assume that (xp,yp,zp), (xq,yq,zq) have already been casted to Integers(n)
    xp, yp, zp = P[0], P[1], P[2]
    xq, yq, zq = Q[0], Q[1], Q[2]

    if ((xp, yp, zp) == (0, 1, 0)):
        return Q
    if ((xq, yq, zq) == (0, 1, 0)):
        return P
    if (xp == xq and yp == -yq):
        return (0, 1, 0)

    d, a1, lamb = 0, 0, 0
    if (xp != xq):
        d, a1, a2 = xgcd(int(xp-xq), int(n))
        if (d != 1):
            raise Exception(int(d))
        lamb = (yp-yq)*a1
    else:
        d, a1, a2 = xgcd(int(2*yp), int(n))
        if (d != 1):
            raise Exception(int(d))
        lamb = (3*pow(xp, 2)+a)*a1

    xr = pow(lamb, 2)-xp-xq
    yr = lamb * (xp-xr) - yp
    return (xr, yr, 1)


def double_add(a, b, P, k, n):
    Zn = Integers(n)
    R = (Zn(0), Zn(1), Zn(0))
    P = (Zn(P[0]), Zn(P[1]), Zn(P[2]))
    nb = k.nbits()-1
    for i in range(nb, -1, -1):
        R = add(a, b, R, R, n)
        if k & (1 << (i)):
            R = add(a, b, R, P, n)
    return R




def ecm_factor(n, B, curves, anomalous=False, w=30, test=False):
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
    ks = []
    # for anomalous curves we just multiply the point by n
    if anomalous:
        ks.append(Integer(n))
    # for other curves we use the standard ECM method
    else:
        k = 1
        # populate list ks with all pi^alphai
        for p in prime_range(B):
            ki = p
            k *= p
            while ki < B:
                ki *= p
                if (ki <= B):
                    k *= p
            ks.append(k)
            k = 1
    # try to factor with all the input curves
    for i in range(len(curves)):
        P = curves[i][1]
        if len(P) == 2:
            P.append(1)
        a, b = curves[i][0][0] % n, curves[i][0][1] % n
        if (test):
            print(f"Trying curve {i} ...")
        t = time.time()
        for ki in ks:
            try:
                P = double_add(a, b, P, ki, n)
            # we found a non trivial factor
            except Exception as e:
                if type(e.args[0]) == int:
                    return e.args[0]
                else:
                    raise e

            # max time exceeded, try another curve
            if (time.time()-t > w):
                break
    return None