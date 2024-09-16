import re
import time

##################################################################################################
#This code re-adapted from https://hxp.io/blog/86/hxp-CTF-2021-f_cktoring-writeup/
##################################################################################################


# Sage can't invert in ℤ[x]/(n,f(x)); we use this workaround
def inv(f):
    ff = f.lift().change_ring(QQ)
    gg = f.parent().modulus().change_ring(QQ)
    hh = xgcd(ff,gg)            # compute xgcd over ℚ
    assert hh[0] == 1          
    return f.parent()(hh[1])    # reduce back

# X1 is P-Q, X2 is P, X3 is Q, returns 2P, P+Q
def xDBLADD(a,b,X1,X2,X3):
    Z1 = Z2 = Z3 = 1
    if X1 == (): X1,Z1 = 1,0
    if X2 == (): X2,Z2 = 1,0
    if X3 == (): X3,Z3 = 1,0
    X4 = (X2^2-a*Z2^2)^2-8*b*X2*Z2^3
    Z4 = 4*(X2*Z2*(X2^2+a*Z2^2)+b*Z2^4)
    R = 2*(X2*Z3+X3*Z2)*(X2*X3+a*Z2*Z3)+4*b*Z2^2*Z3^2
    S = (X2*Z3-X3*Z2)^2
    X5 = R-S*X1
    Z5 = S
    return (X4*inv(Z4) if Z4 else ()), (X5*inv(Z5) if Z5 else ())

def xMUL(a,b,n,P):
    Q,PQ = (),P
    for bit in map(int,f'{n:b}'):
        if bit:
            P,Q = xDBLADD(a,b,PQ,P,Q)
        else:
            Q,P = xDBLADD(a,b,PQ,Q,P)
    return Q

def ECM_anomalous(n, d):
    print(f"trying d = {d}")
    n = ZZ(n)
    R.<x> = Zmod(n)[]
    S.<J> = R.quotient(hilbert_class_polynomial(d))  # J is a symbolic root of H_{-d} modulo n
    try:
        a0, b0 = 27*J*inv(1728-J)/4, -27*J*inv(1728-J)/4  # y^2 = x^3 + a(j)*x - a(j) has j-invariant J
        for i in range(1,10):
            P = J.parent().random_element()     # random point
            Q = xMUL(a0,b0,n, P)
    except ZeroDivisionError as e:      # probably found a divisor!
        p = ZZ(re.findall('[0-9]+', e.args[0])[0])
        p = gcd(p,n)
        if 1 < p < n:
            print(f'\x1b[34m{n} = {p} * {n//p}\x1b[0m'); 
            return p
    return None

