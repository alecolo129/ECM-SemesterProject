load("parameters.sage")
load("williams.sage")
load("ecm.sage")
import time

"""
Test the supersingular ECM on the RSA moduli contained in "paremeters.sage"
"""
def test_ecm(curves, start=0, end=len(Ns), w=30):
    count = 0
    elapsed = time.time()
    for i in range(start,end):
        t = time.time()
        N = Ns[i][0]
        B = Ns[i][1]
        d = ecm_factor(n=N, B=B, curves=curves, anomalous = False, w=w, test=True)
        if d:
            print(f"N{i}: {time.time()-t}s, found p={d}\n")
            count+=1
    elapsed = time.time()-elapsed
    print(f"Williams p+1 -> percentage of factoring = {count/(end-start).n()}, elapsed time:{elapsed}, avg: {elapsed/(end-start)}")
    print("----------------------------------------------------------------------------------------------")

"""
Test Williams p+1 on the RSA moduli contained in "paremeters.sage"
"""
def test_williams(rep = 5, start=0, end=len(Ns)):
    count = 0
    elapsed = time.time()
    for i in range(start, end):
        t = time.time()
        B = Ns[i][1]
        N = Ns[i][0]
        d = williams_factor(N,B,rep)
        if d:
            print(f"N{i}: {time.time()-t}s, found p={d}\n")
            count+=1
    elapsed = time.time()-elapsed
    print(f"Williams p+1 -> percentage of factoring = {count/(end-start).n()}, elapsed time:{elapsed}, avg: {elapsed/(end-start)}")
    print("----------------------------------------------------------------------------------------------")

if __name__ == "__main__":
    test_williams(10,0,2)
    #test_ecm(Ns, CURVES_SUPERSINGULAR, start=0, end=2)

