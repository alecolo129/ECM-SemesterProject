load("parameters_anomalous.sage")
load("ecm_anomalous.sage")
load("ecm.sage")
import time

"""
Test ecm with known j invariants on the RSA moduli contained in "parameters_anomalous.sage"
"""
def test_ecm_known_j(Ns, curves):
    print("----------------------------------------------------------------------------------------------")
    count = 0
    elapsed = time.time()
    for i in range(len(Ns)):
        print(f"Factoring N{i}:")
        t = time.time()
        d = ecm_factor(n=Ns[i], B=0, curves = curves, anomalous=True)
        if d:
            print(f'\x1b[34m{Ns[i]} = {d} * {Ns[i]//d}\x1b[0m')
            count+=1
    elapsed = time.time()-elapsed
    print(f"ECM anomalous known j -> percentage of factoring = {(count/len(Ns)).n()}, elapsed time:{elapsed}, avg: {elapsed/len(Ns)}")
    print("----------------------------------------------------------------------------------------------")


"""
Test ecm with unknown j invariants on the RSA moduli contained in "parameters_anomalous.sage"
"""
def test_ecm_unknown_j(rep = 5, start=0, end=len(Ns_unknown)):
    print("----------------------------------------------------------------------------------------------")
    count = 0
    elapsed = time.time()
    for i in range(start, end):
        print(f"Factoring N{i}:")
        t = time.time()
        N = Ns_unknown[i][0]
        D = Ns_unknown[i][1]
        d = ECM_anomalous(N, D)
        if d:
            count+=1
    elapsed = time.time()-elapsed
    print(f"ECM anomalous unknown j -> percentage of factoring = {count/(end-start).n()}, elapsed time:{elapsed}, avg: {elapsed/(end-start)}")
    print("----------------------------------------------------------------------------------------------")

if __name__ == "__main__":
    test_ecm_known_j(Ns_known_11, CURVES_ANOMALOUS_11)
    test_ecm_unknown_j()

