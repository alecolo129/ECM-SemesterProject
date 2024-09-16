import sys
load("parameters_anomalous.sage")

def gen(bits, D=None):
    if D!= None: 
        if D>0 or D%4!=1:
            print("D has to be negative and satisfy D = 1 mod 4") 
            return None
        else:
            b = (-D+1)//4
            while True:
                m = randint(2^bits,2^(bits+1))
                p = -D*m*(m+1)+b
                if is_prime(p):
                    return (p,D)
    else:
        while True:
            for D in Ds:
                m = randint(2^bits,2^(bits+1))
                p = -D*m*(m+1)+(-D+1)//4
                if is_prime(p):
                    return (p,D)


if __name__ == "__main__":
    if len(sys.argv) not in [2, 3]:
        raise ValueError("Usage: sage generate.sage <size> or sage generate.sage <size> <discriminant>")
    else:
        bits = int(sys.argv[1])
        prime_bits = bits//2
        D = None
        if len(sys.argv) > 2: 
            D = int(sys.argv[2])
        for i in range(10):
            p, d = gen(prime_bits//2, D)
            print(f"({p*random_prime(2^(prime_bits-1),2^(prime_bits))}, {d})")