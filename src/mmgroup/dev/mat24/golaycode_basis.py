from mmgroup.dev.mat24 import Mat24

basis = Mat24.basis[12:]


def print_vector(v, final =  0):    
    for j in range(24):
        print((b >> j) & 1, end = "")
        if j %4 == 3: print(" & ", end = " ")
    if not final: print(r"\\")
    else: print("")

print(r"\text{\scriptsize $\begin{matrix}")
for i, b in enumerate(basis):
    print_vector(b, final = i == len(basis) - 1)
print(r" \end{matrix}$}")
