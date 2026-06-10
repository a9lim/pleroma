# Brute-force squares in Z/2^k vs is_square_mod_two_power
def is_square_mod_two_power(a, k):
    if a == 0: return True
    u, v = a, 0
    while u % 2 == 0:
        u //= 2; v += 1
    if v % 2 != 0: return False
    m = k - v
    if m in (0, 1): return True
    if m == 2: return u % 4 == 1
    return u % 8 == 1

for k in range(1, 12):
    M = 2**k
    squares = set((x*x) % M for x in range(M))
    for a in range(M):
        assert is_square_mod_two_power(a, k) == (a in squares), (a, k)
print("is_square_mod_two_power matches brute force for k=1..11")

# Qp capped-relative addition spot check: v(2+6)=3 in Z_2 (2+6=8)
# code: lo=2=2^1*1, hi=6=2^1*3, d=0, b=(1+3)=4 -> normalized(4, 1) -> unit 1, val 3. True v_2(8)=3 ok.
print("v_2(2+6): expected 3, code gives 1 +", (4).bit_length()-1, "= 3")

# 1/3 in Z_5 to 4 digits: 3^{-1} mod 625 = ?
inv3 = pow(3, -1, 625)
print("1/3 mod 5^4 =", inv3, "; 3*", inv3, "mod 625 =", (3*inv3)%625)
# 1+2+4+8+... = -1 in Z_2 (sanity for the model: from_i128(-1) unit = 2^k - 1)
print("-1 mod 2^6 =", (-1) % 64)
