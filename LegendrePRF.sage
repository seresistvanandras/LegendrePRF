# Draft Arithmetization of the LegendrePRF
from sage.rings.polynomial.hilbert import first_hilbert_series
from sage.rings.polynomial.hilbert import hilbert_poincare_series
prime = 0x0000000000000000000000000000000000000000000fffffffffffffffffffdd
key = 0x0000000000000000000000000000000000000000000027aaa97c746c22e12d04

legendreSequence = ''
prngByte = 0
cnt = 0

for x in range(2048):
    prngByte+=((kronecker(key+x,prime)+1)/2)*(2^(3-cnt))
    cnt+=1
    if cnt==4:
        legendreSequence+=(str(hex(prngByte))[2:])
        prngByte=0
        cnt=0

# The implementation is checked by the given challenge bits here: https://legendreprf.org/bountyinstances (See Challenge 2)
correctPrngOutput = 'bafca94ade9b5201633be31512efcaec7cbe64cbfd2806e83ca398ee34209e01a14bc727418baa31692ebc91681018527738fea54c9f1d45233fff8de5cce971b2111e012374f10ee3fbca4e276313ba8ffed1f400f1d5e046ffad63b6f48caecb7668b263190d1b0d822397b1fd72cf6a5c24f80af7254240bb432a6bb518588950e82b07e63980bbdd754ce80b39090ba04c52e4e186f42f75e7f9bd097fdf23105f123a7b95101dd053e66d84a2ddbc939815986ca510e29f2864df6a513f143800f79cc62bfed7d4b2ba0a128090ac2e7a2b4857bb703cd425f941d8e47c80ef243a7705f05beef1c4a0d2dabb1cbc22bca06bd935ca25a7237ffee5bfcc'
print("Our own LegendrePRF implementation is correctly implemented",legendreSequence==correctPrngOutput)

# We need a single quadratic non-residue (denoted here by r) for the arithmetization of the given Legendre sequence
r=2
while kronecker(r,prime)!=-1:
    r+=1

print("Smallest quadratic non-residue mod",prime,"is",r)


n = 20 #length of the Legendre PRF output we want to arithmetize
lIdeal = []
R = PolynomialRing(GF(prime),n,'x'); R
monomial = [0]*n #constant term in the multivariate polynomial ring
term1 = []
term2 = []
for x in range(n-1):
    if kronecker(key+x,prime)==1 and kronecker(key+x+1,prime)==1:
        term1 = [0]*n
        term2 = [0]*n
        term1[x+1] = 2
        term2[x] = 2
        quadreq = R({tuple(term1):1,tuple(term2):-1,tuple(monomial):-1})
        lIdeal.append(quadreq)
        print('11',quadreq)
    elif kronecker(key+x,prime)==1 and kronecker(key+x+1,prime)==-1:
        term1 = [0]*n
        term2 = [0]*n
        term1[x+1] = 2
        term2[x] = 2
        quadreq = R({tuple(term1):1,tuple(term2):-r,tuple(monomial):-r})
        lIdeal.append(quadreq)
        print('10',quadreq)
    elif kronecker(key+x,prime)==-1 and kronecker(key+x+1,prime)==1:
        term1 = [0]*n
        term2 = [0]*n
        term1[x+1] = 2
        term2[x] = 2
        quadreq = R({tuple(term1):r,tuple(term2):-1,tuple(monomial):-r})
        lIdeal.append(quadreq)
        print('01',quadreq)
    else:
        term1 = [0]*n
        term2 = [0]*n
        term1[x+1] = 2
        term2[x] = 2
        quadreq = R({tuple(term1):1,tuple(term2):-1,tuple(monomial):-r})
        lIdeal.append(quadreq)
        print('00',quadreq)
print(R)
myIdeal = Ideal(lIdeal)
gb = myIdeal.groebner_basis()
print(gb)
print("Maximum degree in the Gröbner basis",max(f.degree() for f in gb))
print("Hilbert polynomial of the ideal",myIdeal.homogenize().hilbert_polynomial())
print("We have got back a reduced Gröbner-basis:", gb == gb.reduced())
print("Hilbert series of the ideal",first_hilbert_series(myIdeal))
print("Degree of regularity:",hilbert_poincare_series(myIdeal).numerator().degree())
