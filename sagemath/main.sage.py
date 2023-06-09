

# This file was *autogenerated* from the file main.sage
from sage.all_cmdline import *   # import sage library

_sage_const_1 = Integer(1); _sage_const_91 = Integer(91); _sage_const_57 = Integer(57); _sage_const_110 = Integer(110); _sage_const_67 = Integer(67); _sage_const_191 = Integer(191); _sage_const_117 = Integer(117); _sage_const_273 = Integer(273); _sage_const_172 = Integer(172); _sage_const_356 = Integer(356); _sage_const_215 = Integer(215); _sage_const_216 = Integer(216); _sage_const_137 = Integer(137); _sage_const_250 = Integer(250); _sage_const_159 = Integer(159); _sage_const_305 = Integer(305); _sage_const_192 = Integer(192); _sage_const_372 = Integer(372); _sage_const_239 = Integer(239); _sage_const_2 = Integer(2); _sage_const_3 = Integer(3); _sage_const_0 = Integer(0); _sage_const_6 = Integer(6)
import sys
import argparse

pXXX = ['p182', 'p217', 'p377', 'p546', 'p697', 'p434', 'p503', 'p610', 'p751']
# -------------------------------
def arguments(args=sys.argv[_sage_const_1 :]):
    parser = argparse.ArgumentParser(description="Parses command.")
    parser.add_argument("--prime", type=str, help="prime field characteristic,", required=True, choices=pXXX)
    parser.add_argument("--strategies", help="optimization using strategies,", action='store_true')
    parser.add_argument("--shortcut", help="shortcut optimization,", action='store_true')

    if len(sys.argv) == _sage_const_1 :
        parser.print_help(sys.stderr)
        sys.exit(_sage_const_1 )

    options = parser.parse_args(args)
    return options

import public_values_aux
from public_values_aux import *

# SIKEpXXX parameters
pxxx = arguments(sys.argv[_sage_const_1 :]).prime
strategies = arguments(sys.argv[_sage_const_1 :]).strategies
shortcut = arguments(sys.argv[_sage_const_1 :]).shortcut

if pxxx == 'p182':
    a, b = _sage_const_91 , _sage_const_57 
elif pxxx == 'p217':
    a, b = _sage_const_110 , _sage_const_67 
elif pxxx == 'p377':
    a, b = _sage_const_191 , _sage_const_117 
elif pxxx == 'p546':
    a, b = _sage_const_273 , _sage_const_172 
elif pxxx == 'p697':
    a, b = _sage_const_356 , _sage_const_215 
elif pxxx == 'p434':
    a, b = _sage_const_216 , _sage_const_137 
elif pxxx == 'p503':
    a, b = _sage_const_250 , _sage_const_159 
elif pxxx == 'p610':
    a, b = _sage_const_305 , _sage_const_192 
elif pxxx == 'p751':
    a, b = _sage_const_372 , _sage_const_239 
else:
    exit(-_sage_const_1 )

if not shortcut:
    load('castryck_decru_attack.sage')
else:
    load('castryck_decru_shortcut.sage')


# Set the prime, finite fields and starting curve
# with known endomorphism
p = _sage_const_2 **a*_sage_const_3 **b - _sage_const_1 
public_values_aux.p = p

Fp2 = GF(p**_sage_const_2 , modulus=x**_sage_const_2 +_sage_const_1 , names=('i',)); (i,) = Fp2._first_ngens(1)
R = PolynomialRing(Fp2, names=('x',)); (x,) = R._first_ngens(1)

E_start = EllipticCurve(Fp2, [_sage_const_0 ,_sage_const_6 ,_sage_const_0 ,_sage_const_1 ,_sage_const_0 ])
E_start.set_order((p+_sage_const_1 )**_sage_const_2 ) # Speeds things up in Sage

# Generation of the endomorphism 2i
two_i = generate_distortion_map(E_start)

# Generate public torsion points, for SIKE implementations
# these are fixed but to save loading in constants we can
# just generate them on the fly
P2, Q2, P3, Q3 = generate_torsion_points(E_start, a, b)
check_torsion_points(E_start, a, b, P2, Q2, P3, Q3)

# Generate Bob's key pair
bob_private_key, EB, PB, QB = gen_bob_keypair(E_start, b, P2, Q2, P3, Q3)
solution = Integer(bob_private_key).digits(base=_sage_const_3 )

print(f"Running the attack against SIDHp{p.bit_length()} parameters, which has a prime: 2^{a}*3^{b} - 1")
print(f"If all goes well then the following digits should be found: {solution}")


def RunAttack(num_cores):
        return CastryckDecruAttack(E_start, P2, Q2, EB, PB, QB, two_i, num_cores=num_cores, strategies=strategies)


if __name__ == '__main__' and '__file__' in globals():
    if '--parallel' in sys.argv:
        # Set number of cores for parallel computation
        num_cores = os.cpu_count()
        print(f"Performing the attack in parallel using {num_cores} cores")
    else:
        num_cores = _sage_const_1 
    recovered_key = RunAttack(num_cores)


