# Computing Quotient Groups of Smooth Order with Applications to Isogenies over Higher-Dimensional Abelian Varieties

This code corresponds with the work titled: "Computing Quotient Groups of Smooth Order with Applications to Isogenies over Higher-Dimensional Abelian Varieties".
The preprint version is available at [eprint 2023/508](https://eprint.iacr.org/2023/508).


## Requirements

The **Magma** Computer Algebra System, and the **SageMath** library.

## Description

The **SageMath**-language code takes as a baseline the public code from [eprint 2022/1283](https://eprint.iacr.org/2022/1283).

The **Magma**language code takes as a baseline the following two public codes:

The (2^n,2^n)-isogeny implementation from [eprint 2022/990](https://eprint.iacr.org/2022/990), and
The (3^n,3^n)-isogeny implementation from [eprint 2023/376](https://eprint.iacr.org/2023/376).

## Examples

The **SageMath**-language code allows testing different over prime fields and allows to decide using [or not] the strategies.
We automatize using the __shortcut__ as an input argument.
Adding more prime fields from the file `main.sage` is straightforward and intuitive.
```bash
cd sagemath
% sage main.sage --help
usage: main.sage.py [-h] --prime {p182,p217,p377,p546,p697,p434,p503,p610,p751} [--strategies] [--shortcut]

Parses command.

options:
 -h, --help      show this help message and exit
 --prime {p182,p217,p377,p546,p697,p434,p503,p610,p751}
            prime field characteristic,
 --strategies     optimization using strategies,
 --shortcut      shortcut optimization,

```


### Without strategies

```bash
% sage main.sage --prime p182 --shortcut 
Running the attack against SIDHp182 parameters, which has a prime: 2^91*3^57 - 1
If all goes well then the following digits should be found: [1, 1, 0, 2, 2, 0, 0, 2, 1, 2, 2, 1, 2, 0, 2, 0, 1, 1, 2, 0, 0, 0, 1, 2, 1, 1, 1, 1, 2, 0, 0, 0, 1, 1, 2, 1, 0, 0, 1, 2, 1, 1, 0, 2, 0, 1, 2, 0, 1, 0, 2, 0, 2, 1, 1, 1]
Determination of first 2 ternary digits. We are working with 2^91-torsion.
Strategy technique:	False

Testing digits: [0, 0]
Testing digits: [1, 0]
Testing digits: [2, 0]
Testing digits: [0, 1]
Testing digits: [1, 1]
Computing image of 3-adic torsion in split factor CB
Glue-and-split! These are most likely the secret digits.
Bob's secret key revealed as: 266441314757113948437860404
In ternary, this is: [1, 1, 0, 2, 2, 0, 0, 2, 1, 2, 2, 1, 2, 0, 2, 0, 1, 1, 2, 0, 0, 0, 1, 2, 1, 1, 1, 1, 2, 0, 0, 0, 1, 1, 2, 1, 0, 0, 1, 2, 1, 1, 0, 2, 0, 1, 2, 0, 1, 0, 2, 0, 2, 1, 1, 1]
Altogether this took 7.650943994522095 seconds.
```

### With strategies

```bash
% sage main.sage --prime p182 --strategies --shortcut
Running the attack against SIDHp182 parameters, which has a prime: 2^91*3^57 - 1
If all goes well then the following digits should be found: [1, 1, 1, 1, 0, 0, 0, 1, 2, 2, 2, 1, 1, 0, 2, 0, 0, 2, 2, 0, 2, 2, 0, 2, 0, 2, 1, 0, 2, 1, 2, 0, 2, 2, 1, 0, 2, 2, 2, 1, 1, 2, 0, 1, 0, 0, 2, 1, 2, 0, 0, 2, 1, 1, 1, 0, 1]
Determination of first 2 ternary digits. We are working with 2^91-torsion.
Strategy technique:	True

Testing digits: [0, 0]
Testing digits: [1, 0]
Testing digits: [2, 0]
Testing digits: [0, 1]
Testing digits: [1, 1]
Computing image of 3-adic torsion in split factor CB
Glue-and-split! These are most likely the secret digits.
Bob's secret key revealed as: 611853354437883654938185282
In ternary, this is: [1, 1, 1, 1, 0, 0, 0, 1, 2, 2, 2, 1, 1, 0, 2, 0, 0, 2, 2, 0, 2, 2, 0, 2, 0, 2, 1, 0, 2, 1, 2, 0, 2, 2, 1, 0, 2, 2, 2, 1, 1, 2, 0, 1, 0, 0, 2, 1, 2, 0, 0, 2, 1, 1, 1, 0, 1]
Altogether this took 7.04803204536438 seconds.
```

### The **Magma** code

Just run
```bash
cd magma/22_isogeny
magma 22_test_p171.m

cd magma/33_isogeny
magma SIKEp751_attack.m
```


## Remarks

We already did a [Pull Request](https://github.com/jack4818/Castryck-Decru-SageMath/pull/27) to the repository from [eprint 2022/1283](https://eprint.iacr.org/2022/1283).
