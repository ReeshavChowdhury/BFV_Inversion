# BFV slot inversion

This program encrypts one number under the BFV homomorphic encryption
scheme and computes its multiplicative inverse, without decrypting in
the middle. It runs two circuits on the same ciphertext and checks that
both decrypt to the correct inverse.

Microsoft SEAL 4.1 must already be installed. The source tree contains
no compiled programs. The commands below create them.

## Requirements

- CMake 3.16 or newer
- A C++17 compiler
- Microsoft SEAL 4.1 (SEAL 4.0 is used if 4.1 is not installed)

`find_package(SEAL)` must succeed. SEAL is often installed under
`/usr/local`.

## Build

From this directory:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

That produces:

- `build/bfv_inversion` — one correctness test
- `build/bfv_bench` — a timing test

## Check the result

The default is the parameter set used for the main experiment:

- plaintext prime 163841
- polynomials of degree 131072
- 16384 numbers packed in one ciphertext
- ciphertext modulus made of 24 primes, 50 bits each

This takes a few minutes.

```bash
./build/bfv_inversion 12345
```

`12345` is the secret value. Any integer from 1 to 163840 is accepted.
A correct run prints the following. `FMC` is the first circuit
(Frobenius Monomial Circuit) and `NT` is the second (Norm Tower).

```text
m         = 12345
m^{-1}    = 143190
m*m^{-1}  = 1
FMC result = 143190  CORRECT
[depth   21]  FMC
m * FMC(m) mod p = 1  (=1 OK)
NT  result = 143190  CORRECT
[depth   22]  NT
```

Depth 21 and depth 22 count encrypted multiplications along the longest
path of each circuit. The program exits with status 0 only when both
answers equal 143190 and the check product equals 1.

A smaller parameter set finishes in under a minute. It uses prime
40961, degree 32768, and 4096 packed numbers. The inverse of 12345
modulo 40961 is 15243. The depths are 19 and 20.

```bash
./build/bfv_inversion 12345 small
```

## Time the two circuits

Keys are generated first and are not included in the times. The number
is how many times each circuit is repeated. The median is the middle run.

Main experiment, three repeats:

```bash
./build/bfv_bench 3
```

Smaller parameter set, three repeats:

```bash
./build/bfv_bench 3 small
```

`./build/bfv_bench` alone is the same as `./build/bfv_bench 3`.

Help text:

```bash
./build/bfv_inversion --list
```

## What is computed

One BFV ciphertext holds many independent numbers, called slots. A
single circuit operation acts on every slot. Here each slot is inverted.

Two circuits are implemented.

The Frobenius Monomial Circuit has multiplicative depth 21 on the main
parameter set (19 on the smaller one). It evaluates

\[
m^{-1} = m^{p-2} \cdot \kappa_p(m^{p-1}) \cdot \kappa_{p^2}(m^{p-1}) \cdots
\]

up to \(d-1\) applications of the Frobenius map \(\kappa_p\), which
sends a slot \(m\) to \(m^p\). Those maps are key-switches. They do not
count as multiplications. For the prime 163841,

\[
p-1 = 2^{17}+2^{15}, \qquad p-2 = 2^{17}+2^{15}-1,
\]

so both powers are products of values produced by one sequence of
squarings. In that sequence the next square and the running product do
not depend on each other, and they are computed at the same time.

The Norm Tower has depth 22 on the main parameter set (20 on the smaller
one). It sends the slot down to the base field by a product of Frobenius
images (the field norm), inverts that base-field element, and multiplies
back by the remaining images. The product of the remaining images runs
at the same time as the base-field inversion.

After every multiplication the ciphertext modulus loses one prime, so
later steps are smaller. Each thread uses its own SEAL memory pool.

Both parameter sets are inside a 128-bit security estimate for the
modulus they use.

| | Main experiment | Smaller test (`small`) |
|---|---|---|
| Plaintext prime | 163841 | 40961 |
| Polynomial degree | 131072 | 32768 |
| Numbers in one ciphertext | 16384 | 4096 |
| Ciphertext modulus | 24 primes of 50 bits | 22 primes of 32 bits |
| First circuit | depth 21 | depth 19 |
| Second circuit | depth 22 | depth 20 |
