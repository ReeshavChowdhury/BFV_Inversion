# BFV slot inversion

This program encrypts a number with the BFV encryption scheme and computes
its inverse without decrypting in the middle. It runs the two methods from
the paper on the same encrypted value and checks that both decrypt to the
right answer.

The source folder does not contain compiled programs. The commands below
build them. Microsoft SEAL 4.1 must already be installed.

## Requirements

- CMake 3.16 or newer
- A C++17 compiler
- Microsoft SEAL 4.1. If SEAL 4.1 is missing, the build uses SEAL 4.0.

SEAL is usually installed under `/usr/local`, so that CMake can find it.

## Build

From this directory:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

This creates two programs:

- `build/bfv_inversion` checks that the answer is correct
- `build/bfv_bench` measures how long each method takes

## Check the answer

The default settings are the ones used in the paper:

- plaintext prime 163841
- polynomials of degree 131072
- 16384 numbers packed into one ciphertext
- ciphertext modulus made from 24 primes of 50 bits each

The run takes a few minutes.

```bash
./build/bfv_inversion 12345
```

Here 12345 is the secret number. Any whole number from 1 to 163840 is
allowed. A correct run prints:

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

FMC is the first method in the paper, the Frobenius Monomial Circuit.
NT is the second method, the Norm Tower. The depth is how many encrypted
multiplications sit on the longest path. The direct way of computing the
inverse, by raising the number to a huge power, needs depth 139 on these
parameters. The Norm Tower needs depth 22. The Frobenius Monomial Circuit
needs depth 21, which the paper shows is the smallest possible depth.

The program exits with status 0 only when both answers are 143190 and the
check product is 1.

A smaller test finishes in under a minute. It uses prime 40961, degree
32768, and 4096 packed numbers. The inverse of 12345 is then 15243, and
the two depths are 19 and 20.

```bash
./build/bfv_inversion 12345 small
```

## Measure the running time

Key generation is finished before the timer starts. The number after the
program name is how many times each method is repeated. The median is the
middle of those repeats.

Paper parameters, three repeats:

```bash
./build/bfv_bench 3
```

Smaller parameters, three repeats:

```bash
./build/bfv_bench 3 small
```

Running `./build/bfv_bench` with no arguments is the same as
`./build/bfv_bench 3`.

On the machine used for the paper's timing table (a 10-core Intel Core
i9-10900, Microsoft SEAL 4.1.2), three repeats of the main experiment gave:

| Method | Run 1 | Run 2 | Run 3 | Middle value |
|---|---:|---:|---:|---:|
| Frobenius Monomial Circuit, first stage | 31.0 s | 30.4 s | 30.5 s | 30.5 s |
| Frobenius Monomial Circuit, second stage | 1.0 s | 1.0 s | 1.0 s | 1.0 s |
| Frobenius Monomial Circuit, total, depth 21 | 32.0 s | 31.3 s | 31.5 s | 31.5 s |
| Norm Tower, total, depth 22 | 45.8 s | 45.6 s | 45.9 s | 45.8 s |

The first method inverts all 16384 packed numbers in 31.5 seconds, about
521 numbers per second. The Norm Tower takes 45.8 seconds, about 358
numbers per second. The first method is about 31 percent faster. Almost
all of its time, 30.5 of the 31.5 seconds, is the first stage. The second
stage, which multiplies the eight pieces together, takes 1.0 second. The
Norm Tower is slower because it first multiplies eight full-size
ciphertexts before it can start the same kind of power computation.

On the smaller parameter set the same three-repeat measurement gave 5.5
seconds for the Frobenius Monomial Circuit and 7.9 seconds for the Norm
Tower.

```bash
./build/bfv_inversion --list
```

prints these two choices again.

## What the two methods do

One ciphertext holds many numbers at once, and one operation changes all
of them. Both methods turn every packed number into its inverse.

The encryption scheme can raise a packed number to the p-th power, and
repeat that, at no extra multiplication depth. Those images are the raw
material for both methods. Only one extra key is needed for this map.
Later images are made by applying the same key again.

The Frobenius Monomial Circuit builds the inverse from those images with
21 encrypted multiplications on the longest path. For the prime 163841
the two hard powers are produced from one shared sequence of squarings.
The next square and the running product do not depend on each other, so
they are computed at the same time. The eight results are then multiplied
together in a tree of depth 3.

The Norm Tower first multiplies all the images together. That product
lands in the smaller base field, where inversion is cheaper, and a final
multiplication puts the result back in the original slot. That final
multiplication is why this method has depth 22, one more than the other.
The extra product that rebuilds the answer runs at the same time as the
base-field inversion.

After every multiplication the ciphertext drops one prime from its
modulus, so later steps handle a smaller object. The main experiment uses
24 primes because depth 22, plus one check multiplication, needs 23
working levels, and one extra prime is reserved for key switching. Each
worker thread uses its own scratch memory.

Both parameter sets stay inside a 128-bit security estimate for the
modulus they use.

| | Main experiment | Smaller test |
|---|---|---|
| Plaintext prime | 163841 | 40961 |
| Polynomial degree | 131072 | 32768 |
| Numbers in one ciphertext | 16384 | 4096 |
| Ciphertext modulus | 24 primes of 50 bits | 22 primes of 32 bits |
| Frobenius Monomial Circuit | depth 21, 31.5 seconds | depth 19, 5.5 seconds |
| Norm Tower | depth 22, 45.8 seconds | depth 20, 7.9 seconds |
