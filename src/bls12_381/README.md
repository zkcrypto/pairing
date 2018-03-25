# BLS12-381

This is an implementation of the BLS12-381 pairing-friendly elliptic curve construction.

## BLS12 Parameterization

BLS12 curves are parameterized by a value *x* such that the base field modulus *q* and subgroup *r* can be computed by:

* q = (x - 1)<sup>2</sup> ((x<sup>4</sup> - x<sup>2</sup> + 1) / 3) + x
* r = (x<sup>4</sup> - x<sup>2</sup> + 1)

Given primes *q* and *r* parameterized as above, we can easily construct an elliptic curve over the prime field F<sub>*q*</sub> which contains a subgroup of order *r* such that *r* | (*q*<sup>12</sup> - 1), giving it an embedding degree of 12. Instantiating its sextic twist over an extension field F<sub>q<sup>2</sup></sub> gives rise to an efficient bilinear pairing function between elements of the order *r* subgroups of either curves, into an order *r* multiplicative subgroup of F<sub>q<sup>12</sup></sub>.

In zk-SNARK schemes, we require F<sub>r</sub> with large 2<sup>n</sup> roots of unity for performing efficient fast-fourier transforms. As such, guaranteeing that large 2<sup>n</sup> | (r - 1), or equivalently that *x* has a large 2<sup>n</sup> factor, gives rise to BLS12 curves suitable for zk-SNARKs.

Due to recent research, it is estimated by many that *q* should be approximately 384 bits to target 128-bit security. Conveniently, *r* is approximately 256 bits when *q* is approximately 384 bits, making BLS12 curves ideal for 128-bit security. It also makes them ideal for many zk-SNARK applications, as the scalar field can be used for keying material such as embedded curve constructions.

Many curves match our descriptions, but we require some extra properties for efficiency purposes:

* *q* should be smaller than 2<sup>383</sup>, and *r* should be smaller than 2<sup>255</sup>, so that the most significant bit is unset when using 64-bit or 32-bit limbs. This allows for cheap reductions.
* F<sub>q<sup>12</sup></sub> is typically constructed using towers of extension fields. As a byproduct of [research](https://eprint.iacr.org/2011/465.pdf) for BLS curves of embedding degree 24, we can identify subfamilies of BLS12 curves (for our purposes, where x mod 72 = {16, 64}) that produce efficient extension field towers and twisting isomorphisms.
* We desire *x* of small Hamming weight, to increase the performance of the pairing function.

## BLS12-381 Instantiation

The BLS12-381 construction is instantiated by `x = -0xd201000000010000`, which produces the largest `q` and smallest Hamming weight of `x` that meets the above requirements. This produces:

* q = `0x1a0111ea397fe69a4b1ba7b6434bacd764774b84f38512bf6730d2a0f6b0f6241eabfffeb153ffffb9feffffffffaaab` (381 bits)
* r = `0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001` (255 bits)

Our extension field tower is constructed as follows:

1. F<sub>q<sup>2</sup></sub> is constructed as F<sub>q</sub>(u) / (u<sup>2</sup> - β) where β = -1.
2. F<sub>q<sup>6</sup></sub> is constructed as F<sub>q<sup>2</sup></sub>(v) / (v<sup>3</sup> - ξ) where ξ = u + 1
3. F<sub>q<sup>12</sup></sub> is constructed as F<sub>q<sup>6</sup></sub>(w) / (w<sup>2</sup> - γ) where γ = v

Now, we instantiate the elliptic curve E(F<sub>q</sub>) : y<sup>2</sup> = x<sup>3</sup> + 4, and the elliptic curve E'(F<sub>q<sup>2</sup></sub>) : y<sup>2</sup> = x<sup>3</sup> + 4(u + 1).

The group G<sub>1</sub> is the *r* order subgroup of E, which has cofactor (x - 1)<sup>2</sup> / 3. The group G<sub>2</sub> is the *r* order subgroup of E', which has cofactor (x<sup>8</sup> - 4x<sup>7</sup> + 5x<sup>6</sup> - 4x<sup>4</sup> + 6x<sup>3</sup> - 4x<sup>2</sup> - 4x + 13) / 9.

### Generators

The generators of G<sub>1</sub> and G<sub>2</sub> are computed by finding the lexicographically smallest valid `x`-coordinate, and its lexicographically smallest `y`-coordinate and scaling it by the cofactor such that the result is not the point at infinity.

#### G1

```
x = 3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507
y = 1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569
```

#### G2

```
x = 3059144344244213709971259814753781636986470325476647558659373206291635324768958432433509563104347017837885763365758*u + 352701069587466618187139116011060144890029952792775240219908644239793785735715026873347600343865175952761926303160
y = 927553665492332455747201965776037880757740193453592970025027978793976877002675564980949289727957565575433344219582*u + 1985150602287291935568054521177171638300868978215655730859378665066344726373823718423869104263333984641494340347905
```

### Serialization

* Fq elements are encoded in big-endian form. They occupy 48 bytes in this form.
* Fq2 elements are encoded in big-endian form, meaning that the Fq element c0 + c1 * u is represented by the Fq element c1 followed by the Fq element c0. This means Fq2 elements occupy 96 bytes in this form.
* The group G1 uses Fq elements for coordinates. The group G2 uses Fq2 elements for coordinates.
* G1 and G2 elements can be encoded in uncompressed form (the x-coordinate followed by the y-coordinate) or in compressed form (just the x-coordinate). G1 elements occupy 96 bytes in uncompressed form, and 48 bytes in compressed form. G2 elements occupy 192 bytes in uncompressed form, and 96 bytes in compressed form.

The most-significant three bits of a G1 or G2 encoding should be masked away before the coordinate(s) are interpreted. These bits are used to unambiguously represent the underlying element:

* The most significant bit, when set, indicates that the point is in compressed form. Otherwise, the point is in uncompressed form.
* The second-most significant bit indicates that the point is at infinity. If this bit is set, the remaining bits of the group element's encoding should be set to zero.
* The third-most significant bit is set if (and only if) this point is in compressed form _and_ it is not the point at infinity _and_ its y-coordinate is the lexicographically largest of the two associated with the encoded x-coordinate.

### Hashing to the Curve

* Elements of Fq are encoded by taking the output of a BLAKE2b digest and interpreting it as a big-endian integer. The integer is reduced (mod q), making it uniform in the field with negligible bias.
* Elements of Fq2 are encoded by appending "_c0" to the supplied BLAKE2b preimage, to compute the encoding to Fq for c0, and by appending "_c1" to the same BLAKE2b preimage, to compute the encoding to Fq for c1.
* Elements of E and E' are encoded from an element *t* by taking the first valid abscissa (for the b in each respective curve):
    * x_1 = (-1 + sqrt(-3))/2 - (sqrt(-3) * t^2)/(1 + b + t^2)
    * x_2 = (-1 - sqrt(-3))/2 + (sqrt(-3) * t^2)/(1 + b + t^2)
    * x_3 = 1 - (1 + b + t^2)^2 / (3 * t^2)
* In this encoding, we always map t=0 to the point at infinity. For the encoding to E, We also map:
    * `t = 0x019cfaba0c258165d092f6bca9a081871e62a126c499340dc71c0e9527f923f3b299592a7a9503066cc5362484d96dd7` to the fixed generator. (See "Generators" above.)
    * `t = 0x186417302d5a65347a88b0f999ab2b504614aa5e2eebdeb1a014c40bceb7d2306c12a6d436befcf94d39c9db7b263cd4` to the negative of the fixed generator. (See "Generators" above.)
* In this encoding, the y-coordinate of the resulting point is chosen to be the lexicographically largest if (and only if) the input `t` is also lexicographically larger than its negative.
* Hashing to G1, given a message `msg`, involves hashing to Fq by supplying the BLAKE2b preimage `msg | "G1_0"` to create `t0`, and hashing to Fq by supplying the BLAKE2b preimage `msg | "G1_1"` to create `t1`. `t0` and `t1` are used to encode into E via the above encoding, the resulting points are added together, and that result is multiplied by the cofactor.
* Hashing to G2, given a message `msg`, involves hashing to Fq2 by supplying the BLAKE2b preimage `msg | "G2_0"` to create `t0`, and hashing to Fq2 by supplying the BLAKE2b preimage `msg | "G2_1"` to create `t1`. `t0` and `t1` are used to encode into E' via the above encoding, the resulting points are added together, and that result is multiplied by the cofactor.
