const JSCurve = require('./JSCurve');
const JSField2 = require('./JSField2');
const BigInteger = require('./BigNumber').BigInteger;

/*
 * JSCurve2.js
 * Sextic twist of Barreto-Naehrig (BN) pairing-friendly elliptic curves.
 * Copyright (C) Jorge H. C. Lucema.
 */


/*
    * Build the standard sextic twist E': y'^2 = x'^3 + b' of the standard BN curve E: y^2 = x^3 + b.
    * @param {JSCurve} E given BN curve.
    * @returns the desired curve, or null if E does not have a suitable twist of the above form.
    */
var JSCurve2 = function (E) {

    if (E instanceof JSCurve) {
        this.E = E;
        this.Fp2_0 = E.bn.Fp2_0;
        this.Fp2_1 = E.bn.Fp2_1;
        this.Fp2_i = E.bn.Fp2_i;
        const JSPoint2 = require('./JSPoint2');
        this.infinity = new JSPoint2(this);
        if (E.b.intValue() === 3) {
            this.bt = new JSField2(E.bn, E.b).multiplyV(); // b' = b*(1 + i), standard non-square non-cube
            this.xt = this.Fp2_1; // standard x-coord
            this.yt = this.xt.multiply(this.xt).multiply(this.xt).add(this.bt).sqrt();
        } else {
            this.bt = this.Fp2_1.subtract(this.Fp2_i); // b' = 1 + i, standard non-square non-cube
            this.xt = this.Fp2_i.negate();
            this.yt = this.Fp2_1;
        }
        this.Gt = new JSPoint2(this, this.xt, this.yt);
        this.Gt = this.Gt.multiply(E.bn.ht).normalize();
        this.pp16Gt = new Array(Math.round((this.E.bn.n.bitLength() + 3)/4));
        for (var i = 0; i < Math.round((this.E.bn.n.bitLength() + 3)/4); i++) {
            this.pp16Gt[i] = new Array(16);
        }
        this.pp16Gi = this.pp16Gt[0];
        this.pp16Gi[0] = this.infinity;
        this.pp16Gi[1] = this.Gt;
        for (var i = 1, j = 2; i <= 7; i++, j += 2) {
            this.pp16Gi[j] = this.pp16Gi[i].twice(1);
            this.pp16Gi[j+1] = this.pp16Gi[j].add(this.Gt);
        }
        for (var i = 1; i < this.pp16Gt.length; i++) {
            this.pp16Gh = this.pp16Gi;
            this.pp16Gi = this.pp16Gt[i];
            this.pp16Gi[0] = this.pp16Gh[0];
            for (var j = 1; j < 16; j++) {
                this.pp16Gi[j] = this.pp16Gh[j].twice(4);
            }
        }
    }
};

/*
    * Get a random nonzero point on this curve, given a fixed base point.
    * @param {SecureRandom} rand a cryptographically strong PRNG
    * @return a random nonzero point on this curve
    */
JSCurve2.pointFactory = function (rand) {
    do {
        var k = new BigInteger(this.E.bn.n.bitLength(), rand).mod(this.E.bn.n);
    } while (k.signum() === 0);
    return this.Gt.multiply(k);
};

/*
    * Check whether this curve contains a given point
    * (i.e. whether that point satisfies the curve equation)
    * @param {JSPoint2} P the point whose pertinence or not to this curve is to be determined
    * @return true if this curve contains P, otherwise false
    */
JSCurve2.contains = function (P) {
    if (P.E !== this) {
        return false;
    }
    // check the projective equation y^2 = x^3 + bt*z^6,
    // i.e. x^3 + bt*(z^2)^3 - y^2 = 0
    var x  = P.x, y  = P.y, z  = P.z;
    return y.square().equals(x.cube().add(this.bt.multiply(z.square().cube())));
};

/*
    * Compute k*G
    * @param {BigInteger} k scalar by which the base point G is to be multiplied.
    * @return  k*G
    */
JSCurve2.kG = function (k) {
    k = k.mod(this.E.bn.n); // reduce k mod n
    var A = this.infinity;
    for (var i = 0, w = 0; i < this.pp16Gt.length; i++, w >>>= 4) {
        if ((i & 7) === 0) {
            w = k.intValue();
            k = k.shiftRight(32);
        }
        A = A.add(this.pp16Gt[i][w & 0xf]);
    }
    return A;
};

JSCurve2.prototype.pointFactory = JSCurve2.pointFactory;
JSCurve2.prototype.contains = JSCurve2.contains;
JSCurve2.prototype.kG = JSCurve2.kG;

module.exports = JSCurve2;