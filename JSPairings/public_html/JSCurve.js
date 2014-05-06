/*
 * JSCurve.js
 * Barreto-Naehrig (BN) pairing-friendly elliptic curves.
 * Copyright (C) Jorge H. C. Lucema.
 */


var JSCurve = (function () {

    /*
     * Build the standard BN curve BN(u): y^2 = x^3 + b.
     * @param {BNParams} bn BN parameters of the curve
     * @returns the desired curve, or null if the given index does not define suitable parameters
     */
    var JSCurve = function (bn) {

        if (bn instanceof JSParams) {

            this.bn = bn;
            this.b = (bn.b === 3) ? bn._3 : bn._2; // standard curve
            this.infinity = new JSPoint(this); // caveat: must be set *after* p but *before* G!
            this.G = (bn.b === 3) ? new JSPoint(this, bn._1, bn._2) : new JSPoint(this, bn._1.negate(), bn._1); // standard curve
        }
    };

    /*
     * Get a random nonzero point on this curve, given a fixed base point.
     * @param {SecureRandom} rand a cryptographically strong PRNG
     * @returns a random nonzero point on this curve
     */
    JSCurve.pointFactory = function (rand) {
        if (rand instanceof SecureRandom) {
            do {
                var x = new BigInteger(2*this.bn.p.bitLength(), rand).mod(this.bn.p);
                var y = this.bn.sqrt(x.multiply(x).multiply(x).add(this.b));
            } while (y === null);
            return new BNPoint(this, x, y);
        } else {
            throw new Error("Parameter is not a cryptographically strong PRNG");
        }
    };

    /*
     * 
     */
    JSCurve.getCurveParams = function () {
        return this.bn;
    };

    /*
     * 
     */
    JSCurve.getOrder = function () {
        return this.bn.n;
    };

    /*
     * 
     */
    JSCurve.getCurveCoefficient = function () {
        return this.b;
    };

    /*
     * 
     */
    JSCurve.getCurveGenerator = function () {
        return this.G;
    };

    /*
     * Check whether this curve contains a given point
     * @param {BNPoint} P the point whose pertinence or not to this curve is to be determined
     * @returns true if this curve contains P, otherwise false
     */
    JSCurve.contains = function (P) {
        if (P.E !== this) {
            return false;
        }
        // check the projective equation y^2 = x^3 + b*z^6,
        // i.e. x*x^2 + b*z^2*(z^2)^2 - y^2 = 0
        // (the computation below never uses intermediate values larger than 3p^2)
        var x, y, z = new BigInteger();
        var x2, z2, z4, br = new BigInteger();
        x  = P.x,
        y  = P.y,
        z  = P.z,
        x2 = x.multiply(x).mod(this.bn.p),
        z2 = z.multiply(z).mod(this.bn.p),
        z4 = z2.multiply(z2).mod(this.bn.p),
        br = this.b.multiply(z2).mod(this.bn.p);
        return x.multiply(x2).add(br.multiply(z4)).subtract(y.multiply(y)).mod(this.bn.p).signum() === 0;
    };

    JSCurve.prototype.pointFactory = JSCurve.pointFactory;
    JSCurve.prototype.getCurveParams = JSCurve.getCurveParams;
    JSCurve.prototype.getOrder = JSCurve.getOrder;
    JSCurve.prototype.getCurveCoefficient = JSCurve.getCurveCoefficient;
    JSCurve.prototype.getCurveGenerator = JSCurve.getCurveGenerator;
    JSCurve.prototype.contains = JSCurve.contains;

    return JSCurve;
})();
