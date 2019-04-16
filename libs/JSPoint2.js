const JSCurve2 = require('./JSCurve2');
const JSField2 = require('./JSField2');

/*
 * JSPoint2.js
 * Arithmetic in the group of points on the sextic twist a BN elliptic curve over GF(p^2).
 * Copyright (C) Jorge H. C. Lucema.
 */

var pointNotOnCurve = "The given point does not belong to the given elliptic curve";

var JSPoint2 = function (E, x, y, z) {

    if(arguments.length === 1) {
        /*
            * Create an instance of the JSCurve point at infinity on curve E.
            * @param {JSCurve2} E the elliptic curve where the created point is located.
            */
        if (E instanceof JSCurve2) {
            this.E = E;
            /*
                * the point at infinity is represented as (1, 1, 0) after IEEE Std 1363:2000
                * (notice that this triple satisfies the projective curve equation y^2 = x^3 + b.z^6)
                */
            this.x = E.Fp2_1;
            this.y = E.Fp2_1;
            this.z = E.Fp2_0;
        }
        /*
            * Create a clone of a given point.
            * @param {JSPoint2} Q the point to be cloned.
            */
        if (E instanceof JSPoint2) {
            var Q = E;
            this.E = Q.E;
            this.x = Q.x;
            this.y = Q.y;
            this.z = Q.z;
        }
    }
    if(arguments.length === 3) {
        /*
            * Create a normalized twist point from given affine coordinates and a curve
            * @param {JSCurve2} E the underlying elliptic curve
            * @param {JSField2} x the affine x-coordinate.
            * @param {JSField2} y the affine y-coordinate.
            */
        if ((x instanceof JSField2) && (y instanceof JSField2)) {
            this.E = E;
            this.x = x;
            this.y = y;
            this.z = E.Fp2_1; // normalized
            if (!E.contains(this)) {
                throw new Error(pointNotOnCurve);
            }
        }
        /*
            * Create an JSCurve point from a given affine x-coordinate, a y-bit and a curve
            * @param {JSCurve2} E the underlying elliptic curve.
            * @param {JSField2} x the affine x-coordinate.
            * @param {Integer} y the least significant bit of the y-coordinate.
            */
        else if (!(y instanceof JSField2)) {
            var yBit = y;
            this.E = E;
            this.x = x;
            if (this.x.isZero()) {
                throw new Error(pointNotOnCurve);
            } else {
                this.y = this.x.cube().add(this.E.bt).sqrt();
                if (this.y === null) {
                    throw new Error(pointNotOnCurve);
                }
                if (this.y.re.testBit(0) !== ((yBit & 1) === 1)) {
                    this.y = this.y.negate();
                }
            }
            this.z = this.E.Fp2_1; // normalized
        }
        /*
            * Create an JSCurve point from a given x-trit, an affine y-coordinate, and a curve
            * @param {JSCurve2} E the underlying elliptic curve.
            * @param {Integer} x the least significant trit of the x-coordinate.
            * @param {JSField2} y the affine y-coordinate.
            */
        else if (!(x instanceof JSField2)) {
            var xTrit = x;
            this.E = E;
            this.y = y;
            if (this.y.isZero()) {
                throw new Error(pointNotOnCurve); // otherwise the curve order would not be prime
            } else {
                this.x = this.y.square().subtract(this.E.bt).cbrt();
                if (this.x === null) {
                    throw new Error(pointNotOnCurve);
                }
                // either x, zeta*x, or zeta^2*x is the desired x-coordinate:
                if (this.x.re.mod(this.E.E.bn._3).intValue() !== xTrit) {
                    var zeta = this.E.E.bn.zeta; // shorthand
                    this.x = this.x.multiply(zeta);
                    if (this.x.re.mod(this.E.E.bn._3).intValue() !== xTrit) {
                        this.x = this.x.multiply(zeta);
                        if (this.x.re.mod(this.E.E.bn._3).intValue() !== xTrit) {
                            throw new Error(pointNotOnCurve);
                        }
                    }
                }
            }
            this.z = this.E.Fp2_1; // normalized
        }
    }
    /*
        * Create an JSCurve point from given projective coordinates and a curve.
        * @param {} E the underlying elliptic curve.
        * @param {} x the affine x-coordinate.
        * @param {} y the affine y-coordinate.
        * @param {} z the affine z-coordinate.
        */
    if (arguments.length === 4) {
        this.E = E;
        this.x = x;
        this.y = y;
        this.z = z;
    }
};

/*
    * performing arithmetic operations on elliptic curve points
    * generally implies knowing the nature of these points (more precisely,
    * the nature of the finite field to which their coordinates belong),
    * hence they are done by the underlying elliptic curve.
    */

/*
    * Check whether this is the point at infinity (i.e. the JSCurve group zero element).
    * @return true if this is the point at infinity, otherwise false.
    */
JSPoint2.isZero = function () {
    return this.z.isZero();
};

/*
    * Compare this point to a given object.
    * @param {Object} Q the elliptic curve point to be compared to this.
    * @returns true if this point and Q are equal, otherwise false.
    */
JSPoint2.equals = function (Q) {
    if (!(Q instanceof JSPoint2) || !this.isOnSameCurve(Q)) {
        return false;
    }
    var P = Q;
    if (this.z.isZero() || P.isZero()) {
        return this.z.equals(P.z);
    } else {
        // x/z^2 = x'/z'^2 <=> x*z'^2 = x'*z^2.
        // y/z^3 = y'/z'^3 <=> y*z'^3 = y'*z^3,
        var z2 = this.z.square();
        var z3 = this.z.multiply(z2);
        var pz2 = P.z.square();
        var pz3 = P.z.multiply(pz2);
        return this.x.multiply(pz2).equals(P.x.multiply(z2)) &&
            this.y.multiply(pz3).equals(P.y.multiply(z3));
    }
};

/*
    * Check whether Q lays on the same curve as this point.
    * @param {JSPoint2} Q an elliptic curve point.
    * @returns true if Q lays on the same curve as this point, otherwise false.
    */
JSPoint2.isOnSameCurve = function (Q) {
    return this.E.E.bn === Q.E.E.bn; // singleton comparison
};

/*
    * Compute a random point on the same curve as this.
    * @param {SecureRandom} rand a cryptographically strong pseudo-random number generator.
    * @returns a random point on the same curve as this.
    */
JSPoint2.randomize = function (rand) {
    return this.E.pointFactory(rand);
};

/*
    * Normalize this point.
    * @returns a normalized point equivalent to this.
    */
JSPoint2.normalize = function () {
    if (this.z.isZero() || this.z.isOne()) {
        return this; // already normalized
    }
    var zinv = this.z.inverse(), zinv2 = zinv.square(), zinv3 = zinv.multiply(zinv2);
    return new JSPoint2(this.E, this.x.multiply(zinv2), this.y.multiply(zinv3), this.E.Fp2_1);
};

/*
    * Compute -this.
    * @returns -this.
    */
JSPoint2.negate = function () {
    return new JSPoint2(this.E, this.x, this.y.negate(), this.z);
};

/*
    * Check if a point equals -this.
    * @param {JSPoint2} P
    * @returns true if P=-this
    */
JSPoint2.opposite = function (P) {
    return this.equals(P.negate());
};

/*
    * Compute this + Q.
    * @param {JSPoint2} Q an elliptic curve point.
    * @returns this + Q.
    */
JSPoint2.add = function (Q) {
    if (this.isZero()) {
        return Q;
    }
    if (Q.isZero()) {
        return this;
    }
    //*
    // EFD addition formulas:
    // <http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html>
    var Fp2_1 = this.E.E.bn.Fp2_1;
    var X1 = this.x, Y1 = this.y, Z1 = this.z, X2 = Q.x, Y2 = Q.y, Z2 = Q.z, Z1Z1 = Fp2_1, Z2Z2 = Fp2_1, U1 = this.x, U2 = Q.x, S1 = this.y, S2 = Q.y, H, I, J, R, V, X3, Y3, Z3;
    var Z1is1 = Z1.isOne();
    var Z2is1 = Z2.isOne();
    if (!Z1is1) {
        Z1Z1 = Z1.square(); // Z1Z1 = Z1^2
        U2 = X2.multiply(Z1Z1); // U2 = X2*Z1Z1
        S2 = Y2.multiply(Z1).multiply(Z1Z1); // S2 = Y2*Z1*Z1Z1
    }
    if (!Z2is1) {
        Z2Z2 = Z2.square(); // Z2Z2 = Z2^2
        U1 = X1.multiply(Z2Z2); // U1 = X1*Z2Z2
        S1 = Y1.multiply(Z2).multiply(Z2Z2); // S1 = Y1*Z2*Z2Z2
    }
    if (U1.equals(U2)) {
        if (S1.equals(S2)) {
            return this.twice(1);
        } else {
            return this.E.infinity;
        }
    }
    H = U2.subtract(U1); // H = U2-U1
    I = H.twice(1).square(); // I = (2*H)^2
    J = H.multiply(I); // J = H*I
    R = S2.subtract(S1).twice(1); // r = 2*(S2-S1)
    V = U1.multiply(I); // V = U1*I
    X3 = R.square().subtract(J).subtract(V.twice(1)); // X3 = r^2-J-2*V
    Y3 = R.multiply(V.subtract(X3)).subtract(S1.multiply(J).twice(1)); // Y3 = r*(V-X3)-2*S1*J
    if (Z2is1) {
        if (Z1is1) {
            Z3 = H.twice(1); // Z3 = 2*H
        } else {
            Z3 = Z1.multiply(H).twice(1); // Z3 = ((Z1+1)^2-Z1Z1-1)*H
        }
    } else {
        if (Z1is1) {
            Z3 = Z2.multiply(H).twice(1); // Z3 = ((1+Z2)^2-1-Z2Z2)*H
        } else {
            Z3 = Z1.add(Z2).square().subtract(Z1Z1).subtract(Z2Z2).multiply(H);// Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2)*H
        }
    }
    return new JSPoint2(this.E, X3, Y3, Z3);
};

/*
    * Left-shift this point by a given distance n, i.e. compute (2^^n)*this.
    * @param {Integer} n the shift amount.
    * @returns (2^^n)*this.
    */
JSPoint2.twice = function (n) {
    if (this.isZero()) return this;
    // EDF doubling formulas:
    // <http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l>
    // <http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-mdbl-2007-bl>
    var A, B, C, S, M, X = this.x, Y = this.y, Z = this.z;
    while (n-- > 0) {
        A = X.square(); // A = X1^2
        B = Y.square(); // B = Y1^2
        C = B.square(); // C = B^2
        S = X.add(B).square().subtract(A).subtract(C).twice(1); // S = 2*((X1+B)^2-A-C)
        M = A.multiply(new Number(3)); // M = 3*A
        X = M.square().subtract(S.twice(1)); // X3 = M^2-2*S
        Z = Y.multiply(Z).twice(1); // Z3 = 2*Y1*Z1
        Y = M.multiply(S.subtract(X)).subtract(C.twice(3)); // Y3 = M*(S-X3)-8*C
    }
    return new JSPoint2(this.E, X, Y, Z);
};

/*
    * Compute k*this
    * This method implements the GLV strategy.
    * @param {BigInteger} k scalar by which this point is to be multiplied.
    * @returns k*this.
    */
JSPoint2.multiply = function (k) {
    var bn = this.E.E.bn;
    var P = this.normalize();
    if (k.signum() < 0) {
        k = k.negate();
        P = P.negate();
    }
    var halfn = bn.n.shiftRight(1);
    var w = new Array(4);
    for (var i=0; i<4; i++) {
        w[i] = k.multiply(bn.latInv[i]);
        if (w[i].mod(bn.n).compareTo(halfn) <= 0) {
            w[i] = w[i].divide(bn.n);
        } else {
            w[i] = w[i].divide(bn.n).add(bn._1);
        }
    }
    var u = new Array(4);
    for (var j=0; j<4; j++) {
        u[j] = bn._0;
        for (var i=0; i<4; i++) {
            u[j] = u[j].add(w[i].multiply(bn.latRed[i][j]));
        }
        u[j] = u[j].negate();
    }
    u[0] = u[0].add(k);
    var Q = P.frobex(1);
    var R = P.frobex(2);
    var S = P.frobex(3);
    return this.simultaneous(u[0], P, u[1], Q, u[2], R, u[3], S);
};

JSPoint2.simultaneous = function (kP, P, kQ, Q, kR, R, kS, S) {
    var hV = new Array(16);
    P = P.normalize();
    if (kP.signum() < 0) {
        kP = kP.negate(); P = P.negate();
    }
    Q = Q.normalize();
    if (kQ.signum() < 0) {
        kQ = kQ.negate(); Q = Q.negate();
    }
    R = R.normalize();
    if (kR.signum() < 0) {
        kR = kR.negate(); R = R.negate();
    }
    S = S.normalize();
    if (kS.signum() < 0) {
        kS = kS.negate(); S = S.negate();
    }
    hV[0] = this.E.infinity;
    hV[1] = P; hV[2] = Q; hV[4] = R; hV[8] = S;
    for (var i = 2; i < 16; i <<= 1) {
        for (var j = 1; j < i; j++) {
            hV[i + j] = hV[i].add(hV[j]);
        }
    }
    var t = Math.max(Math.max(kP.bitLength(), kQ.bitLength()), Math.max(kR.bitLength(), kS.bitLength()));
    var V = this.E.infinity;
    for (var i=t-1; i>=0; i--) {
        var j = (kS.testBit(i)?8:0) | (kR.testBit(i)?4:0) | (kQ.testBit(i) ?   2 : 0) | (kP.testBit(i)?1:0);
        V = V.twice(2).add(hV[j]);
    }
    return V.normalize();
};

/*
    * 
    */
JSPoint2.frobex = function (k) {
    if (!this.z.isOne()) {
        throw new Error("Logic Error!");
    }
    var bn = this.E.E.bn;
    switch (k) {
        case 1:
            return (bn.b === 3) ?
                new JSPoint2(this.E, this.x.multiplyI().conjugate().multiply(bn.zeta), this.y.multiplyV().conjugate().multiply(bn.zeta0sigma), this.z) :
                new JSPoint2(this.E, this.x.conjugate().multiplyI().multiply(bn.zeta), this.y.conjugate().multiplyV().multiply(bn.zeta1sigma), this.z);
        case 2:
            return new JSPoint2(this.E, this.x.multiply(bn.zeta1).negate(), this.y.negate(), this.z);
        case 3:
            return (bn.b === 3) ?
                new JSPoint2(this.E, this.x.multiplyI().conjugate(), this.y.multiplyV().conjugate().multiply(bn.zeta0sigma).negate(), this.z) :
                new JSPoint2(this.E, this.x.conjugate().multiplyI(), this.y.conjugate().multiplyV().multiply(bn.zeta1sigma).negate(), this.z);
        default:
            return null;
    }
};

JSPoint2.prototype.isZero = JSPoint2.isZero;
JSPoint2.prototype.equals = JSPoint2.equals;
JSPoint2.prototype.isOnSameCurve = JSPoint2.isOnSameCurve;
JSPoint2.prototype.randomize = JSPoint2.randomize;
JSPoint2.prototype.normalize = JSPoint2.normalize;
JSPoint2.prototype.negate = JSPoint2.negate;
JSPoint2.prototype.opposite = JSPoint2.opposite;
JSPoint2.prototype.add = JSPoint2.add;
JSPoint2.prototype.twice = JSPoint2.twice;
JSPoint2.prototype.multiply = JSPoint2.multiply;
JSPoint2.prototype.simultaneous = JSPoint2.simultaneous;
JSPoint2.prototype.frobex = JSPoint2.frobex;

module.exports = JSPoint2;