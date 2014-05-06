/*
 * JSPoint.js
 * Arithmetic in the group of points on a BN elliptic curve over GF(p).
 * Copyright (C) Jorge H. C. Lucema.
 */

var JSPoint = (function () {

    var JSPoint = function (E, x, y, z) {

        this.pp16P = null;
        this.preComp = null;

        if (arguments.length === 1) {
            /*
             * Create an instance of the BNCurve point at infinity on curve E.
             * @param {BNCurve} E the elliptic curve where the created point is located.
             */
            if (E instanceof JSCurve) {
                this.E = E;
                /*
                 * the point at infinity is represented as (1, 1, 0) after IEEE Std 1363:2000
                 * (notice that this triple satisfies the projective curve equation y^2 = x^3 + b.z^6)
                 */
               // The projective x-coordinate
               this.x = E.bn._1;
               // The projective y-coordinate
               this.y = E.bn._1;
               // The projective z-coordinate
               this.z = E.bn._0;
            }
            /*
             * Create a clone of a given point.
             * @param {JSPoint} Q the point to be cloned.
             */
            if (E instanceof JSPoint) {
                var Q = E;
                this.E = Q.E;
                this.x = Q.x;
                this.y = Q.y;
                this.z = Q.z;
            }
        }

        if (arguments.length === 3) {
            /*
             * Create a normalized BNCurve point from given affine coordinates and a curve
             * @param {BNCurve} E the underlying elliptic curve.
             * @param {BigInteger} x the affine x-coordinate (mod p).
             * @param {BigInteger} y the affine y-coordinate (mod p).
             */
            if ((x instanceof BigInteger) && (y instanceof BigInteger)) {
                this.E = E;
                var p = E.bn.p; // shorthand
                this.x = x.mod(p);
                this.y = y.mod(p);
                this.z = E.bn._1; // normalized
            }
            /*
             * Create an BNCurve point from a given affine x-coordinate, a y-bit, and a curve
             * @param {BNCurve} E the underlying elliptic curve.
             * @param {BigInteger} x the affine x-coordinate (mod p).
             * @param {Integer} y the least significant bit of the y-coordinate.
             */
            else if (x instanceof BigInteger) {
                var yBit = y;
                this.E = E;
                var p = E.bn.p; // shorthand
                this.x = x.mod(p);
                if (x.signum() === 0) throw new Error("The given point does not belong to the given elliptic curve"); // otherwise the curve order would not be prime
                else {
                    this.y = E.bn.sqrt(x.multiply(x).multiply(x).add(E.b).mod(p));
                    if (this.y === null) throw new Error("The given point does not belong to the given elliptic curve");
                    if (this.y.testBit(0) !== ((yBit & 1) === 1)) this.y = p.subtract(y);
                }
                this.z = E.bn._1; // normalized
            }
            /*
             * Create an BNCurve point from a given x-trit, an affine y-coordinate, and a curve
             * @param {BNCurve} E the underlying elliptic curve.
             * @param {Integer} x the least significant trit of the x-coordinate.
             * @param {BigInteger} y the affine y-coordinate (mod p).
             */
            else if (y instanceof BigInteger) {
                var xTrit = x;
                this.E = E;
                var p = E.bn.p; // shorthand
                this.y = y.mod(p);
                if (y.signum() === 0) throw new Error("The given point does not belong to the given elliptic curve"); // otherwise the curve order would not be prime
                else {
                    this.x = E.bn.cbrt(y.multiply(y).subtract(E.b).mod(p));
                    if (this.x === null) throw new Error("The given point does not belong to the given elliptic curve");
                    // either x, zeta*x, or zeta^2*x is the desired x-coordinate:
                    if (this.x.mod(E.bn._3).intValue() !== xTrit) {
                        var zeta = E.bn.zeta; // shorthand
                        this.x = zeta.multiply(x).mod(p);
                        if (this.x.mod(E.bn._3).intValue() !== xTrit) {
                            this.x = zeta.multiply(x).mod(p);
                            if (this.x.mod(E.bn._3).intValue() !== xTrit) throw new Error("The given point does not belong to the given elliptic curve");
                        }
                    }
                }
                this.z = E.bn._1; // normalized
            }
        }

        /*
         * Create an BNCurve point from given projective coordinates and a curve.
         * @param {BNCurve} E the underlying elliptic curve.
         * @param {BigInteger} x the affine x-coordinate (mod p).
         * @param {BigInteger} y the affine y-coordinate (mod p).
         * @param {BigInteger} z the affine z-coordinate (mod p).
         */
        if (arguments.length === 4) {
            this.E = E;
            this.x = x;
            this.y = y;
            this.z = z;
        }

    };

    /*
     * Check whether this is the point at infinity (i.e. the BNCurve group zero element).
     * @returns true if this is the point at infinity, otherwise false.
     */
    JSPoint.isZero = function () {
        return this.z.signum() === 0;
    };

    /*
     * Compare this point to a given object.
     * @param {Object} Q the elliptic curve point to be compared to this.
     * @returns true if this point and Q are equal, otherwise false.
     */
    JSPoint.equals = function (Q) {
        if (!(Q instanceof JSPoint && this.isOnSameCurve(Q))) return false;
        var P = Q;
        if (this.z.signum() === 0 || P.z.signum() === 0) return this.z.equals(P.z);
        var p = this.E.bn.p; // shorthand
        var z2 = this.z.multiply(this.z).mod(p), z3 = this.z.multiply(z2).mod(p),
            pz2 = P.z.multiply(P.z).mod(p), pz3 = P.z.multiply(pz2).mod(p);
        return this.x.multiply(pz2).subtract(P.x.multiply(z2)).mod(p).signum() === 0 && this.y.multiply(pz3).subtract(P.y.multiply(z3)).mod(p).signum() === 0;
    };

    /*
     * Check whether Q lays on the same curve as this point.
     * @param {JSPoint} Q an elliptic curve point.
     * @returns true if Q lays on the same curve as this point, otherwise false.
     */
    JSPoint.isOnSameCurve = function (Q) {
        return this.E.bn === Q.E.bn; // singleton comparison
    };

    /*
     * Compute a random point on the same curve as this.
     * @param {SecureRandom} rand a cryptographically strong pseudo-random number generator.
     * @returns a random point on the same curve as this.
     */
    JSPoint.randomize = function (rand) {
        if (rand instanceof SecureRandom) return this.E.pointFactory(rand);
    };

    /*
     * Normalize this point.
     * @returns a normalized point equivalent to this.
     */
    JSPoint.normalize = function () {
        if (this.z.signum() === 0 || this.z.compareTo(this.E.bn._1) === 0) return this; // already normalized
        var p = this.E.bn.p; // shorthand
        var zinv = null;
        try {
            zinv = this.z.modInverse(p);
        } catch (a) {
            //
        }
        var zinv2 = zinv.multiply(zinv); // mod p
        return new JSPoint(this.E, this.x.multiply(zinv2).mod(p), this.y.multiply(zinv).multiply(zinv2).mod(p), this.E.bn._1);
    };

    /*
     * 
     * @returns true if z is one
     */
    JSPoint.isNormal = function () {
        return (this.z.signum() === 0 || this.z.compareTo(E.bn._1) === 0);
    };

    /*
     * Compute -this.
     * @return -this.
     */
    JSPoint.negate = function () {
        return new JSPoint(this.E, this.x, (this.y.signum() !== 0) ? this.E.bn.p.subtract(this.y) : this.y, this.z);
    };

    /*
     * Check if a point equals -this.
     * @param {JSPoint} P
     * @returns true if P=-this
     */
    JSPoint.opposite = function (P) {
        if (!isOnSameCurve(P)) return false;
        if (this.z.signum() === 0 || P.isZero()) return this.z.compareTo(P.z) === 0;
        var p = E.bn.p; // shorthand
        var z2 = this.z.multiply(this.z);
        var z3 = this.z.multiply(z2).mod(p);
        var pz2 = P.z.multiply(P.z);
        var pz3 = P.z.multiply(pz2).mod(p);
        return this.x.multiply(pz2).subtract(P.x.multiply(z2)).mod(p).signum() === 0 &&
                this.y.multiply(pz3).add(P.y.multiply(z3)).mod(p).signum() === 0;
    };

    /*
     * Compute this + Q.
     * @param {JSPoint} Q an elliptic curve point.
     * @return this + Q.
     */
    JSPoint.add = function (Q) {
        if (this.isZero()) return Q;
        if (Q.isZero()) return this;
        //*
        // EFD addition formulas:
        // <http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html>
        var p = this.E.bn.p;
        var X1 = this.x, Y1 = this.y, Z1 = this.z, U1 = this.x, S1 = this.y, Z1Z1;
        var X2 = Q.x, Y2 = Q.y, Z2 = Q.z, U2 = Q.x, S2 = Q.y, Z2Z2, H, I, J, R, V, X3, Y3, Z3;
        var Z1is1 = (Z1.compareTo(this.E.bn._1) === 0);
        var Z2is1 = (Z2.compareTo(this.E.bn._1) === 0);
        if (!Z1is1) {
            Z1Z1 = Z1.multiply(Z1).mod(p); // Z1Z1 = Z1^2
            U2 = X2.multiply(Z1Z1).mod(p); // U2 = X2*Z1Z1
            S2 = Y2.multiply(Z1).multiply(Z1Z1).mod(p); // S2 = Y2*Z1*Z1Z1
        }
        if (!Z2is1) {
            Z2Z2 = Z2.multiply(Z2).mod(p); // Z2Z2 = Z2^2
            U1 = X1.multiply(Z2Z2).mod(p); // U1 = X1*Z2Z2
            S1 = Y1.multiply(Z2).multiply(Z2Z2).mod(p); // S1 = Y1*Z2*Z2Z2
        }
        if (U1.compareTo(U2) === 0)
            if (S1.compareTo(S2) === 0) return this.twice(1);
            else return this.E.infinity;
        H = U2.subtract(U1); // H = U2-U1
        I = H.shiftLeft(1);
        I = I.multiply(I).mod(p); // I = (2*H)^2
        J = H.multiply(I);//.mod(p); // J = H*I
        R = S2.subtract(S1).shiftLeft(1); // r = 2*(S2-S1)
        V = U1.multiply(I);//.mod(p); // V = U1*I
        X3 = R.multiply(R).subtract(J).subtract(V.shiftLeft(1)).mod(p); // X3 = r^2-J-2*V
        Y3 = R.multiply(V.subtract(X3)).subtract(S1.multiply(J).shiftLeft(1)).mod(p); // Y3 = r*(V-X3)-2*S1*J
        if (Z2is1)
            if (Z1is1) Z3 = H.shiftLeft(1).mod(p); // Z3 = 2*H
            else Z3 = Z1.multiply(H).shiftLeft(1).mod(p); // Z3 = ((Z1+1)^2-Z1Z1-1)*H = (Z1+H)^2-Z1Z1-H^2
        else
            if (Z1is1) Z3 = Z2.multiply(H).shiftLeft(1).mod(p); // Z3 = ((1+Z2)^2-1-Z2Z2)*H = (H+Z2)^2-H^2-Z2Z2
            else {
                Z3 = Z1.add(Z2);
                Z3 = Z3.multiply(Z3).subtract(Z1Z1).subtract(Z2Z2).multiply(H).mod(p);// Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2)*H
            }
        return new JSPoint(this.E, X3, Y3, Z3);
    };

    /*
     * Compute this - Q.
     * @param {JSPoint} Q an elliptic curve point.
     * @return this - Q.
     */
    JSPoint.subtract = function (Q) {
        return this.add(Q.negate());
    };

    /*
     * Left-shift this point by a given distance n, i.e. compute (2^^n)*this.
     * @param {Integer} n the shift amount.
     * @return (2^^n)*this.
     */
    JSPoint.twice = function (n) {
        if (this.isZero()) return this;
        // EDF doubling formulas:
        // <http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l>
        // <http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-mdbl-2007-bl>
        var A, B, C, S, M, X = this.x, Y = this.y, Z = this.z;
        var p = this.E.bn.p; // shorthand
        while (n-- > 0) {
            A = X.multiply(X); // A = X1^2 (modular reduction is irrelevant)
            B = Y.multiply(Y).mod(p); // B = Y1^2
            C = B.multiply(B); // C = B^2 (modular reduction is irrelevant)
            S = X.add(B);
            S = S.multiply(S).subtract(A).subtract(C).shiftLeft(1).mod(p); // S = 2*((X1+B)^2-A-C)
            M = A.multiply(this.E.bn._3).mod(p); // M = 3*A
            X = M.multiply(M).subtract(S.shiftLeft(1)).mod(p); // X3 = M^2-2*S
            Z = Y.multiply(Z).shiftLeft(1).mod(p); // Z3 = 2*Y1*Z1
            Y = M.multiply(S.subtract(X)).subtract(C.shiftLeft(3)).mod(p); // Y3 = M*(S-X3)-8*C
        }
        return new JSPoint(this.E, X, Y, Z);
    };

    /*
     * Compute k*this
     * This method implements the GLV strategy when no precomputed table is available,
     * otherwise the quaternary window method with precomputation.
     * @param {BigInteger} k scalar by which this point is to be multiplied
     * @return  k*this
     */
    JSPoint.multiply = function (k) {
        if(this.pp16P === null) {
            if(k.compareTo(BigInteger.ONE) === 0) {
                return this;
            }
            var bn = this.E.bn;
            var P = this.normalize();
            if (k.signum() < 0) {
                k = k.negate();
                P = P.negate();
            }
            var r = bn.u.shiftLeft(1).add(this.E.bn._1); // 2*u + 1
            var t = bn.u.multiply(this.E.bn._3).add(this.E.bn._1).multiply(bn.u.shiftLeft(1)); // (3*u + 1)*2*u = 6*u^2 + 2*u
            //
            var halfn = bn.n.shiftRight(1);
            var kr = k.multiply(r);
            if (kr.mod(bn.n).compareTo(halfn) <= 0) kr = kr.divide(bn.n);
            else kr = kr.divide(bn.n).add(this.E.bn._1);
            var kt = k.multiply(t);
            if (kt.mod(bn.n).compareTo(halfn) <= 0) kt = kt.divide(bn.n);
            else kt = kt.divide(bn.n).add(this.E.bn._1);
            //
            var sr = k.subtract(kr.multiply(r).add(kt.multiply(t.add(r))));
            var st = kr.multiply(t).subtract(kt.multiply(r));
            //
            var Y = new JSPoint(this.E, P.x.multiply(bn.zeta), P.y, P.z);
            return P.simultaneous(sr, st, Y);
        } else {
            k = k.mod(this.E.bn.n);
            var A = this.E.infinity;
            for (var i = 0, w = 0; i < this.pp16P.length; i++, w >>>= 4) {
                if ((i & 7) === 0) {
                    w = k.intValue();
                    k = k.shiftRight(32);
                }
                A = A.add(this.pp16P[i][w & 0xf]);
            }
            return A;
        }
    };

    JSPoint.simultaneous = function (ks, kr, Y, Q, kR, R, kS, S) {
        /*
         * Compute ks*this + kr*Y.  This is useful in the verification part of several signature algorithms,
         * and (hopely) faster than two scalar multiplications.
         * @param {BigInteger} ks scalar by which this point is to be multiplied.
         * @param {BigInteger} kr scalar by which Y is to be multiplied.
         * @param {JSPoint} Y a curve point.
         * @returns ks*this + kr*Y
         */
        if (arguments.length === 3) {
            if (this.pp16P === null) {
                var hV = new Array(16);
                var P = this.normalize();
                Y = Y.normalize();
                if (ks.signum() < 0) {
                    ks = ks.negate();
                    P = P.negate();
                }
                if (kr.signum() < 0) {
                    kr = kr.negate();
                    Y = Y.negate();
                }
                hV[0] = this.E.infinity;
                hV[1] = P;
                hV[2] = Y;
                hV[3] = P.add(Y);
                for (var i = 4; i < 16; i += 4) {
                    hV[i] = hV[i >> 2].twice(1);
                    hV[i + 1] = hV[i].add(hV[1]);
                    hV[i + 2] = hV[i].add(hV[2]);
                    hV[i + 3] = hV[i].add(hV[3]);
                }
                var t = Math.max(kr.bitLength(), ks.bitLength());
                var R = this.E.infinity;
                for (var i = (((t + 1) >> 1) << 1) - 1; i >= 0; i -= 2) {
                    var j = (kr.testBit(i  ) ? 8 : 0) |
                            (ks.testBit(i  ) ? 4 : 0) |
                            (kr.testBit(i-1) ? 2 : 0) |
                            (ks.testBit(i-1) ? 1 : 0);
                    R = R.twice(2).add(hV[j]);
                }
                return R;
            } else R = this.multiply(ks).add(Y.multiply(ks));
        }
        /*
         * 
         */
        if (arguments.length === 8) {
            var kP = ks, P = kr, kQ = Y;
            var hV = new Array(16);
            P = P.normalize();
            if (kP.signum() < 0) {
                kP = kP.negate();
                P = P.negate();
            }
            Q = Q.normalize();
            if (kQ.signum() < 0) {
                kQ = kQ.negate();
                Q = Q.negate();
            }
            R = R.normalize();
            if (kR.signum() < 0) {
                kR = kR.negate();
                R = R.negate();
            }
            S = S.normalize();
            if (kS.signum() < 0) {
                kS = kS.negate();
                S = S.negate();
            }
            hV[0] = this.E.infinity;
            hV[1] = P;
            hV[2] = Q;
            hV[4] = R;
            hV[8] = S;
            for (var i = 2; i < 16; i <<= 1) {
                for (var j = 1; j < i; j++) {
                    hV[i + j] = hV[i].add(hV[j]);
                }
            }
            var t = Math.max(Math.max(kP.bitLength(), kQ.bitLength()), Math.max(kR.bitLength(), kS.bitLength()));
            var V = this.E.infinity;
            for (var i = t - 1; i >= 0; i--) {
                var j = (kS.testBit(i) ?   8 : 0) |
                        (kR.testBit(i) ?   4 : 0) |
                        (kQ.testBit(i) ?   2 : 0) |
                        (kP.testBit(i) ?   1 : 0);
                V = V.twice(1).add(hV[j]);
            }
            return V;
        }
    };

    JSPoint.getSerializedTable = function () {
        if (this.pp16P === null) {
            var length = Math.floor((this.E.bn.n.bitLength() + 3)/4);
            this.pp16P = new Array(length);
            this.pp16Pserial = new Array(length);
            for (var i = 0; i < length; i++) {
                this.pp16P[i] = new Array(16);
                this.pp16Pserial[i] = new Array(16);
            }
            var P = this.normalize();
            var pp16Pi = this.pp16P[0];
            pp16Pi[0] = this.E.infinity;
            pp16Pi[1] = P;
            for (var i = 1, j = 2; i <= 7; i++, j += 2) {
                pp16Pi[j] = pp16Pi[i].twice(1).normalize();
                pp16Pi[j+1] = pp16Pi[j].add(P).normalize();
            }
            for (var i = 1; i < this.pp16P.length; i++) {
                var pp16Ph = pp16Pi;
                pp16Pi = this.pp16P[i];
                pp16Pi[0] = pp16Ph[0];
                for (var j = 1; j < 16; j++) {
                    pp16Pi[j] = pp16Ph[j].twice(4).normalize();
                }
            }
            for (var i = 0; i < this.pp16P.length; i++) {
                for (var j = 0; j < this.pp16P[i].length; j++) {
                    this.pp16Pserial[i][j] = this.pp16P[i][j].toByteArray(4);
                }
            }
        }
        return this.pp16Pserial;
    };

    JSPoint.preCompute = function () {
        if (this.preComp === null) {
            // width=8 => 2^(8-1)-1=127
            this.preComp = new Array(127);
            this.preComp[0] = this;
            // Compute 2*P
            var twiceP = this.twice(1);
            for (var i = 1; i < 127; i++) {
                // Compute the new ECPoints for the precomputation array.
                // The values 1, 3, 5, ..., 2^(width-1)-1 times p are computed.
                this.preComp[i] = twiceP.add(this.preComp[i - 1]);
            }
        }
    };

    JSPoint.preCompute2 = function () {
        if (this.preComp === null) {
            var l = this.E.bn.n.bitLength();
            var w;
            if (l < 13) w = 2;
            else if (l < 41) w = 3;
            else if (l < 121) w = 5;
            else if (l < 337) w = 7;
            else if (l < 897) w = 7;
            else if (l < 2305) w = 8;
            else w = 8;
            var a = Math.ceil(l/w);
            var v;
            if (a%2 === 0) {
                v = Math.ceil(a/2);
            } else if (a%3 === 0) {
                v = Math.ceil(a/3);
            } else if (a%5 === 0) {
                v = Math.ceil(a/5);
            } else if (a%7 === 0) {
                v = Math.ceil(a/7);
            } else if (a%9 === 0) {
                v = Math.ceil(a/9);
            } else if (a%11 === 0) {
                v = Math.ceil(a/11);
            } else if (a%13 === 0) {
                v = Math.ceil(a/13);
            } else {
                v = 1;
            }
            // si a es divisible por 2, osea: a%2=0, entonces v=a/2
            // caso contrario probar para a%3, a%5, ... cuando cero entonces la division es v
            var b = Math.ceil(a/v);
            var s = new Array(w);
            for (var i=0; i<w; i++) s[i] = Math.pow(2,i);
            var d = new Array(Math.pow(2,w-2));
            for (var i=0; i<d.length; i++) d[i] = 2*i+1;
            var P = this.normalize();
            this.preComp = new Array(v);
            for (var i=0; i<v; i++) this.preComp[i] = new Array();
            for (var i=0; i<s.length; i++) {
                for (var j=0; j<d.length; j++) {
                    this.preComp[0][s[i]*d[j]] = P.multiply(new BigInteger((s[i]*d[j]).toString(),10));
                }
            }
            for (var i=0; i<s.length; i++) {
                for (var j=0; j<d.length; j++) {
                    for (var k=1; k<v; k++) {
                        this.preComp[k][s[i]*d[j]] = this.preComp[0][s[i]*d[j]].twice(k*w*b);
                    }
                }
            }
        }
    };

    JSPoint.toByteArray = function (formFlags) {
        var len = Math.floor((this.E.bn.p.bitLength() + 7)/8);
        var resLen = 1, pc = 0;
        var P = this.normalize();
        var osX = null, osY = null;
        if (!P.isZero()) {
            osX = P.x.toByteArray();
            resLen += len;
            if ((formFlags & 2) !== 0) {
                pc |= 2 | (P.y.testBit(0) ? 1 : 0);
            }
            if ((formFlags & 4) !== 0) {
                pc |= 4;
                osY = P.y.toByteArray();
                resLen += len;
            }
        }
        var buf = new Uint8Array(resLen);
        for (var i = 0; i < buf.length; i++) {
            buf[i] = 0;
        }
        buf[0] = pc;
        if (osX !== null) {
            if (osX.length <= len) {
                JSPoint.arrayCopy(osX,0,buf,1+len-osX.length,osX.length);
            } else {
                JSPoint.arrayCopy(osX,1,buf,1,len);
            }
        }
        if (osY !== null) {
            if (osY.length <= len) {
                JSPoint.arrayCopy(osY,0,buf,1+2*len-osY.length,osY.length);
            } else {
                JSPoint.arrayCopy(osY,1,buf,1+len,len);
            }
        }
        return buf;
    };

    JSPoint.arrayCopy = function (src, srcPos, dest, destPos, length) {
        for (var i = srcPos, j = destPos; i < length + srcPos; i++, j++) {
            if (dest[j] !== "undefined") {
                dest[j] = src[i];
            } else {
                throw "array index out of bounds exception";
            }
        }
    };

    JSPoint.windowNAF = function (w, k) {
        // Clamp the window width in the range [2,8]
        if (w < 2 || w > 8) throw new Error("'width' must be in the range [2,8]");
        //The window NAF is at most 1 element longer than the binary representation of the integer k.
        var wnaf = new Int8Array(this.E.bn.n.bitLength() + 1);
        // 2^width as short and BigInteger
        var pow2wB = (1 << w);
        var pow2wBI = new BigInteger(pow2wB.toString(), 10);
        var i = 0;
        // The actual length of the WNAF
        var length = 0;
        // while k >= 1
        while (k.signum() > 0) {
            // if k is odd
            if (k.testBit(0)) {
                // ki = k mod 2^width
                var remainder = k.mod(pow2wBI);
                // if ki > 2^(width - 1) - 1
                if (remainder.testBit(w - 1)) {
                    wnaf[i] = (remainder.intValue() - pow2wB);
                } else {
                    wnaf[i] = remainder.intValue();
                }
                // wnaf[i] is now in [-2^(width-1), 2^(width-1)-1] => ki
                // k = k - ki
                k = k.subtract(new BigInteger(wnaf[i].toString(), 10));
                length = i;
            } else {
                wnaf[i] = 0;
            }
            // k = k/2
            k = k.shiftRight(1);
            i++;
        }
        length++;
        // Reduce the WNAF array to its actual length
//        var wnafShort = new Int8Array(length);
//        JSPoint.arrayCopy(wnaf,0,wnafShort,0,length);
        if (length < this.E.bn.n.bitLength()) {
            var wnafLong = new Int8Array(this.E.bn.n.bitLength());
            JSPoint.arrayCopy(wnaf,0,wnafLong,0,this.E.bn.n.bitLength());
            return wnafLong;
        }
        return wnaf;
    };

    JSPoint.multiplyNAF = function (k) {
        // floor(log2(k))
        var m = k.bitLength();
        // width of the Window NAF
        var w;
        // Required length of precomputation array
//        var reqPreqComp;
        // Determine optimal width and corresponding length of precomputation
        // array based on literature values
        if (m < 13) w = 2;
        else if (m < 41) w = 3;
        else if (m < 121) w = 5;
        else if (m < 337) w = 7;
        else if (m < 897) w = 7;
        else if (m < 2305) w = 8;
        else w = 8;
        // Compute the Window NAF of the desired width
        var wnaf = this.windowNAF(w,k);
//        var preComp = new Array(reqPreqComp);
//        preComp[0] = this;
//        // Compute 2*P
//        var twiceP = this.twice(1);
//        for (var i = 1; i < reqPreqComp; i++) {
//            // Compute the new ECPoints for the precomputation array.
//            // The values 1, 3, 5, ..., 2^(width-1)-1 times p are computed.
//            preComp[i] = twiceP.add(preComp[i - 1]);
//        }
        // Q point infinite
        var Q = this.E.infinity;
        var l = wnaf.length;
        // Apply the Window NAF to P using the precomputed ECPoint values.
        for (var i = l - 1; i >= 0; i--) {
            Q = Q.twice(1);
            if (wnaf[i] !== 0) {
                if (wnaf[i] > 0) {
                    Q = Q.add(this.preComp[(wnaf[i] - 1)/2]);
                } else {
                    // wnaf[i] < 0
                    Q = Q.subtract(this.preComp[(-wnaf[i] - 1)/2]);
                }
            }
        }
        return Q;
    };

    JSPoint.proposedMultiply = function (k) {
        // floor(log2(k))
        var m = this.E.bn.n.bitLength();
        // width of the Window NAF
        var w;
        // Determine optimal width and corresponding length of precomputation
        // array based on literature values
        if (m < 13) w = 2;
        else if (m < 41) w = 3;
        else if (m < 121) w = 5;
        else if (m < 337) w = 7;
        else if (m < 897) w = 7;
        else if (m < 2305) w = 8;
        else w = 8;
        // Compute the Window NAF of the desired width
        var wnafLong = this.windowNAF(w,k);
        var l = wnafLong.length;
        var wnaf, v;
        // si a es divisible por 2, osea: a%2=0, entonces v=a/2
        // caso contrario probar para a%3, a%5, ... cuando cero entonces la division es v
        var a = Math.ceil(l/w);
        if (a%2 === 0) {
            v = Math.ceil(a/2);
        } else if (a%3 === 0) {
            v = Math.ceil(a/3);
        } else if (a%5 === 0) {
            v = Math.ceil(a/5);
        } else if (a%7 === 0) {
            v = Math.ceil(a/7);
        } else if (a%9 === 0) {
            v = Math.ceil(a/9);
        } else if (a%11 === 0) {
            v = Math.ceil(a/11);
        } else if (a%13 === 0) {
            v = Math.ceil(a/13);
        } else {
            v = 1;
        }
        var b = Math.ceil(a/v);
        var zeros = a*w-l;
        if (zeros !== 0) {
            wnaf = new Int8Array(l+zeros);
            JSPoint.arrayCopy(wnafLong,0,wnaf,0,l+zeros);
        } else {
            wnaf = wnafLong;
        }
        var e = new Array(w);
        for (var k=0; k<w; k++) e[k]=new Array(a);
        var x=0;
        for (var j=0; j<a; j++) {
            for (var i=0; i<w; i++) {
                e[i][j] = wnaf[x];
                x++;
            }
        }
        // Q point infinite
        var Q = this.E.infinity;
        // Apply Proposed Method
        var temp;
        for (var t = b - 1; t >= 0; t--) {
            Q = Q.twice(w);
            for (var j=v-1; j>=0; j--) {
                temp=0;
                for (k=0; k<w; k++) {
                    if(e[k][j*b+t] !== 0) temp+=e[k][j*b+t]*Math.pow(2,k);
                }
                if (temp>0) {
                    Q=Q.add(this.preComp[j][temp]);
                } else if (temp<0) {
                    Q=Q.subtract(this.preComp[j][-temp]);
                }
            }
        }
        return Q;
    };

    JSPoint.prototype.isZero = JSPoint.isZero;
    JSPoint.prototype.equals = JSPoint.equals;
    JSPoint.prototype.isOnSameCurve = JSPoint.isOnSameCurve;
    JSPoint.prototype.randomize = JSPoint.randomize;
    JSPoint.prototype.normalize = JSPoint.normalize;
    JSPoint.prototype.isNormal = JSPoint.isNormal;
    JSPoint.prototype.negate = JSPoint.negate;
    JSPoint.prototype.opposite = JSPoint.opposite;
    JSPoint.prototype.add = JSPoint.add;
    JSPoint.prototype.subtract = JSPoint.subtract;
    JSPoint.prototype.twice = JSPoint.twice;
    JSPoint.prototype.multiply = JSPoint.multiply;
    JSPoint.prototype.simultaneous = JSPoint.simultaneous;
    JSPoint.prototype.getSerializedTable = JSPoint.getSerializedTable;
    JSPoint.prototype.preCompute = JSPoint.preCompute;
    JSPoint.prototype.preCompute2 = JSPoint.preCompute2;
    JSPoint.prototype.toByteArray = JSPoint.toByteArray;
    JSPoint.prototype.windowNAF = JSPoint.windowNAF;
    JSPoint.prototype.multiplyNAF = JSPoint.multiplyNAF;
    JSPoint.prototype.proposedMultiply = JSPoint.proposedMultiply;

    return JSPoint;
})();
