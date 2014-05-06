/*
 * JSField6.js
 * Arithmetic in the finite extension field GF(p^6) with p = 3 (mod 4).
 * Copyright (C) Jorge H. C. Lucema.
 */

var differentFields = "Operands are in different finite fields";

var JSField6 = (function () {

    var JSField6 = function (bn, k, v1, v2) {
        /*
         * 
         */
        if (arguments.length === 1) {
            this.bn = bn;
            this.v = new Array(3);
            this.v[0] = this.v[1] = this.v[2] = bn.Fp2_0;
        }
        if (arguments.length === 2) {
            /*
             * 
             */
            if (k instanceof BigInteger) {
                this.bn = bn;
                this.v = new Array(3);
                this.v[0] = new JSField2(bn, k);
                this.v[1] = this.v[2] = bn.Fp2_0;
            }
            /*
             * 
             */
            else if (k instanceof JSField2) {
                var v0 = k;
                this.bn = bn;
                this.v = new Array(3);
                this.v[0] = v0;
                this.v[1] = this.v[2] = bn.Fp2_0;
            }
            /*
             * 
             */
            else {
                var rand = k;
                this.bn = bn;
                this.v = new Array(3);
                this.v[0] = new JSField2(bn, rand);
                this.v[1] = new JSField2(bn, rand);
                this.v[2] = new JSField2(bn, rand);
            }
        }
        /*
         * 
         */
        if (arguments.length === 4) {
            var v0 = k;
            this.bn = bn;
            this.v = new Array(3);
            this.v[0] = v0;
            this.v[1] = v1;
            this.v[2] = v2;
        }
    };

    /*
     * 
     * @param {} rand
     * @returns 
     */
    JSField6.randomize = function (rand) {
        return new JSField6(this.bn, rand);
    };

    /*
     * 
     * @returns 
     */
    JSField6.isZero = function () {
        return this.v[0].isZero() && this.v[1].isZero() && this.v[2].isZero();
    };

    /*
     * 
     * @returns 
     */
    JSField6.isOne = function () {
        return this.v[0].isOne() && this.v[1].isZero() && this.v[2].isZero();
    };

    /*
     * 
     * @param {} o
     * @returns 
     */
    JSField6.equals = function (o) {
        if (!(o instanceof JSField6)) {
            return false;
        }
        var w = o;
        return this.bn === w.bn && // singleton comparison
            this.v[0].equals(w.v[0]) && this.v[1].equals(w.v[1]) && this.v[2].equals(w.v[2]);
    };

    /*
     * 
     * @returns 
     */
    JSField6.negate = function () {
        return new JSField6(this.bn, this.v[0].negate(), this.v[1].negate(), this.v[2].negate());
    };

    /*
     * Compute this^((p^2)^m), the m-th conjugate of this over GF(p^2).
     * @param {} m
     * @returns 
     */
    JSField6.conjugate = function (m) {
        switch (m) {
            case 0:
                return this;
            case 1:
                return new JSField6(this.bn, this.v[0], this.v[1].multiply(this.bn.zeta1).negate(), this.v[2].multiply(this.bn.zeta0));
            case 2:
                return new JSField6(this.bn, this.v[0], this.v[1].multiply(this.bn.zeta0), this.v[2].multiply(this.bn.zeta1).negate());
            default: // only to make the compiler happy
        }
    };

    /*
     * 
     * @param {} w
     * @returns 
     */
    JSField6.add = function (w) {
        if (this.bn !== w.bn) { // singleton comparison
            throw new Error(differentFields);
        }
        return new JSField6(this.bn, this.v[0].add(w.v[0]), this.v[1].add(w.v[1]), this.v[2].add(w.v[2]));
    };

    /*
     * 
     * @param {} w
     * @returns 
     */
    JSField6.subtract = function (w) {
        if (this.bn !== w.bn) { // singleton comparison
            throw new Error(differentFields);
        }
        return new JSField6(this.bn, this.v[0].subtract(w.v[0]), this.v[1].subtract(w.v[1]), this.v[2].subtract(w.v[2]));
    };

    /*
     * 
     * @param {} k
     * @returns 
     */
    JSField6.twice = function (k) {
        return new JSField6(this.bn, this.v[0].twice(k), this.v[1].twice(k), this.v[2].twice(k));
    };

    /*
     * 
     * @returns 
     */
    JSField6.halve = function () {
        return new JSField6(this.bn, this.v[0].halve(), this.v[1].halve(), this.v[2].halve());
    };

    /*
     * 
     * @param {} w
     * @returns 
     */
    JSField6.multiply = function (w) {
        if (w instanceof JSField6) {
            if (w === this) {
                return square();
            }
            if (this.bn !== w.bn) { // singleton comparison
                throw new Error(differentFields);
            }
            if (this.isOne() || w.isZero()) {
                return w;
            }
            if (this.isZero() || w.isOne()) {
                return this;
            }
            var d00 = this.v[0].multiply(w.v[0]);
            var d11 = this.v[1].multiply(w.v[1]);
            var d22 = this.v[2].multiply(w.v[2]);
            var d01 = this.v[0].add(this.v[1]).multiply(w.v[0].add(w.v[1])).subtract(d00.add(d11));
            var d02 = this.v[0].add(this.v[2]).multiply(w.v[0].add(w.v[2])).subtract(d00.add(d22));
            var d12 = this.v[1].add(this.v[2]).multiply(w.v[1].add(w.v[2])).subtract(d11.add(d22));
            if (this.bn.b === 3) {
                return new JSField6(this.bn, d12.divideV().add(d00), d22.divideV().add(d01), d02.add(d11));
            } else {
                return new JSField6(this.bn, d12.multiplyV().add(d00), d22.multiplyV().add(d01), d02.add(d11));
            }
        } else if (w instanceof JSField2) {
            if (this.bn !== w.bn) { // singleton comparison
                throw new Error(differentFields);
            }
            if (w.isOne()) {
                return this;
            }
            return new JSField6(this.bn, this.v[0].multiply(w), this.v[1].multiply(w), this.v[2].multiply(w));
        }
    };

    /*
     * 
     * @returns 
     */
    JSField6.multiplyConj = function () {
        if (this.isOne() || this.isZero()) {
            return this;
        }
        if (this.bn.b === 3) {
            return new JSField6(this.bn, this.v[0].square().subtract(this.v[1].multiply(this.v[2]).divideV()),
                this.v[2].square().divideV().subtract(this.v[0].multiply(this.v[1])).multiply(this.bn.zeta0),
                this.v[0].multiply(this.v[2]).subtract(this.v[1].square()).multiply(this.bn.zeta1));
        } else {
            // (v0^2 - v1*v2*xi) + (v2^2*xi - v0*v1)*zeta*w + (v0*v2 - v1^2)*(zeta+1)*w^2 =
            return new JSField6(this.bn, this.v[0].square().subtract(this.v[1].multiply(this.v[2]).multiplyV()),
                this.v[2].square().multiplyV().subtract(this.v[0].multiply(this.v[1])).multiply(this.bn.zeta0),
                this.v[0].multiply(this.v[2]).subtract(this.v[1].square()).multiply(this.bn.zeta1));
        }
    };

    /*
     * Complete the norm evaluation.
     * @param {} k
     * @returns 
     */
    JSField6.normCompletion = function (k) {
        var d00 = this.v[0].multiply(k.v[0]);
        var d12 = this.v[1].multiply(k.v[2]).add(this.v[2].multiply(k.v[1]));
        if (this.bn.b === 3) {
            return d12.divideV().add(d00);
        } else {
            return d12.multiplyV().add(d00);
        }
    };

    /*
     * 
     * @returns 
     */
    JSField6.square = function () {
        if (this.isZero() || this.isOne()) {
            return this;
        }
        var a0 = this.v[0];
        var a1 = this.v[1];
        var a2 = this.v[2];
        // Chung-Hasan SQR3 for F_{p^6} over F_{p^2}
        /*
        c0 = S0 = a0^2,
        S1 = (a2 + a1 + a0)^2,
        S2 = (a2 - a1 + a0)^2,
        c3 = S3 = 2*a1*a2,
        c4 = S4 = a2^2,
        T1 = (S1 + S2)/2,
        c1 = S1 - T1 - S3,
        c2 = T1 - S4 - S0.
        */
       var c0 = a0.square();
       var S1 = a2.add(a1).add(a0).square();
       var S2 = a2.subtract(a1).add(a0).square();
       var c3 = a1.multiply(a2).twice(1);
       var c4 = a2.square();
       var T1 = S1.add(S2).halve();
       var c1 = S1.subtract(T1).subtract(c3);
       var c2 = T1.subtract(c4).subtract(c0);
       // z^3 = xi
        // z^4 = xi*z
        // c4^*xi*z + c3*xi + c2*z^2 + c1*z + c0
        // c2^*z^2 + (c4^*xi + c1)*z + (c3*xi + c0)
        if (this.bn.b === 3) {
            c0 = c0.add(c3.divideV());
            c1 = c1.add(c4.divideV());
        } else {
            c0 = c0.add(c3.multiplyV());
            c1 = c1.add(c4.multiplyV());
        }
        return new JSField6(this.bn, c0, c1, c2);
    };

    /*
     * 
     */
    JSField6.multiplyV = function () {
        // (a0, a1, a2) -> (a2*xi, a0, a1)
        return new JSField6(this.bn, this.v[2].multiplyV(), this.v[0], this.v[1]);
    };

    /*
     * 
     */
    JSField6.divideV = function () {
        // (a0, a1, a2) -> (a2/xi, a0, a1)
        return new JSField6(this.bn, this.v[2].divideV(), this.v[0], this.v[1]);
    };

    JSField6.prototype.randomize = JSField6.randomize;
    JSField6.prototype.isZero = JSField6.isZero;
    JSField6.prototype.isOne = JSField6.isOne;
    JSField6.prototype.equals = JSField6.equals;
    JSField6.prototype.negate = JSField6.negate;
    JSField6.prototype.conjugate = JSField6.conjugate;
    JSField6.prototype.add = JSField6.add;
    JSField6.prototype.subtract = JSField6.subtract;
    JSField6.prototype.twice = JSField6.twice;
    JSField6.prototype.halve = JSField6.halve;
    JSField6.prototype.multiply = JSField6.multiply;
    JSField6.prototype.multiplyConj = JSField6.multiplyConj;
    JSField6.prototype.normCompletion = JSField6.normCompletion;
    JSField6.prototype.square = JSField6.square;
    JSField6.prototype.multiplyV = JSField6.multiplyV;
    JSField6.prototype.divideV = JSField6.divideV;

    return JSField6;
})();
