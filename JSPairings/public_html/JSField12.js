/*
 * JSField12.js
 * Arithmetic in the finite extension field GF(p^12) with p = 3 (mod 4).
 * Copyright (C) Jorge H. C. Lucema.
 */

var differentFields = "Operands are in different finite fields";

var JSField12 = (function () {

    var JSField12 = function (bn, k) {
        /*
         * 
         */
        if (arguments.length === 1) {
            var f = bn;
            this.bn = f.bn;
            this.v = new Array(6);
            for (var i = 0; i < 6; i++) {
                this.v[i] = f.v[i];
            }
        }
        if (arguments.length === 2) {
            /*
             * 
             */
            if (k instanceof BigInteger) {
                this.bn = bn;
                this.v = new Array(6);
                this.v[0] = new JSField2(bn, k); // caveat: no modular reduction!
                for (var i = 1; i < 6; i++) {
                    this.v[i] = new JSField2(bn);
                }
            }
            /*
             * 
             */
            else if (k instanceof Array) {
                this.bn = bn;
                this.v = k;
            }
            /*
             * 
             */
            else {
                var rand = k;
                this.bn = bn;
                this.v = new Array(6);
                for (var i = 0; i < 6; i++) {
                    this.v[i] = new JSField2(bn, rand);
                }
            }
        }
    };

    /*
     * 
     * @param {} rand
     * @returns 
     */
    JSField12.randomize = function (rand) {
        return new JSField12(this.bn, rand);
    };

    /*
     * 
     * @returns 
     */
    JSField12.isZero = function () {
        return this.v[0].isZero() && this.v[1].isZero() && this.v[2].isZero() &&
                this.v[3].isZero() && this.v[4].isZero() && this.v[5].isZero();
    };

    /*
     * 
     * @returns 
     */
    JSField12.isOne = function () {
        return this.v[0].isOne() && this.v[1].isZero() && this.v[2].isZero() &&
                this.v[3].isZero() && this.v[4].isZero() && this.v[5].isZero();
    };

    /*
     * 
     * @param {} o
     * @returns 
     */
    JSField12.equals = function (o) {
        if (!(o instanceof JSField12)) {
            return false;
        }
        var w = o;
        return this.bn === w.bn && // singleton comparison
            this.v[0].equals(w.v[0]) &&
            this.v[1].equals(w.v[1]) &&
            this.v[2].equals(w.v[2]) &&
            this.v[3].equals(w.v[3]) &&
            this.v[4].equals(w.v[4]) &&
            this.v[5].equals(w.v[5]);
    };

    /*
     * 
     * @returns 
     */
    JSField12.negate = function () {
        var w = new Array(6);
        for (var i = 0; i < 6; i++) {
            w[i] = this.v[i].negate();
        }
        return new JSField12(this.bn, w);
    };

    /*
     * 
     * @returns 
     */
    JSField12.frobenius = function () {
        var w = new Array(6);
        w[0] = this.v[0].conjugate();
        if (this.bn.b === 3) {
            w[1] = this.v[1].conjugate().multiplyV().multiply(this.bn.sigma);
            w[2] = this.v[2].conjugate().multiply(this.bn.zeta0).multiplyI().negate();
            w[3] = this.v[3].multiplyV().conjugate().multiply(this.bn.zeta0sigma);
            w[4] = this.v[4].conjugate().multiply(this.bn.zeta1);
            w[5] = this.v[5].conjugate().multiplyV().multiply(this.bn.zeta1sigma);
        } else {
            w[1] = this.v[1].multiplyV().conjugate().multiply(this.bn.zeta0sigma).negate();
            w[2] = this.v[2].conjugate().multiply(this.bn.zeta0).multiplyI();
            w[3] = this.v[3].conjugate().multiplyV().multiply(this.bn.zeta1sigma);
            w[4] = this.v[4].conjugate().multiply(this.bn.zeta1);
            w[5] = this.v[5].multiplyV().conjugate().multiply(this.bn.sigma);
        }
        return new JSField12(this.bn, w);
    };

    /*
     * Compute this^((p^2)^m), the m-th conjugate of this over GF(p^2).
     * @param {} m
     * @returns 
     */
    JSField12.conjugate = function (m) {
        switch (m) {
            case 0:
                return this;
            case 1:
                var w = new Array(6);
                w[0] = this.v[0];
                w[1] = this.v[1].multiply(this.bn.zeta0).negate();
                w[2] = this.v[2].multiply(this.bn.zeta1).negate();
                w[3] = this.v[3].negate();
                w[4] = this.v[4].multiply(this.bn.zeta0);
                w[5] = this.v[5].multiply(this.bn.zeta1);
                return new JSField12(this.bn, w);
            case 2:
                var w = new Array(6);
                w[0] = this.v[0];
                w[1] = this.v[1].multiply(this.bn.zeta1).negate();
                w[2] = this.v[2].multiply(this.bn.zeta0);
                w[3] = this.v[3];
                w[4] = this.v[4].multiply(this.bn.zeta1).negate();
                w[5] = this.v[5].multiply(this.bn.zeta0);
                return new JSField12(this.bn, w);
            case 3:
                var w = new Array(6);
                w[0] = this.v[0];
                w[1] = this.v[1].negate();
                w[2] = this.v[2];
                w[3] = this.v[3].negate();
                w[4] = this.v[4];
                w[5] = this.v[5].negate();
                return new JSField12(this.bn, w);
            case 4:
                var w = new Array(6);
                w[0] = this.v[0];
                w[1] = this.v[1].multiply(this.bn.zeta0);
                w[2] = this.v[2].multiply(this.bn.zeta1).negate();
                w[3] = this.v[3];
                w[4] = this.v[4].multiply(this.bn.zeta0);
                w[5] = this.v[5].multiply(this.bn.zeta1).negate();
                return new JSField12(this.bn, w);
            case 5:
                var w = new Array(6);
                w[0] = this.v[0];
                w[1] = this.v[1].multiply(this.bn.zeta1);
                w[2] = this.v[2].multiply(this.bn.zeta0);
                w[3] = this.v[3].negate();
                w[4] = this.v[4].multiply(this.bn.zeta1).negate();
                w[5] = this.v[5].multiply(this.bn.zeta0).negate();
                return new JSField12(this.bn, w);
            default: // only to make the compiler happy
        }
    };

    /*
     * 
     * @param {} k
     * @returns 
     */
    JSField12.add = function (k) {
        if (this.bn !== k.bn) { // singleton comparison
            throw new Error(differentFields);
        }
        var w = new Array(6);
        for (var i = 0; i < 6; i++) {
            w[i] = this.v[i].add(k.v[i]);
        }
        return new JSField12(this.bn, w);
    };

    /*
     * 
     * @param {} k
     * @returns 
     */
    JSField12.subtract = function (k) {
        if (this.bn !== k.bn) { // singleton comparison
            throw new Error(differentFields);
        }
        var w = new Array(6);
        for (var i = 0; i < 6; i++) {
            w[i] = this.v[i].subtract(k.v[i]);
        }
        return new JSField12(this.bn, w);
    };

    /*
     * 
     * @param {} k
     * @returns 
     */
    JSField12.multiply = function (k) {
        if (k instanceof JSField12) {
            if (k === this) {
                return this.square();
            }
            if (this.bn !== k.bn) { // singleton comparison
                throw new Error(differentFields);
            }
            if (this.isOne() || k.isZero()) {
                return k;
            }
            if (this.isZero() || k.isOne()) {
                return this;
            }
            JSField2.modoff();
            var w = new Array(6);
            if (k.v[2].isZero() && k.v[4].isZero() && k.v[5].isZero()) {
                if (this.v[2].isZero() && this.v[4].isZero() && this.v[5].isZero()) {
                    var d00 = this.v[0].multiply(k.v[0]);
                    var d11 = this.v[1].multiply(k.v[1]);
                    var d33 = this.v[3].multiply(k.v[3]);
                    var s01 = this.v[0].add(this.v[1]);
                    var t01 = k.v[0].add(k.v[1]);
                    var u01 = d00.add(d11);
                    var z01 = s01.multiply(t01);
                    var d01 = z01.subtract(u01);
                    var d13 = this.v[1].add(this.v[3]).multiply(k.v[1].add(k.v[3])).subtract(d11.add(d33));
                    u01 = u01.add(d01);
                    var d03 = s01.add(this.v[3]).multiply(t01.add(k.v[3])).subtract(u01.add(d33).add(d13));
                    var d05 = z01.subtract(u01);
                    if (this.bn.b === 3) {
                        w[0] = d33.divideV().add(d00);
                    } else {
                        w[0] = d33.multiplyV().add(d00);
                    }
                    w[1] = d01;
                    w[2] = d11;
                    w[3] = d03;
                    w[4] = d13;
                    w[5] = d05;
                } else {
                    var d00 = this.v[0].multiply(k.v[0]);
                    var d11 = this.v[1].multiply(k.v[1]);
                    var d33 = this.v[3].multiply(k.v[3]);
                    var s01 = this.v[0].add(this.v[1]);
                    var t01 = k.v[0].add(k.v[1]);
                    var u01 = d00.add(d11);
                    var d01 = s01.multiply(t01).subtract(u01);
                    var d02 = this.v[0].add(this.v[2]).multiply(k.v[0]).subtract(d00);
                    var d04 = this.v[0].add(this.v[4]).multiply(k.v[0]).subtract(d00);
                    var d13 = this.v[1].add(this.v[3]).multiply(k.v[1].add(k.v[3])).subtract(d11.add(d33));
                    var d15 = this.v[1].add(this.v[5]).multiply(k.v[1]).subtract(d11);
                    var s23 = this.v[2].add(this.v[3]);
                    var d23 = s23.multiply(k.v[3]).subtract(d33);
                    var d35 = this.v[3].add(this.v[5]).multiply(k.v[3]).subtract(d33);
                    u01 = u01.add(d01);
                    var u23 = d33.add(d23);
                    var d03 = s01.add(s23).multiply(t01.add(k.v[3])).subtract(u01.add(u23).add(d02).add(d13));
                    var s45 = this.v[4].add(this.v[5]);
                    var d05 = s01.add(s45).multiply(t01).subtract(u01.add(d04).add(d15));
                    var d25 = s23.add(s45).multiply(k.v[3]).subtract(u23.add(d35));
                    if (this.bn.b === 3) {
                        w[0] = d15.add(d33).divideV().add(d00);
                        w[1] = d25.divideV().add(d01);
                        w[2] = d35.divideV().add(d02).add(d11);
                    } else { // preferred representation:
                        w[0] = d15.add(d33).multiplyV().add(d00);
                        w[1] = d25.multiplyV().add(d01);
                        w[2] = d35.multiplyV().add(d02).add(d11);
                    }
                    w[3] = d03;
                    w[4] = d04.add(d13);
                    w[5] = d05.add(d23);
                }
            } else if (k.v[1].isZero() && k.v[4].isZero() && k.v[5].isZero()) {
                var d00 = this.v[0].multiply(k.v[0]);
                var d22 = this.v[2].multiply(k.v[2]);
                var d33 = this.v[3].multiply(k.v[3]);
                var s01 = this.v[0].add(this.v[1]);
                var d01 = s01.multiply(k.v[0]).subtract(d00);
                var d02 = this.v[0].add(this.v[2]).multiply(k.v[0].add(k.v[2])).subtract(d00.add(d22));
                var d04 = this.v[0].add(this.v[4]).multiply(k.v[0]).subtract(d00);
                var d13 = this.v[1].add(this.v[3]).multiply(k.v[3]).subtract(d33);
                var s23 = this.v[2].add(this.v[3]);
                var t23 = k.v[2].add(k.v[3]);
                var u23 = d22.add(d33);
                var d23 = s23.multiply(t23).subtract(u23);
                var d24 = this.v[2].add(this.v[4]).multiply(k.v[2]).subtract(d22);
                var d35 = this.v[3].add(this.v[5]).multiply(k.v[3]).subtract(d33);
                var u01 = d00.add(d01);
                var d03 = s01.add(s23).multiply(k.v[0].add(t23)).subtract(u01.add(u23).add(d02).add(d13).add(d23));
                var s45 = this.v[4].add(this.v[5]);
                var d05 = s01.add(s45).multiply(k.v[0]).subtract(u01.add(d04));
                var d25 = s23.add(s45).multiply(t23).subtract(u23.add(d23).add(d24).add(d35));
                if (this.bn.b === 3) {
                    w[0] = d24.add(d33).divideV().add(d00);
                    w[1] = d25.divideV().add(d01);
                    w[2] = d35.divideV().add(d02);
                } else { // preferred representation:
                    w[0] = d24.add(d33).multiplyV().add(d00);
                    w[1] = d25.multiplyV().add(d01);
                    w[2] = d35.multiplyV().add(d02);
                }
                w[3] = d03;
                w[4] = d04.add(d13).add(d22);
                w[5] = d05.add(d23);
            } else {
                var d00 = this.v[0].multiply(k.v[0]);
                var d11 = this.v[1].multiply(k.v[1]);
                var d22 = this.v[2].multiply(k.v[2]);
                var d33 = this.v[3].multiply(k.v[3]);
                var d44 = this.v[4].multiply(k.v[4]);
                var d55 = this.v[5].multiply(k.v[5]);
                var s01 = this.v[0].add(this.v[1]);
                var t01 = k.v[0].add(k.v[1]);
                var u01 = d00.add(d11);
                var d01 = s01.multiply(t01).subtract(u01);
                var d02 = this.v[0].add(this.v[2]).multiply(k.v[0].add(k.v[2])).subtract(d00.add(d22));
                var d04 = this.v[0].add(this.v[4]).multiply(k.v[0].add(k.v[4])).subtract(d00.add(d44));
                var d13 = this.v[1].add(this.v[3]).multiply(k.v[1].add(k.v[3])).subtract(d11.add(d33));
                var d15 = this.v[1].add(this.v[5]).multiply(k.v[1].add(k.v[5])).subtract(d11.add(d55));
                var  s23 = this.v[2].add(this.v[3]);
                var t23 = k.v[2].add(k.v[3]);
                var u23 = d22.add(d33);
                var d23 = s23.multiply(t23).subtract(u23);
                var d24 = this.v[2].add(this.v[4]).multiply(k.v[2].add(k.v[4])).subtract(d22.add(d44));
                var d35 = this.v[3].add(this.v[5]).multiply(k.v[3].add(k.v[5])).subtract(d33.add(d55));
                var s45 = this.v[4].add(this.v[5]);
                var t45 = k.v[4].add(k.v[5]);
                var u45 = d44.add(d55);
                var d45 = s45.multiply(t45).subtract(u45);
                u01 = u01.add(d01);
                u23 = u23.add(d23);
                u45 = u45.add(d45);
                var d03 = s01.add(s23).multiply(t01.add(t23)).subtract(u01.add(u23).add(d02).add(d13));
                var d05 = s01.add(s45).multiply(t01.add(t45)).subtract(u01.add(u45).add(d04).add(d15));
                var d25 = s23.add(s45).multiply(t23.add(t45)).subtract(u23.add(u45).add(d24).add(d35));
                if (this.bn.b === 3) {
                    w[0] = d15.add(d24).add(d33).divideV().add(d00); //w[0] =  c6/(1+i)+c0
                    w[1] = d25.divideV().add(d01);                   //w[1] =  c7/(1+i)+c1
                    w[2] = d35.add(d44).divideV().add(d02).add(d11); //w[2] =  c8/(1+i)+c2
                    w[3] = d45.divideV().add(d03);                   //w[3] =  c9/(1+i)+c3
                    w[4] = d55.divideV().add(d04).add(d13).add(d22); //w[4] = c10/(1+i)+c4
                    w[5] = d05.add(d23);                             //w[5] =  c5
                } else { // preferred representation:
                    w[0] = d15.add(d24).add(d33).multiplyV().add(d00);
                    w[1] = d25.multiplyV().add(d01);
                    w[2] = d35.add(d44).multiplyV().add(d02).add(d11);
                    w[3] = d45.multiplyV().add(d03);
                    w[4] = d55.multiplyV().add(d04).add(d13).add(d22);
                    w[5] = d05.add(d23);
                }
            }
            JSField2.modon();
            JSField2.offsetmod(6);
            return new JSField12(this.bn, w);
        } else if (k instanceof BigInteger) {
            var w = new Array(6);
            for (var i = 0; i < 6; i++) {
                w[i] = this.v[i].multiply(k);
            }
            return new JSField12(this.bn, w);
        } else if (k instanceof JSField2) {
            var w = new Array(6);
            for (var i = 0; i < 6; i++) {
                w[i] = this.v[i].multiply(k);
            }
            return new JSField12(this.bn, w);
        } else {
            var dr = new JSField6(this.bn, this.v[0], this.v[2], this.v[4]).multiply(k);
            var di = new JSField6(this.bn, this.v[1], this.v[3], this.v[5]).multiply(k);
            var m = new Array(6);
            m[0] = dr.v[0];
            m[1] = di.v[0];
            m[2] = dr.v[1];
            m[3] = di.v[1];
            m[4] = dr.v[2];
            m[5] = di.v[2];
            return new JSField12(this.bn, m);
        }
    };

    /*
     * 
     * @param {} h
     * @returns 
     */
    JSField12.decompress = function (h) {
        //3~s+3~m = 15 m
        if (!h.v[1].isZero()) {
            if (this.bn.b === 2) {
                h.v[3] = h.v[5].square().multiplyV().add(h.v[2].square().multiply(new Number(3))).subtract(h.v[4].twice(1)).multiply(h.v[1].twice(2).inverse());
                h.v[0] = h.v[3].square().twice(1).add(h.v[1].multiply(h.v[5])).subtract(h.v[4].multiply(h.v[2]).multiply(new Number(3))).multiplyV().add(BigInteger.ONE);
                
            } else {
                h.v[3] = h.v[5].square().divideV().add(h.v[2].square().multiply(new Number(3))).subtract(h.v[4].twice(1)).multiply(h.v[1].twice(2).inverse());
                h.v[0] = h.v[3].square().twice(1).add(h.v[1].multiply(h.v[5])).subtract(h.v[4].multiply(h.v[2]).multiply(new Number(3))).divideV().add(BigInteger.ONE);
            }
        } else {
            h.v[3] = h.v[2].multiply(h.v[5]).twice(1).multiply(h.v[4].inverse());
            h.v[0] = h.v[3].square().twice(1).subtract(h.v[4].multiply(h.v[2]).multiply(new Number(3))).multiplyV().add(BigInteger.ONE);
        }
    };

    /*
     * 
     * @returns 
     */
    JSField12.square = function () {
        if (this.isZero() || this.isOne()) {
            return this;
        }
        JSField2.modoff();
        // 4sF_p^4 + 1mF_p^4 = 4*(3sF_p^2)+3mF_p^2 = 30mF_p
        // Chung-Hasan SQR3
        var a0 =  new JSField4(this.bn, this.v[0], this.v[3]);
        var a1 =  new JSField4(this.bn, this.v[1], this.v[4]);
        var a2 =  new JSField4(this.bn, this.v[2], this.v[5]);
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

        JSField2.modon();
        JSField2.offsetmod(6);
        var v = new Array(6);
        v[0] = c0.re;
        v[1] = c1.re;
        v[2] = c2.re;
        v[3] = c0.im;
        v[4] = c1.im;
        v[5] = c2.im;
        return new JSField12(this.bn, v);
    };

    /*
     * 
     * @returns 
     */
    JSField12.compressedSquare = function () {
        // Karabina's technique to square a element of Cyclotomic Subgroup
        var h = new JSField12(this.bn.Fp12_0);
        if (this.bn.b === 2) {
            var A23 = this.v[1].add(this.v[4]).multiply(this.v[1].add(this.v[4].multiplyV()));
            var A45 = this.v[2].add(this.v[5]).multiply(this.v[2].add(this.v[5].multiplyV()));
            var B45 = this.v[2].multiply(this.v[5]);
            var B23 = this.v[1].multiply(this.v[4]);
            h.v[1] = this.v[1].add(B45.multiplyV().multiply(new Number(3))).twice(1);
            h.v[4] = A45.subtract(B45.add(B45.multiplyV())).multiply(new Number(3)).subtract(this.v[4].twice(1));
            h.v[2] = A23.subtract(B23.add(B23.multiplyV())).multiply(new Number(3)).subtract(this.v[2].twice(1));
            h.v[5] = this.v[5].add(B23.multiply(new Number(3))).twice(1);
        } else {
            var A23 = this.v[1].add(this.v[4]).multiply(this.v[1].add(this.v[4].divideV()));
            var A45 = this.v[2].add(this.v[5]).multiply(this.v[2].add(this.v[5].divideV()));
            var B45 = this.v[2].multiply(this.v[5]);
            var B23 = this.v[1].multiply(this.v[4]);
            h.v[1] = this.v[1].add(B45.divideV().multiply(new Number(3))).twice(1);
            h.v[4] = A45.subtract(B45.add(B45.divideV())).multiply(new Number(3)).subtract(this.v[4].twice(1));
            h.v[2] = A23.subtract(B23.add(B23.divideV())).multiply(new Number(3)).subtract(this.v[2].twice(1));
            h.v[5] = this.v[5].add(B23.multiply(new Number(3))).twice(1);
        }
        //decompress(h);
        return h;
    };

    /*
     * 
     * @returns 
     */
    JSField12.uniSquare = function () {
        JSField2.modoff();
        //* 18 mFp
        // Granger/Scott technique to square a element of Cyclotomic Subgroup
        var a0sqr = this.v[0].square();
        var a1sqr = this.v[3].square();
        var b0sqr = this.v[1].square();
        var b1sqr = this.v[4].square();
        var c0sqr = this.v[2].square();
        var c1sqr = this.v[5].square();
        var a0, a1, b0, b1, c0, c1;
        if (this.bn.b === 3) {
            //a0 = 3*(a0^2 + V^{-1}*a1^2) - 2*a0
            a0 = a1sqr.divideV().add(a0sqr).multiply(new Number(3)).subtract(this.v[0].twice(1));
            //a1 = 3*[(a0+a1)^2 - a0^2 - a1^2] + 2*a1
            a1 = this.v[0].add(this.v[3]).square().subtract(a0sqr).subtract(a1sqr).multiply(new Number(3)).add(this.v[3].twice(1));
            //b0 = 3V*[(c0 + c1)^2 - c_0^2 - c1^2] + 2*b0
            b0 = this.v[2].add(this.v[5]).square().subtract(c0sqr).subtract(c1sqr).multiply(new Number(3)).divideV().add(this.v[1].twice(1));
            //b1 = 3*(c0^2 + V^{-1}*c1^2) - 2*b1
            b1 = c0sqr.add(c1sqr.divideV()).multiply(new Number(3)).subtract(this.v[4].twice(1));
            //c0 = 3*(b0^2 + V^{-1}*b1^2) - 2*c0
            c0 = b1sqr.divideV().add(b0sqr).multiply(new Number(3)).subtract(this.v[2].twice(1));
            //c1 = 3*[(b0 + b1)^2 - b0^2 - b1^2] + 2*c1
            c1 = this.v[1].add(this.v[4]).square().subtract(b0sqr).subtract(b1sqr).multiply(new Number(3)).add(this.v[5].twice(1));
        } else {
            //A = 3a^2-2conj(a)
            //a0 = 3*(a0^2 + V*a1^2) - 2*a0
            a0 = a1sqr.multiplyV().add(a0sqr).multiply(new Number(3)).subtract(this.v[0].twice(1));
            //a1 = 3*[(a0+a1)^2 - a0^2 - a1^2] + 2*a1
            a1 = this.v[0].add(this.v[3]).square().subtract(a0sqr).subtract(a1sqr).multiply(new Number(3)).add(this.v[3].twice(1));
            //B = 3V*c^2-2conj(b)
            //b0 = 3V*[(c0 + c1)^2 - c_0^2 - c1^2] + 2*b0
            b0 = this.v[2].add(this.v[5]).square().subtract(c0sqr).subtract(c1sqr).multiply(new Number(3)).multiplyV().add(this.v[1].twice(1));
            //b1 = 3*(c0^2 + V*c1^2) - 2*b1
            b1 = c0sqr.add(c1sqr.multiplyV()).multiply(new Number(3)).subtract(this.v[4].twice(1));
            //C = 3b^2 - 2conj(c)
            //c0 = 3*(b0^2 + V*b1^2) - 2*c0
            c0 = b1sqr.multiplyV().add(b0sqr).multiply(new Number(3)).subtract(this.v[2].twice(1));
            //c1 = 3*[(b0 + b1)^2 - b0^2 - b1^2] + 2*c1
            c1 = this.v[1].add(this.v[4]).square().subtract(b0sqr).subtract(b1sqr).multiply(new Number(3)).add(this.v[5].twice(1));
        }
        JSField2.modon();
        JSField2.offsetmod(6);
        var m = new Array(6);
        m[0] = a0;
        m[1] = b0;
        m[2] = c0;
        m[3] = a1;
        m[4] = b1;
        m[5] = c1;
        return new JSField12(this.bn, m);
    };

    /*
     * 
     * @returns 
     */
    JSField12.multiplyV = function () {
        // (a0, a1, a2, a3, a4, a5) -> (a4*xi, a5*xi, a0, a1, a2, a3)
        var m = new Array(6);
        m[0] = this.v[4].multiplyV();
        m[1] = this.v[5].multiplyV();
        m[2] = this.v[0];
        m[3] = this.v[1];
        m[4] = this.v[2];
        m[5] = this.v[3];
        return new JSField12(this.bn, m);
    };

    /*
     * 
     * @returns 
     */
    JSField12.divideV = function () {
        // (a0, a1, a2, a3, a4, a5) -> (a4/xi, a5/xi, a0, a1, a2, a3)
        var m = new Array(6);
        m[0] = this.v[4].divideV();
        m[1] = this.v[5].divideV();
        m[2] = this.v[0];
        m[3] = this.v[1];
        m[4] = this.v[2];
        m[5] = this.v[3];
        return new JSField12(this.bn, m);
    };

    /*
     * 
     * @returns 
     */
    JSField12.norm6 = function () {
        // (a + bz)(a - bz) = a^2 - b^2*z^2 = a^2 - b^2*xi
        var re = new JSField6(this.bn, this.v[0], this.v[2], this.v[4]);
        var im = new JSField6(this.bn, this.v[1], this.v[3], this.v[5]);
        if (this.bn.b === 3) {
            return re.square().subtract(im.square().divideV());
        } else {
            return re.square().subtract(im.square().multiplyV());
        }
    };

    /*
     * 
     * @returns 
     */
    JSField12.inverse = function () {
        JSField2.modoff();

        // Total equiv: 22+15+9+9+36 = 91 mFp
        var l = this.norm6();                       //l = f^{1+q^3}                 - 22 mFp
        var m = l.multiplyConj().conjugate(1); //m = (l^q)*(l^{q^2}) in F_q^3  - 15 mFp
        var e = l.normCompletion(m);           //e = l*m                       -  9 mFp
        var d = m.multiply(e.inverse());       //d = m*e^{-1}                  -  9 mFp
        var c = this.conjugate(3).multiply(d);     //c = d*f^{q^3}                 - 36 mFp

        JSField2.modon();
        JSField2.offsetmod(6);

        return c;
    };

    /*
     * 
     * @param {} k
     * @returns 
     */
    JSField12.plainExp = function (k) {
        /*
         * This method is likely to be very fast, because k is very sparse
         */
        var w = this;
        for (var i = k.bitLength()-2; i >= 0; i--) {
            w = w.square();
            if (k.testBit(i)) {
                w = w.multiply(this);
            }
        }
        return w;
    };

    /*
     * 
     */
    JSField12.uniExp = function (k) {
        var w = new JSField12(this);
        for (var i = k.bitLength()-2; i >= 0; i--) {
            w = w.compressedSquare();//uniSquare();
            if (k.testBit(i)) {
                this.decompress(w);
                w = w.multiply(this);
            }
        }

        return w;
    };

    /*
     * 
     */
    JSField12.finExp = function () {
        var f = this;
        // Compute the easy part
        f = f.conjugate(3).multiply(f.inverse()); // f = f^(p^6 - 1)
        f = f.conjugate(1).multiply(f); // f = f^(p^2 + 1)
        var fconj = f.conjugate(3);
        if (this.bn.u.signum() >= 0) {
            var fu  = fconj.uniExp(this.bn.u);            //fu  = f^{p^6*u}
            var fu2 = fu.conjugate(3).uniExp(this.bn.u);  //fu2 = f^{p^12*u^2}
            var fu3 = fu2.conjugate(3).uniExp(this.bn.u); //fu3 = f^{p^18*u^3}
        } else {
            var fu = f.uniExp(this.bn.u.negate());        //fu  = f^{-u}
            var fu2 = fu.uniExp(this.bn.u.negate());      //fu2 = f^{u^2}
            var fu3 = fu2.uniExp(this.bn.u.negate());     //fu3 = f^{-u^3}
        }

        var fp = f.frobenius();
        var fp2 = fp.frobenius();
        var fp3 = fp2.frobenius();

        var fup = fu.frobenius();
        var fu2p = fu2.frobenius();
        var fu3p = fu3.frobenius();
        var fu2p2 = fu2.conjugate(1);

        var y0 = fp.multiply(fp2).multiply(fp3);
        var y1 = fconj;
        var y2 = fu2p2;
        var y3 = fup;
        var y4 = fu.multiply(fu2p.conjugate(3));
        var y5 = fu2.conjugate(3);
        var y6 = fu3.multiply(fu3p);


        var T0 = y6.uniSquare().multiply(y4).multiply(y5);
        var T1 = y3.multiply(y5).multiply(T0).uniSquare();
        T0 = T0.multiply(y2);
        T1 = T1.multiply(T0).uniSquare();
        T0 = T1.multiply(y1).uniSquare();
        T1 = T1.multiply(y0);
        T0 = T0.multiply(T1);
        f = T0;

        return f;
    };

    /*
     * 
     */
    JSField12.exp = function (k) {
        return this.plainExp(k);
    };

    JSField12.prototype.randomize = JSField12.randomize;
    JSField12.prototype.isZero = JSField12.isZero;
    JSField12.prototype.isOne = JSField12.isOne;
    JSField12.prototype.equals = JSField12.equals;
    JSField12.prototype.negate = JSField12.negate;
    JSField12.prototype.frobenius = JSField12.frobenius;
    JSField12.prototype.conjugate = JSField12.conjugate;
    JSField12.prototype.add = JSField12.add;
    JSField12.prototype.subtract = JSField12.subtract;
    JSField12.prototype.multiply = JSField12.multiply;
    JSField12.prototype.decompress = JSField12.decompress;
    JSField12.prototype.square = JSField12.square;
    JSField12.prototype.compressedSquare = JSField12.compressedSquare;
    JSField12.prototype.uniSquare = JSField12.uniSquare;
    JSField12.prototype.multiplyV = JSField12.multiplyV;
    JSField12.prototype.divideV = JSField12.divideV;
    JSField12.prototype.norm6 = JSField12.norm6;
    JSField12.prototype.inverse = JSField12.inverse;
    JSField12.prototype.plainExp = JSField12.plainExp;
    JSField12.prototype.uniExp = JSField12.uniExp;
    JSField12.prototype.finExp = JSField12.finExp;
    JSField12.prototype.exp = JSField12.exp;

    return JSField12;
})();
