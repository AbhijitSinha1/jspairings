/*
 * JSParams.js
 * Parameters for Barreto-Naehrig (BN) pairing-friendly elliptic curves.
 * Copyright (C) Jorge H. C. Lucema.
 */

var invalidParams = "The specified parameters do not properly define a suitable BN curve";

var JSParams = (function () {

    var JSParams = function (fieldBits) {

        if (arguments.length === 1) {
            if (typeof fieldBits !== "number") {
                throw new Error(invalidParams + ": Parameter is not number");
            }

            /*
             * Convenient BigInteger constants
             */
            this._2 = new BigInteger("2");
            this._3 = new BigInteger("3");
            this._4 = new BigInteger("4");
            this._5 = new BigInteger("5");
            this._6 = new BigInteger("6");
            this._7 = new BigInteger("7");
            this._8 = new BigInteger("8");
            this._9 = new BigInteger("9");
            this._24 = new BigInteger("24");
            this._0 = BigInteger.ZERO;
            this._1 = BigInteger.ONE;

            this.m = fieldBits;
            this.b = 3; // default assumption; corrected below on demand
            switch (fieldBits) {
                case 34:
                    this.b = 2; this.u = new BigInteger("-10000101", 2); // Hamming weight 3
                    break;
                case 48:
                    this.u = new BigInteger("-11001011001", 2); // Hamming weight 6
                    break;
                case 56:
                    this.u = new BigInteger("1011001111011", 2); // Hamming weight 9
                    break;
                case 64:
                    this.u = new BigInteger("110010000111111", 2); // Hamming weight 9
                    break;
                case 72:
                    this.u = new BigInteger("10110000111001011", 2); // Hamming weight 9
                    break;
                case 80:
                    this.u = new BigInteger("1101000010001011011", 2); // Hamming weight 9
                    break;
                case 88:
                    this.u = new BigInteger("-110000011000001110001", 2); // Hamming weight 8
                    break;
                case 96:
                    this.u = new BigInteger("11010000000000000010111", 2); // Hamming weight 7
                    break;
                case 104:
                    this.u = new BigInteger("1101000000000000000100011", 2); // Hamming weight 6
                    break;
                case 112:
                    this.u = new BigInteger("-110000001100001000000000001", 2); // Hamming weight 6
                    break;
                case 120:
                    this.u = new BigInteger("11000000100000000100100000011", 2); // Hamming weight 7
                    break;
                case 128:
                    this.u = new BigInteger("-1100111000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 136:
                    this.u = new BigInteger("-110000100000000000000001100000001", 2); // Hamming weight 6
                    break;
                case 144:
                    this.u = new BigInteger("-10110010000000010000000000000000001", 2); // Hamming weight 6
                    break;
                case 152:
                    this.u = new BigInteger("-1100100001000000100000000000000000001", 2); // Hamming weight 6
                    break;
                case 158: ////////////////////////////////
                    this.b = 2; this.u = new BigInteger("100000000000000100000000000000000100011", 2); // Hamming weight 5
                    //u = new BigInteger("-100000010000000100000000010000000000001", 2); // Hamming weight 5
                    break;
                case 160:
                    //u = new BigInteger("110100001000000000000100010000000000011", 2); // *** ISO, Hamming weight 8
                    this.u = new BigInteger("-110010001000000010000000000000001000001", 2); // Hamming weight 7
                    break;
                case 168:
                    this.u = new BigInteger("-11001000000000000000000010000001000000001", 2); // Hamming weight 6
                    break;
                case 176:
                    this.u = new BigInteger("-1100100100000000000000000010000000000000001", 2); // Hamming weight 6
                    break;
                case 184:
                    this.u = new BigInteger("-110000001100000000000000000000001000000000001", 2); // Hamming weight 6
                    break;
                case 190: ////////////////////////////////
                    this.b = 2; this.u = new BigInteger("-10000000010000100100000000000000000000000000001", 2); // Hamming weight 5
                    //u = new BigInteger("10000000010000000000000000000000000000001000011", 2); // Hamming weight 5
                    break;
                case 192:
                    //u = new BigInteger("11000000000000000001000000000000000010000010011", 2); // *** ISO, Hamming weight 7
                    this.u = new BigInteger("-11000000000000000000010010000000000010000000001", 2); // Hamming weight 6
                    break;
                case 200:
                    this.u = new BigInteger("-1101000000000000000000001000000000000010000000001", 2); // Hamming weight 6
                    break;
                case 208:
                    this.u = new BigInteger("110000000000000000000000000000000000000000100000011", 2); // Hamming weight 5
                    break;
                case 216:
                    this.u = new BigInteger("-11000000000000000010000000000000000000000000000000001", 2); // Hamming weight 4
                    break;
                case 222: ////////////////////////////////
                    this.b = 2; this.u = new BigInteger("1000010000000000010000000000000000000000000000000000011", 2); // Hamming weight 5
                    //u = new BigInteger("1000000000000000000000000000000000000000100100000000011", 2); // Hamming weight 5
                    break;
                case 224:
                    //u = new BigInteger("1100000000000000000000100000001000000000000001000000011", 2); // *** ISO, Hamming weight 7
                    this.u = new BigInteger("-1100000100000000000000000010000000100000000000000000001", 2); // Hamming weight 6
                    break;
                case 232:
                    this.u = new BigInteger("-110000000100000000100000000000000000000000000010000000001", 2); // Hamming weight 6
                    break;
                case 240:
                    this.u = new BigInteger("-11000100000000000000000000000010000000000000000000100000001", 2); // Hamming weight 6
                    break;
                case 248:
                    this.u = new BigInteger("-1100010000001000000000100000000000000000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 254: ////////////////////////////////
                    this.b = 2; this.u = new BigInteger("-100000010000000000000000000000000000000000000000000000000000001", 2); // Hamming weight 3
                    //u = new BigInteger("-100000010000000000000000000000000000000001000000000001000000001", 2); // Hamming weight 5
                    break;
                case 256:
                    this.u = new BigInteger("110000010000000000000000000000000000000000001000000000001000011", 2); // *** ISO, Hamming weight 7
                    //u = new BigInteger("-110000100000100000000001000000000000000000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 264:
                    this.u = new BigInteger("11000000000000000001000000000000000000000000000000000100000000011", 2); // Hamming weight 6
                    break;
                case 272:
                    this.u = new BigInteger("1100000100000000000010000000000000000000000000000000000000001000011", 2); // Hamming weight 7
                    break;
                case 280:
                    this.u = new BigInteger("-110001000000000000000000000000100000100000000000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 288:
                    this.u = new BigInteger("11000000000000000000000000000000000100000001000000000000000000000000011", 2); // Hamming weight 6
                    break;
                case 296:
                    this.u = new BigInteger("-1100000000000100000000000000100000000000000000000000000000000000000010001", 2); // Hamming weight 6
                    break;
                case 304:
                    this.u = new BigInteger("110000000000000100000000000000000000000000000000000000000001000000000000011", 2); // Hamming weight 6
                    break;
                case 312:
                    this.u = new BigInteger("-11000000000000001000000000000000000000000001000010000000000000000000000000001", 2); // Hamming weight 6
                    break;

                case 318: ////////////////////////////////
                    this.b = 2; this.u = new BigInteger("1000000000000000100000000000000000000000000000000000000000000000000000000000011", 2); // Hamming weight 4
                    //u = new BigInteger("-1000000000000000100000000000000000000000000000000000000000000000001000000010001", 2); // Hamming weight 5
                    break;
                case 320:
                    this.u = new BigInteger("-1100000001000000000000000000000000000000000010001000000000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 328:
                    this.u = new BigInteger("-110000000000100000000000000000000000000000000000000010000000000100000000000000001", 2); // Hamming weight 6
                    break;
                case 336:
                    this.u = new BigInteger("-11000000000000000000000000000000000000010000000000000000000100000000000000000000001", 2); // Hamming weight 5
                    break;
                case 344:
                    this.u = new BigInteger("-1100100000000000000000000000000000000000010000001000000000000000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 352:
                    this.u = new BigInteger("-110000100000000000000000000000000000001000000000000000000000010000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 360:
                    this.u = new BigInteger("-11000100001000000000000000000000000000000000000000001000000000000000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 368:
                    this.u = new BigInteger("-1100000000000000001010000000000000000000000001000000000000000000000000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 376:
                    this.u = new BigInteger("-110000000010000000000100000000000000000000000000000000000000000000000000000000000000100000001", 2); // Hamming weight 6
                    break;
                case 382: ////////////////////////////////
                    this.b = 2; this.u = new BigInteger("-10000000000000000010001000000000000000000000000000000000000000000000000000000000000000000000001", 2); // Hamming weight 4
                    break;
                case 384:
                    //u = new BigInteger("11001000000000000010000000000000000000000000000000000000000000000100000000000000000000000000011", 2); // *** ISO, Hamming weight 7
                    this.u = new BigInteger("-11000000000000000000000000000000000001000000000000000000000000000000000000000001000000000000001", 2); // Hamming weight 5
                    break;
                case 392:
                    this.u = new BigInteger("-1100100001000000000000000000000000000000000000000000000000000000000000000000100000000000000000001", 2); // Hamming weight 6
                    break;
                case 400:
                    this.u = new BigInteger("110000000000000000000000000000000000000000000000000000000000000000000000000000000100000001000000011", 2); // Hamming weight 6
                    break;
                case 408:
                    this.u = new BigInteger("-11000000000000010000000000000000000000000000000000000000000000000000000000000000000010000000000010001", 2); // Hamming weight 6
                    break;
                case 416:
                    this.u = new BigInteger("-1100100000000000000000000000000000010000000000000000000000000000001000000000000000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 424:
                    this.u = new BigInteger("-110000000000000000000000000000100000000000010000000000000000000000001000000000000000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 432:
                    this.u = new BigInteger("-11000000000000000000000000000000000000000000000000010010000000000000000000000000000000000010000000000000001", 2); // Hamming weight 6
                    break;
                case 440:
                    this.u = new BigInteger("-1100100000000000000000000000000000000000001000000000000000000000000000010000000000000000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 446: ////////////////////////////////
                    this.b = 2; this.u = new BigInteger("-100000000000000000000000000000000000000000000010000000001000000000000000000000000000000000000000000000000000001", 2); // Hamming weight 4
                    break;
                case 448:
                    this.u = new BigInteger("110000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000001000000000000011", 2); // Hamming weight 6
                    break;
                case 456:
                    this.u = new BigInteger("-11000000000000000000000000000100000000000000000000000000000000000000000000010000000000000000000000000000000000001", 2); // Hamming weight 5
                    break;
                case 464:
                    this.u = new BigInteger("-1100100000000000000000000000000000011000000000000000000000000000000000000000000000000000000000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 472:
                    this.u = new BigInteger("-110000001000000000000000000000000000000000000000000000000000000000000000010000000000000000000000100000000000000000001", 2); // Hamming weight 6
                    break;
                case 480:
                    this.u = new BigInteger("-11000000000000000100000000010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001", 2); // Hamming weight 5
                    break;
                case 488:
                    this.u = new BigInteger("-1100000001000000000000000000000000000000000000000010000000000000000000000000000000001000000000000000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 496:
                    this.u = new BigInteger("-110010000000000000000000000100000000000000000000000000000000000000100000000000000000000000000000000000000000000000000000001", 2); // Hamming weight 6
                    break;
                case 504:
                    this.u = new BigInteger("-11010000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000001", 2); // Hamming weight 5
                    break;
                case 512:
                    //u = new BigInteger("1100001000000000000000000000000000000000000000000000000000000000000000000000000000010000000000000000000000001000000000000000011", 2); // *** ISO, Hamming weight 7
                    this.u = new BigInteger("-1100001000000000000000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000010000000001", 2); // Hamming weight 6
                    break;
                default:
                    throw new Error(invalidParams + ": Field size in bits must be a multiple of 8 between 48 and 512");

            }
            //p = 36*u^4 + 36*u^3 + 24*u^2 + 6*u + 1 = (((u + 1)*6*u + 4)*u + 1)*6*u + 1
            this.p = this.u.add(this._1).multiply(this._6.multiply(this.u)).add(this._4).multiply(this.u).add(this._1).multiply(this._6.multiply(this.u)).add(this._1);
            this.t = this._6.multiply(this.u).multiply(this.u).add(this._1); // 6*u^2 + 1
            this.ht = this.p.subtract(this._1).add(this.t);
            //n = 36*u^4 + 36*u^3 + 18*u^2 + 6*u + 1
            this.n = this.p.add(this._1).subtract(this.t);
            //zeta = 18*u^3 + 18*u^2 + 9*u + 1;
            this.zeta = this._9.multiply(this.u).multiply(this.u.shiftLeft(1).multiply(this.u.add(this._1)).add(this._1)).add(this._1);
            this.zeta0 = this.zeta;
            this.zeta1 = this.zeta.add(this._1);
            //rho = |36*u^3 + 18*u^2 + 6*u + 1| = |6*u*(3*u*(2*u + 1) + 1) + 1)|;
            this.rho = this._6.multiply(this.u).multiply(this._3.multiply(this.u).multiply(this.u.shiftLeft(1).add(this._1)).add(this._1)).add(this._1);
            //*
            if (this.rho.signum() < 0) {
                this.rho = this.rho.negate();
            }
            //*/
            this.optOrd = this._6.multiply(this.u).add(this._2);
            //*
            if (this.optOrd.signum() < 0) {
                this.optOrd = this.optOrd.negate();
            }
            //*/
            this.sqrtExponent = this.p.add(this._1).shiftRight(2); // (p + 1)/4
            this.cbrtExponent = this.p.add(this.p).add(this._1).divide(this._9); // (2*p + 1)/9
            this.sqrtExponent2 = this.p.multiply(this.p).add(this._7).shiftRight(4); // (p^2 + 7)/16
            this.cbrtExponent2 = this.p.multiply(this.p).add(this._2).divide(this._9); // (p^2 + 2)/9
            this.sigma = this.p.subtract(this._4).modPow(this.p.subtract(this._1).subtract(this.p.add(this._5).divide(this._24)), this.p); // (-1/4)^((p+5)/24)
            this.zeta0sigma = this.zeta0.multiply(this.sigma).mod(this.p);
            this.zeta1sigma = this.zeta1.multiply(this.sigma).mod(this.p);

            //  2*sigma^2 = -zeta
            // -2*sigma^3 = -sigma*(2*sigma^2) = sigma*zeta
            // -4*sigma^4 = zeta+1
            // -4*sigma^5 = -4*sigma^4*sigma = (zeta+1)*sigma
            //  8*sigma^6 = -1
            //-16*sigma^9 = -2*sigma^3*(8*sigma^6) = 2*sigma^3
            this.invSqrtMinus2 = this.p.subtract(this._2).modPow(this.p.subtract(this._1).subtract(this.p.add(this._1).shiftRight(2)), this.p); // 1/sqrt(-2) = (-2)^{-(p+1)/4}
            this.sqrtI = new JSField2(this, this.invSqrtMinus2, (this.invSqrtMinus2.signum() !== 0) ? this.p.subtract(this.invSqrtMinus2) : this.invSqrtMinus2, false); // sqrt(i) = (1 - i)/sqrt(-2)
            this.Fp2_0 = new JSField2(this, this._0);
            this.Fp2_1 = new JSField2(this, this._1);
            this.Fp2_i = new JSField2(this, this._0, this._1, false);
            this.Fp12_0 = new JSField12(this, this._0);
            this.Fp12_1 = new JSField12(this, this._1);

            this.latInv = new Array(4);
            this.latInv[0] = this.u.shiftLeft(1).add(this._3).multiply(this.u).add(this._1); // 2*u^2 + 3*u + 1 = (2*u + 3)*u + 1
            this.latInv[1] = this.u.multiply(this._3).add(this._2).multiply(this.u).multiply(this.u).shiftLeft(2).add(this.u); // 12*u^3 + 8*u^2 + u = ((3*u + 2)*u)*4*u + u
            this.latInv[2] = this.u.multiply(this._3).add(this._2).multiply(this.u).multiply(this.u).shiftLeft(1).add(this.u); //  6*u^3 + 4*u^2 + u = ((3*u + 2)*u)*2*u + u
            this.latInv[3] = this.u.multiply(this.u).shiftLeft(1).add(this.u).negate(); // -(2*u^2 + u)

            this.latRed = new Array(4);
            for (var i = 0; i < 4; i++) {
                this.latRed[i] = new Array(4);
            }
            /*
                u+1,       u,      u,  -2*u,
                2*u+1,    -u, -(u+1),    -u,
                2*u,   2*u+1,  2*u+1, 2*u+1,
                u-1,   4*u+2, -2*u+1,   u-1
             */

            this.latRed[0][0] = this.u.add(this._1);
            this.latRed[0][1] = this.u;
            this.latRed[0][2] = this.u;
            this.latRed[0][3] = this.u.shiftLeft(1).negate();

            this.latRed[1][0] = this.u.shiftLeft(1).add(this._1);
            this.latRed[1][1] = this.u.negate();
            this.latRed[1][2] = this.u.add(this._1).negate();
            this.latRed[1][3] = this.u.negate();

            this.latRed[2][0] = this.u.shiftLeft(1);
            this.latRed[2][1] = this.u.shiftLeft(1).add(this._1);
            this.latRed[2][2] = this.u.shiftLeft(1).add(this._1);
            this.latRed[2][3] = this.u.shiftLeft(1).add(this._1);

            this.latRed[3][0] = this.u.subtract(this._1);
            this.latRed[3][1] = this.u.shiftLeft(2).add(this._2);
            this.latRed[3][2] = this.u.shiftLeft(1).negate().add(this._1);
            this.latRed[3][3] = this.u.subtract(this._1);
        }
    };

    /*
     * Compute the quadratic character of v, i.e. (v/p) for prime p.
     * @param {BigInteger} v
     * @returns J
     */
    JSParams.legendre = function (v) {
        //return v.modPow(p.shiftRight(1), p).add(this._1).compareTo(p) == 0 ? -1 : 1; // v^((p-1)/2) mod p = (v/p) for prime p
        var J = 1;
        var x = new BigInteger(v);
        var y = this.p;
        if (x.signum() < 0) {
            x = x.negate();
            if (y.testBit(0) && y.testBit(1)) { // y = 3 (mod 4)
                J = -J;
            }
        }
        while (y.compareTo(this._1) > 0) {
            x = x.mod(y);
            if (x.compareTo(y.shiftRight(1)) > 0) {
                x = y.subtract(x);
                if (y.testBit(0) && y.testBit(1)) { // y = 3 (mod 4)
                    J = -J;
                }
            }
            if (x.signum() === 0) {
                x = this._1;
                y = this._0;
                J = 0;
                break;
            }
            while (!x.testBit(0) && !x.testBit(1)) { // 4 divides x
                x = x.shiftRight(2);
            }
            if (!x.testBit(0)) { // 2 divides x
                x = x.shiftRight(1);
                if (y.testBit(0) && (y.testBit(1) === !y.testBit(2))) { // y = �3 (mod 8)
                    J = -J;
                }
            }
            if (x.testBit(0) && x.testBit(1) && y.testBit(0) && y.testBit(1)) { // x = y = 3 (mod 4)
                J = -J;
            }
            var t = new BigInteger();
            t = x; x = y; y = t; // switch x and y
        }
        return J;
    };

    /*
     * Compute a square root of v (mod p).
     * @param {BigInteger} v
     * @returns a square root of v (mod p) if one exists, or null otherwise.
     */
    JSParams.sqrt = function (v) {
        if (v.signum() === 0) {
            return this._0;
        }
        // case I: p = 3 (mod 4):
        if (this.p.testBit(1)) {
            var r = v.modPow(this.p.shiftRight(2).add(this._1), this.p);
            // test solution:
            return r.multiply(r).subtract(v).mod(this.p).signum() === 0 ? r : null;
        }
        // case II: p = 5 (mod 8):
        if (this.p.testBit(2)) {
            var twog = v.shiftLeft(1).mod(this.p);
            var gamma = twog.modPow(this.p.shiftRight(3), this.p);
            var i = twog.multiply(gamma).multiply(gamma).mod(this.p);
            var r = v.multiply(gamma).multiply(i.subtract(this._1)).mod(this.p);
            // test solution:
            return r.multiply(r).subtract(v).mod(this.p).signum() === 0 ? r : null;
        }
        // case III: p = 9 (mod 16):
        if (this.p.testBit(3)) {
            var twou = this.p.shiftRight(2); // (p-1)/4
            var s0 = v.shiftLeft(1).modPow(twou, this.p); // (2v)^(2u) mod p
            var s = s0;
            var d = this._1;
            var fouru = twou.shiftLeft(1);
            while (s.add(this._1).compareTo(this.p) !== 0) {
                d = d.add(this._2);
                s = d.modPow(fouru, this.p).multiply(s0).mod(this.p);
            }
            var w = d.multiply(d).multiply(v).shiftLeft(1).mod(this.p);
            var z = w.modPow(this.p.shiftRight(4), this.p); // w^((p-9)/16)
            var i = z.multiply(z).multiply(w).mod(this.p);
            var r = z.multiply(d).multiply(v).multiply(i.subtract(this._1)).mod(this.p);
            // test solution:
            return r.multiply(r).subtract(v).mod(this.p).signum() === 0 ? r : null;
        }
        // case IV: p = 17 (mod 32):
        if (this.p.testBit(4)) {
            var twou = this.p.shiftRight(3); // (p-1)/8
            var s0 = v.shiftLeft(1).modPow(twou, this.p); // (2v)^(2u) mod p
            var s = s0;
            var d = this._1;
            var fouru = twou.shiftLeft(1); // (p-1)/4
            while (s.add(this._1).compareTo(this.p) !== 0) {
                d = d.add(this._2);
                s = d.modPow(fouru, this.p).multiply(s0).mod(this.p);
            }
            var w = d.multiply(d).multiply(v).shiftLeft(1).mod(this.p);
            var z = w.modPow(this.p.shiftRight(5), this.p); // w^((p-17)/32)
            var i = z.multiply(z).multiply(w).mod(this.p);
            var r = z.multiply(d).multiply(v).multiply(i.subtract(this._1)).mod(this.p);
            // test solution:
            return r.multiply(r).subtract(v).mod(this.p).signum() === 0 ? r : null;
        }
        // case V: p = 1 (mod 4, 8, 16, 32):
        if (v.compareTo(this._4) === 0) {
            return this._2;
        }
        var z = v.subtract(this._4).mod(this.p);
        var t = this._1;
        while (this.legendre(z) >= 0) {
            t = t.add(this._1);
            z = v.multiply(t).multiply(t).subtract(this._4).mod(this.p);
        }
        var z = v.multiply(t).multiply(t).subtract(this._2).mod(this.p);
        var r = this.lucas(z, this.p.shiftRight(2)).multiply(t.modInverse(this.p)).mod(this.p);
        // test solution:
        return r.multiply(r).subtract(v).mod(this.p).signum() === 0 ? r : null;
    };

    /*
     * Compute a cube root of v (mod p) where p = 4 (mod 9).
     * @param {BigInteger} v
     * @returns a cube root of v (mod p) if one exists, or null otherwise.
     */
    JSParams.cbrt = function (v) {
        if (this.p.mod(this._9).intValue() !== 4) {
            throw new Error("This implementation is optimized for, and only works with, prime fields GF(p) where p = 4 (mod 9)");
        }
        if (v.signum() === 0) {
            return this._0;
        }
        var r = v.modPow(this.cbrtExponent, this.p); // r = v^{(2p + 1)/9}
        return r.multiply(r).multiply(r).subtract(v).mod(this.p).signum() === 0 ? r : null;
    };

    /*
     * Postl's algorithm to compute V_k(P, 1)
     * @param {BigInteger} P
     * @param {BigInteger} k
     * @returns result almorithm
     */
    JSParams.lucas = function (P, k) {
        var d_1 = P;
        var d_2 = P.multiply(P).subtract(this._2).mod(this.p);
        var l = k.bitLength() - 1; // k = \sum_{j=0}^l{k_j*2^j}
        for (var j = l - 1; j >= 1; j--) {
            if (k.testBit(j)) {
                d_1 = d_1.multiply(d_2).subtract(P).mod(this.p);
                d_2 = d_2.multiply(d_2).subtract(this._2).mod(this.p);
            } else {
                d_2 = d_1.multiply(d_2).subtract(P).mod(this.p);
                d_1 = d_1.multiply(d_1).subtract(this._2).mod(this.p);
            }
        }
        return (k.testBit(0)) ? d_1.multiply(d_2).subtract(P).mod(this.p) : d_1.multiply(d_1).subtract(this._2).mod(this.p);
    };

    /*
     * Random BigBigInteger
     * @param {Integer} k
     * @param {SecureRandom} rnd
     * @returns new BigInteger
     */
    JSParams.randomBigInteger = function (k, rnd) {
        return new BigInteger(k, rnd);
    };

    /*
     * Compute the quadratic character of a, i.e. the legendre symbol of a
     * @param {BigInteger} a
     * @param {JSParams} bn
     * @return 0, 1 or -1
     */
    JSParams.chi_q = function (a, bn) {
        var arg = a.mod(bn.p);
        if (arg.equals(BigInteger.ZERO)) {
            console.log("argument is zero in F_p!");
            return this._0;
        }
        return bn.legendre(arg) === 1 ? this._1 : this._1.negate();
    };

    /*
     * Obtain Module
     */
    JSParams.getModulus = function () {
        return this.p;
    };

    /*
     * Obtain Curve Order
     */
    JSParams.getCurveOrder = function () {
        return this.n;
    };

    JSParams.prototype.legendre = JSParams.legendre;
    JSParams.prototype.sqrt = JSParams.sqrt;
    JSParams.prototype.cbrt = JSParams.cbrt;
    JSParams.prototype.lucas = JSParams.lucas;
    JSParams.prototype.randomBigInteger = JSParams.randomBigInteger;
    JSParams.prototype.chi_q = JSParams.chi_q;
//    JSParams.prototype.SWEncJS = JSParams.SWEncJS;
    JSParams.prototype.getModulus = JSParams.getModulus;
    JSParams.prototype.getCurveOrder = JSParams.getCurveOrder;

    return JSParams;
})();
