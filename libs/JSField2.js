const BigInteger = require('./BigNumber').BigInteger;
const SecureRandom = require('./RNG');

/*
 * JSField2.js
 * Arithmetic in the finite extension field GF(p^2) with p = 3 (mod 4) and p = 4 (mod 9).
 * Copyright (C) Jorge H. C. Lucema.
 */

var differentFields = "Operands are in different finite fields";
var SMALLSIZE = 32;

var JSField2 = function (bn, re, im, reduce) {

    // statistics
    this.addcount = 0;
    this.mulcount = 0;
    this.sqrcount = 0;
    this.modcount = 0;
    this.fpmcount = 0;
    this.docount = false;
    this.modenable = true;

    /*
        * 
        */
    if(arguments.length === 1) {
        this.bn = bn;
        this.re = bn._0;
        this.im = bn._0;
    }
    /*
        * 
        */
    if(arguments.length === 2) {
        /*
            * 
            */
        if (re instanceof BigInteger) {
            this.bn = bn;
            this.re = re; // caveat: no modular reduction!
            this.im = bn._0;
        }
        /*
            * 
            */
        else if (re instanceof SecureRandom) {
            this.bn = bn;
            var rand = re;
            do {
                this.re = new BigInteger(this.bn.p.bitLength(), rand);
            } while (this.re.compareTo(this.bn.p) >= 0);
            do {
                this.im = new BigInteger(this.bn.p.bitLength(), rand);
            } while (this.im.compareTo(this.bn.p) >= 0);
        }
    }
    /*
        * 
        */
    if(arguments.length === 4) {
        this.bn = bn;
        if (reduce) {
            if (this.docount && this.modenable) {
                if (re.signum() < 0 || re.compareTo(bn.p) >= 0) {
                    this.modcount++;
                }
                if (im.signum() < 0 || im.compareTo(bn.p) >= 0) {
                    this.modcount++;
                }
            }
            this.re = re.mod(bn.p);
            this.im = im.mod(bn.p);
        } else {
            this.re = re;
            this.im = im;
        }
    }
};

/*
    * 
    * @param {boolean} enable
    * @returns true
    */
JSField2.counton = function (enable) {
    if (enable) {
        this.docount = true;
    }
};

/*
    * 
    * @param {boolean} enable
    * @returns false
    */
JSField2.countoff = function (enable) {
    if (enable) {
        this.docount = false;
    }
};

/*
    * 
    * @returns true
    */
JSField2.modon = function () {
    this.modenable = true;
};

/*
    * 
    * @returns false
    */
JSField2.modoff = function () {
    this.modenable = false;
};

/*
    * 
    * @returns 0
    */
JSField2.reset = function () {
    this.addcount = 0;
    this.mulcount = 0;
    this.sqrcount = 0;
    this.modcount = 0;
    this.fpmcount = 0;
};

/*
    * 
    * @returns value of addcount
    */
JSField2.getadd = function () {
    return this.addcount;
};

/*
    * 
    * @returns value of mulcount
    */
JSField2.getmul = function () {
    return this.mulcount;
};

/*
    * 
    * @returns value of sqrcount
    */
JSField2.getsqr = function () {
    return this.sqrcount;
};

/*
    * 
    * @returns value of modcount/2
    */
JSField2.getmod = function () {
    return this.modcount/2;
};

/*
    * 
    * @returns value of fpmcount
    */
JSField2.getfpm = function () {
    return this.fpmcount;
};

/*
    * 
    * @param {long} delta
    * @returns increase value of modcount
    */
JSField2.offsetmod = function (delta) {
    if (this.docount) {
        this.modcount += delta;
    }
};

/*
    * Compute a random field element.
    * @param {SecureRandom} rand a cryptographically strong pseudo-random number generator.
    * @returns a random field element.
    */
JSField2.randomize = function (rand) {
    return new JSField2(this.bn, rand);
};

/*
    * 
    * @returns comparison with 0
    */
JSField2.isZero = function () {
    return this.re.signum() === 0 && this.im.signum() === 0;
};

/*
    * 
    * @returns comparison with 1
    */
JSField2.isOne = function() {
    return this.re.compareTo(this.bn._1) === 0 && this.im.signum() === 0;
};

/*
    * 
    * @param {Object} u
    * @returns comparison objects
    */
JSField2.equals = function(u) {
    if (!(u instanceof JSField2)) {
        return false;
    }
    var v = u;
    return this.bn === v.bn && // singleton comparison
        this.re.compareTo(v.re) === 0 &&
        this.im.compareTo(v.im) === 0;
};

/*
    * -(x + yi)
    * @returns a Field2 negate
    */
JSField2.negate = function () {
    return new JSField2(this.bn, (this.re.signum() !== 0) ? this.bn.p.subtract(this.re) : this.re, (this.im.signum() !== 0) ? this.bn.p.subtract(this.im) : this.im, false);
};

/*
    * (x + yi)^p = x - yi
    * @returns a Field2 conjugate
    */
JSField2.conjugate = function () {
    return new JSField2(this.bn, this.re, (this.im.signum() !== 0) ? this.bn.p.subtract(this.im) : this.im, false);
};

JSField2.add = function (v) {
    /*
        * 
        * @param {JSField2} v
        * @returns addition of two points JSField2
        */
    if (v instanceof JSField2) {
        if (this.bn !== v.bn) { // singleton comparison
            throw new Error(differentFields);
        }
        var r = this.re.add(v.re);
        if (this.docount && this.re.bitLength() > SMALLSIZE && v.re.bitLength() > SMALLSIZE) {
            this.addcount++;
        }
        if (r.compareTo(this.bn.p) >= 0) {
            r = r.subtract(this.bn.p);
        }
        var i = this.im.add(v.im);
        if (this.docount && this.im.bitLength() > SMALLSIZE && v.im.bitLength() > SMALLSIZE) {
            this.addcount++;
        }
        if (i.compareTo(this.bn.p) >= 0) {
            i = i.subtract(this.bn.p);
        }
        return new JSField2(this.bn, r, i, false);
    }
    /*
        * 
        * @param {BigInteger} v
        * @returns addition of value v with the point JSField2 this
        */
    else if (v instanceof BigInteger) {
        var s = this.re.add(v);
        if (this.docount && this.re.bitLength() > SMALLSIZE && v.bitLength() > SMALLSIZE) {
            this.addcount++;
        }
        if (s.compareTo(this.bn.p) >= 0) {
            s = s.subtract(this.bn.p);
        }
        return new JSField2(this.bn, s, this.im, false);
    }
};

JSField2.subtract = function (v) {
    /*
        * 
        * @param {JSField2} v
        * @returns difference between this and v
        */
    if (v instanceof JSField2) {
        if (this.bn !== v.bn) { // singleton comparison
            throw new Error(differentFields);
        }
        var r = this.re.subtract(v.re);
        if (this.docount && this.re.bitLength() > SMALLSIZE && v.re.bitLength() > SMALLSIZE) {
            this.addcount++;
        }
        if (r.signum() < 0) {
            r = r.add(this.bn.p);
        }
        var i = this.im.subtract(v.im);
        if (this.docount && this.im.bitLength() > SMALLSIZE && v.im.bitLength() > SMALLSIZE) {
            this.addcount++;
        }
        if (i.signum() < 0) {
            i = i.add(this.bn.p);
        }
        return new JSField2(this.bn, r, i, false);
    }
    /*
        * 
        * @param {JSField2} v
        * @returns subtract v of value real of the this
        */
    else if (v instanceof BigInteger) {
        if (this.docount && this.re.bitLength() > SMALLSIZE && v.bitLength() > SMALLSIZE) {
            this.addcount++;
        }
        var r = this.re.subtract(v);
        if (r.signum() < 0) {
            r = r.add(this.bn.p);
        }
        return new JSField2(this.bn, r, this.im, false);
    }
};

/*
    * 
    * @param {Integer} k
    * @returns 2^k*this
    */
JSField2.twice = function (k) {
    var r = this.re;
    var i = this.im;
    while (k-- > 0) {
        r = r.shiftLeft(1);
        if (r.compareTo(this.bn.p) >= 0) {
            r = r.subtract(this.bn.p);
        }
        i = i.shiftLeft(1);
        if (i.compareTo(this.bn.p) >= 0) {
            i = i.subtract(this.bn.p);
        }
    }
    return new JSField2(this.bn, r, i, false);
};

/*
    * 
    * @returns half of this
    */
JSField2.halve = function () {
    return new JSField2(this.bn,
        (this.re.testBit(0) ? this.re.add(this.bn.p) : this.re).shiftRight(1),
        (this.im.testBit(0) ? this.im.add(this.bn.p) : this.im).shiftRight(1),
        false);
};

JSField2.multiply = function (v) {
    /*
        * (x + yi)(u + vi)
        * @param {JSField2} v
        * @returns (xu - yv) + ((x + y)(u + v) - xu - yv)i
        */
    if (v instanceof JSField2) {
        if (this === v) {
            return this.square();
        }
        if (this.bn !== v.bn) { // singleton comparison
            throw new Error(differentFields);
        }
        if (this.isOne() || v.isZero()) {
            return v;
        }
        if (this.isZero() || v.isOne()) {
            return this;
        }
        //*
        if (this.docount) {
            this.addcount += 5;
            this.mulcount++;
            this.fpmcount += 3;
        }
        //*/
        var re2 = this.re.multiply(v.re); // mod p
        var im2 = this.im.multiply(v.im); // mod p
        var mix = this.re.add(this.im).multiply(v.re.add(v.im)); // mod p;
        return new JSField2(this.bn,
            re2.subtract(im2),
            mix.subtract(re2).subtract(im2),
            true);
    }
    /*
        * (x + yi)v
        * @param {BigInteger} v
        * @returns xv + yvi
        */
    else if (v instanceof BigInteger) {
        return new JSField2(this.bn, this.re.multiply(v), this.im.multiply(v), true);
    }
    /*
        * v(x + yi)
        * @param {Integer} v
        * @returns vx + vyi
        */
    else if (v instanceof Number) {
        var newre = this.re.multiply(new BigInteger(v.toString()));
        while (newre.signum() < 0) {
            newre = newre.add(this.bn.p);
        }
        while (newre.compareTo(this.bn.p) >= 0) {
            newre = newre.subtract(this.bn.p);
        }
        var newim = this.im.multiply(new BigInteger(v.toString()));
        while (newim.signum() < 0) {
            newim = newim.add(this.bn.p);
        }
        while (newim.compareTo(this.bn.p) >= 0) {
            newim = newim.subtract(this.bn.p);
        }
        return new JSField2(this.bn, newre, newim, false);
    } else {
        throw new Error("Incorrect type argument");
    }
};

/*
    * (x + yi)^2
    * @returns (x + y)(x - y) + ((x+y)^2 - (x + y)(x - y))i
    */
JSField2.square = function () {
    if (this.isZero() || this.isOne()) {
        return this;
    }
    if (this.im.signum() === 0) {
        //*
        if (this.docount) {
            this.fpmcount++;
        }
        //*/
        return new JSField2(this.bn,
            this.re.multiply(this.re), this.bn._0, true);
    }
    if (this.re.signum() === 0) {
        //*
        if (this.docount) {
            this.fpmcount++;
        }
        //*/
        return new JSField2(this.bn,
            this.im.multiply(this.im).negate(), this.bn._0, true);
    }
    //*
    if (this.docount) {
        this.addcount += 2;
        this.sqrcount++;
        this.fpmcount += 2;
    }
    //*/
    return new JSField2(this.bn,
        this.re.add(this.im).multiply(this.re.subtract(this.im)),
        this.re.multiply(this.im).shiftLeft(1),
        true);
};

/*
    * (x + yi)^3
    * @returns x(x^2 - 3y^2) + y(3x^2 - y^2)i
    */
JSField2.cube = function () {
    var re2 = this.re.multiply(this.re); // mod p
    var im2 = this.im.multiply(this.im); // mod p
    return new JSField2(this.bn,
        this.re.multiply(re2.subtract(im2.add(im2).add(im2))),
        this.im.multiply(re2.add(re2).add(re2).subtract(im2)),
        true);
};

/*
    * (x + yi)^{-1}
    * @returns (x - yi)/(x^2 + y^2)
    */
JSField2.inverse = function () {
    var d = this.re.multiply(this.re).add(this.im.multiply(this.im)).modInverse(this.bn.p);
    return new JSField2(this.bn, this.re.multiply(d), this.bn.p.subtract(this.im).multiply(d), true);
};

/*
    * (x + yi)i
    * @returns (-y + xi)
    */
JSField2.multiplyI = function () {
    return new JSField2(this.bn, (this.im.signum() !== 0) ? this.bn.p.subtract(this.im) : this.im, this.re, false);
};

/*
    * (x + yi)/i
    * @returns y - xi
    */
JSField2.divideI = function () {
    return new JSField2(this.bn, this.im, (this.re.signum() !== 0) ? bn.p.subtract(this.re) : this.re, false);
};

/*
    * (x + yi)(1 + i)
    * @returns (x - y) + (x + y)i
    */
JSField2.multiplyV = function () {
    var r = this.re.subtract(this.im);
    if (r.signum() < 0) {
        r = r.add(this.bn.p);
    }
    var i = this.re.add(this.im);
    if (i.compareTo(this.bn.p) >= 0) {
        i = i.subtract(this.bn.p);
    }
    return new JSField2(this.bn, r, i, false);
};

/*
    * 
    * @returns 
    */
JSField2.divideV = function () {
    var qre = this.re.add(this.im);
    if (qre.compareTo(this.bn.p) >= 0) {
        qre = qre.subtract(this.bn.p);
    }
    var qim = this.im.subtract(this.re);
    if (qim.signum() < 0) {
        qim = qim.add(this.bn.p);
    }
    return new JSField2(this.bn, (qre.testBit(0) ? qre.add(this.bn.p) : qre).shiftRight(1),
        (qim.testBit(0) ? qim.add(this.bn.p) : qim).shiftRight(1), false);
};

/*
    * 
    * @param {BigInteger} k
    * @returns (e^k)*this
    */
JSField2.exp = function (k) {
    var P = this;
    if (k.signum() < 0) {
        k = k.negate();
        P = P.inverse();
    }
    var e = k.toByteArray();
    var mP = new Array(16);
    mP[0] = this.bn.Fp2_1;
    mP[1] = P;
    for (var m = 1; m <= 7; m++) {
        mP[2*m] = mP[m].square();
        mP[2*m + 1] = mP[2*m].multiply(P);
    }
    var A = mP[0];
    for (var i = 0; i < e.length; i++) {
        var u = e[i] & 0xff;
        A = A.square().square().square().square().multiply(mP[u >>> 4]).square().square().square().square().multiply(mP[u & 0xf]);
    }
    return A;
};

/*
    * Compute a square root of this.
    * @returns  a square root of this if one exists, or null otherwise.
    */
JSField2.sqrt = function () {
    if (this.isZero()) {
        return this;
    }
    var r = this.exp(this.bn.sqrtExponent2); // r = v^{(p^2 + 7)/16}
    var r2 = r.square();
    if (r2.subtract(this).isZero()) {
        return r;
    }
    if (r2.add(this).isZero()) {
        return r.multiplyI();
    }
    r2 = r2.multiplyI();
    //JSField2 sqrtI = new JSField2(bn, bn.invSqrtMinus2, bn.p.subtract(bn.invSqrtMinus2), false); // sqrt(i) = (1 - i)/sqrt(-2)
    r = r.multiply(this.bn.sqrtI);
    if (r2.subtract(this).isZero()) {
        return r;
    }
    if (r2.add(this).isZero()) {
        return r.multiplyI();
    }
    return null;
};

/*
    * Compute a cube root of this.
    * @returns  a cube root of this if one exists, or null otherwise.
    */
JSField2.cbrt = function () {
    if (this.isZero()) {
        return this;
    }
    var r = this.exp(bn.cbrtExponent2); // r = v^{(p^2 + 2)/9}
    return r.cube().subtract(this).isZero() ? r : null;
};

JSField2.prototype.counton = JSField2.counton;
JSField2.prototype.countoff = JSField2.countoff;
JSField2.prototype.modon = JSField2.modon;
JSField2.prototype.modoff = JSField2.modoff;
JSField2.prototype.reset = JSField2.reset;
JSField2.prototype.getadd = JSField2.getadd;
JSField2.prototype.getmul = JSField2.getmul;
JSField2.prototype.getsqr = JSField2.getsqr;
JSField2.prototype.getmod = JSField2.getmod;
JSField2.prototype.getfpm = JSField2.getfpm;
JSField2.prototype.offsetmod = JSField2.offsetmod;
JSField2.prototype.randomize = JSField2.randomize;
JSField2.prototype.isZero = JSField2.isZero;
JSField2.prototype.isOne = JSField2.isOne;
JSField2.prototype.equals = JSField2.equals;
JSField2.prototype.negate = JSField2.negate;
JSField2.prototype.conjugate = JSField2.conjugate;
JSField2.prototype.add = JSField2.add;
JSField2.prototype.subtract = JSField2.subtract;
JSField2.prototype.twice = JSField2.twice;
JSField2.prototype.halve = JSField2.halve;
JSField2.prototype.multiply = JSField2.multiply;
JSField2.prototype.square = JSField2.square;
JSField2.prototype.cube = JSField2.cube;
JSField2.prototype.inverse = JSField2.inverse;
JSField2.prototype.multiplyI = JSField2.multiplyI;
JSField2.prototype.divideI = JSField2.divideI;
JSField2.prototype.multiplyV = JSField2.multiplyV;
JSField2.prototype.divideV = JSField2.divideV;
JSField2.prototype.exp = JSField2.exp;
JSField2.prototype.sqrt = JSField2.sqrt;
JSField2.prototype.cbrt = JSField2.cbrt;

module.exports = JSField2;