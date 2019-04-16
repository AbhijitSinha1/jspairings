const BigInteger = require('./BigNumber').BigInteger;
const SecureRandom = require('./RNG');
const JSField2 = require('./JSField2');

/*
 * JSField4.js
 * Arithmetic in the finite extension field GF(p^4) with p = 3 (mod 4).
 * Copyright (C) Jorge H. C. Lucema.
 */

var differentFields = "Operands are in different finite fields";

var JSField4 = function (bn, re, im) {
    /*
        * 
        */
    if (arguments.length === 1) {
        this.bn = bn;
        this.re = this.im = bn.Fp2_0;
    }
    if (arguments.length === 2) {
        /*
            * 
            */
        if (re instanceof JSField2) {
            this.bn = bn;
            this.re = re;
            this.im = bn.Fp2_0;
        }
        /*
            * 
            */
        else if (re instanceof BigInteger()) {
            var k = re;
            this.bn = bn;
            this.re = new JSField2(bn, k);
            this.im = bn.Fp2_0;
        }
        /*
            * Compute a random field element.
            * @param {JSParams} bn BN parameters
            * @param {SecureRandom} rand a cryptographically strong pseudo-random number generator.
            * @returns a random field element.
            */
        else if (re instanceof SecureRandom()) {
            var rand = re;
            this.bn = bn;
            this.re = new JSField2(bn, rand);
            this.im = new JSField2(bn, rand);
        }
    }
    /*
        * 
        */
    if (arguments.length === 3) {
        this.bn = bn;
        this.re = re;
        this.im = im;
    }
};

/*
    * 
    * @param {SecureRandom} rand a cryptographically strong pseudo-random number generator.
    * @returns a random field element.
    */
JSField4.randomize = function (rand) {
    return new JSField4(this.bn, rand);
};

/*
    * 
    * @returns comparison with 0
    */
JSField4.isZero = function () {
    return this.re.isZero() && this.im.isZero();
};

/*
    * 
    * @returns comparison with 1
    */
JSField4.isOne = function () {
    return this.re.isOne() && this.im.isZero();
};

/*
    * 
    * @param {Object} o
    * @returns comparison objects
    */
JSField4.equals = function (o) {
    if (!(o instanceof JSField4)) {
        return false;
    }
    var w = o;
    return this.bn === w.bn && // singleton comparison
        this.re.equals(w.re) && this.im.equals(w.im);
};

/*
    * 
    * @returns a Field4 negate
    */
JSField4.negate = function () {
    return new JSField4(this.bn, this.re.negate(), this.im.negate());
};

/*
    * 
    * @param {JSField4} w
    * @returns 
    */
JSField4.add = function (w) {
    if (this.bn !== w.bn) { // singleton comparison
        throw new Error(differentFields);
    }
    return new JSField4(this.bn, this.re.add(w.re), this.im.add(w.im));
};

/*
    * 
    * @param {JSField4} w
    * @returns 
    */
JSField4.subtract = function (w) {
    if (this.bn !== w.bn) { // singleton comparison
        throw new Error(differentFields);
    }
    return new JSField4(this.bn, this.re.subtract(w.re), this.im.subtract(w.im));
};

/*
    * 
    * @param {Integer} k
    * @returns 
    */
JSField4.twice = function (k) {
    return new JSField4(this.bn, this.re.twice(k), this.im.twice(k));
};

/*
    * 
    * @returns 
    */
JSField4.halve = function () {
    return new JSField4(this.bn, this.re.halve(), this.im.halve());
};

/*
    * (re + im*v)*v = im*xi + re*v
    * @returns 
    */
JSField4.multiplyV = function () {
    return new JSField4(this.bn, this.im.multiplyV(), this.re);
};

/*
    * 
    * @returns 
    */
JSField4.divideV = function () {
    return new JSField4(this.bn, this.im.divideV(), this.re);
};

JSField4.multiply = function (w) {
    /*
        * 
        * @param {JSField4} w
        * @returns
        */
    if (w instanceof JSField4) {
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
        var b0 = this.re;
        var b1 = this.im;
        var c0 = w.re;
        var c1 = w.im;
        // [b0, b1]*[c0, c1] = [b0*c0 + b1*c1*xi, (b0 + b1)*(c0 + c1) - b0*c0 - b1*c1]
        var b0c0 = b0.multiply(c0);
        var b1c1 = b1.multiply(c1);
        return new JSField4(this.bn,
            (this.bn.b === 3) ? b0c0.add(b1c1.divideV()) :  b0c0.add(b1c1.multiplyV()),
            b0.add(b1).multiply(c0.add(c1)).subtract(b0c0).subtract(b1c1));
    }
    /*
        * 
        */
    if (w instanceof JSField2) {
        if (this.bn !== w.bn) { // singleton comparison
            throw new Error(differentFields);
        }
        if (w.isOne()) {
            return this;
        }
        return new JSField4(this.bn, this.re.multiply(w), this.im.multiply(w));
    }
};

/*
    * 
    * @returns 
    */
JSField4.square = function () {
    if (this.isZero() || this.isOne()) {
        return this;
    }
    var a0 = this.re;
    var a1 = this.im;
    // [a0, a1]^2 = [a0^2 + a1^2*xi, (a0 + a1)^2 - a0^2 - a1^2]
    var a02 = a0.square();
    var a12 = a1.square();
    return new JSField4 (this.bn,
        (this.bn.b === 3) ? a02.add(a12.divideV()) : a02.add(a12.multiplyV()),
        a0.add(a1).square().subtract(a02).subtract(a12));
};

/*
    * (x + ys)^{-1} = (x - ys)/(x^2 - y^2*xi)
    * @returns
    */
JSField4.inverse = function () {
    var d = this.re.square().subtract(this.im.square().multiplyV());
    return new JSField4(this.bn, this.re.multiply(d), this.im.multiply(d).negate());
};

JSField4.prototype.randomize = JSField4.randomize;
JSField4.prototype.isZero = JSField4.isZero;
JSField4.prototype.isOne = JSField4.isOne;
JSField4.prototype.equals = JSField4.equals;
JSField4.prototype.negate = JSField4.negate;
JSField4.prototype.add = JSField4.add;
JSField4.prototype.subtract = JSField4.subtract;
JSField4.prototype.twice = JSField4.twice;
JSField4.prototype.halve = JSField4.halve;
JSField4.prototype.multiplyV = JSField4.multiplyV;
JSField4.prototype.divideV = JSField4.divideV;
JSField4.prototype.multiply = JSField4.multiply;
JSField4.prototype.square = JSField4.square;
JSField4.prototype.inverse = JSField4.inverse;

module.exports = JSField4;