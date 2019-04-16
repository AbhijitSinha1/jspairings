const SecureRandom = require('../libs/RNG');
const JSParams = require('../libs/JSParams');
const JSCurve = require('../libs/JSCurve');
const JSCurve2 = require('../libs/JSCurve2');
const BigInteger = require('../libs/BigNumber').BigInteger;
const JSPairing = require('../libs/JSPairing');

/*
 * JSTests.js
 * Simple tests for Barreto-Naehrig (BN) pairing-friendly elliptic curves.
 * Copyright (C) Jorge H. C. Lucema.
 */

var JSTests = function () {
    this.validBitsRange = new Array(160,192,224,256);
};

JSTests.multiplyFp = function () {
    var rng = new SecureRandom();
    var total = new Array(this.validBitsRange.length);
    for (var i = 0; i < this.validBitsRange.length; i++) {
        total[i] = 0;
        var bn = new JSParams(this.validBitsRange[i]);
        var E = new JSCurve(bn);
        var p = E.pointFactory(rng);
        for (var j = 0; j < 1000; j++) {
            var m = new BigInteger(this.validBitsRange[i], rng);
            var before = Date.now();
            var result = p.multiply(m);
            var after = Date.now();
            total[i] += (after - before);
        }
    }
    return total;
};

JSTests.multiplyFp2 = function () {
    var rng = new SecureRandom();
    var total = new Array(this.validBitsRange.length);
    for (var i = 0; i < this.validBitsRange.length; i++) {
        total[i] = 0;
        var bn = new JSParams(this.validBitsRange[i]);
        var E = new JSCurve(bn);
        var Et = new JSCurve2(E);
        var p = Et.pointFactory(rng);
        for (var j = 0; j < 1000; j++) {
            var m = new BigInteger(this.validBitsRange[i], rng);
            var before = Date.now();
            var result = p.multiply(m);
            var after = Date.now();
            total[i] += (after - before);
        }
    }
    return total;
};

JSTests.precomputedMultiplyFp = function () {
    var rng = new SecureRandom();
    var total = new Array(this.validBitsRange.length);
    for (var i = 0; i < this.validBitsRange.length; i++) {
        total[i] = 0;
        var bn = new JSParams(this.validBitsRange[i]);
        var E = new JSCurve(bn);
        var p = E.pointFactory(rng);
        p.getSerializedTable();
        for (var j = 0; j < 1000; j++) {
            var m = new BigInteger(this.validBitsRange[i], rng);
            var before = Date.now();
            var result = p.multiply(m);
            var after = Date.now();
            total[i] += (after - before);
        }
    }
    return total;
};

JSTests.addFp = function () {
    var rng = new SecureRandom();
    var total = new Array(this.validBitsRange.length);
    for (var i = 0; i < this.validBitsRange.length; i++) {
        total[i] = 0;
        var bn = new JSParams(this.validBitsRange[i]);
        var E = new JSCurve(bn);
        var p = E.pointFactory(rng);
        var q = E.pointFactory(rng);
        for (var j = 0; j < 1000; j++) {
            var before = Date.now();
            var result = p.add(q);
            var after = Date.now();
            total[i] += (after - before);
        }
    }
    return total;
};

JSTests.pairing = function () {
    var rng = new SecureRandom();
    var total = new Array(this.validBitsRange.length);
    for (var i = 0; i < this.validBitsRange.length; i++) {
        total[i] = 0;
        var bn = new JSParams(this.validBitsRange[i]);
        var E = new JSCurve(bn);
        var Et = new JSCurve2(E);
        var pair = new JSPairing(Et);
        for (var j = 0; j < 100; j++) {
            var m = new BigInteger(this.validBitsRange[i], rng);
            var a = Et.Gt.multiply(m);
            var before = Date.now();
            var b = pair.opt(a, E.G);
            var after = Date.now();
            total[i] += (after - before);
        }
    }
    return total;
};

JSTests.test1 = function (pairType) {
    var rng = new SecureRandom();
    for (var i = 0; i < this.validBitsRange.length; i++) {
        var bn = new JSParams(this.validBitsRange[i]);
        var E = new JSCurve(bn);
        var Et = new JSCurve2(E);
        var pair = new JSPairing(Et);
        var Q = Et.Gt;
        var nQ = Q.multiply(bn.n).normalize();
        var tQ = Q.multiply(bn.t.subtract(bn._1)).normalize();
        var pQ = Q.multiply(bn.p).normalize();
        var fbx = Q.frobex(1);
        var p2Q = Q.multiply(bn.p.pow(2)).normalize();
        var fbx2 = Q.frobex(2);
        var p3Q = Q.multiply(bn.p.pow(3)).normalize();
        var fbx3 = Q.frobex(3);
        console.log("n*Q = [(" + nQ.x.re + ", " + nQ.x.im + ") : (" + nQ.y.re + ", " + nQ.y.im + ") : (" + nQ.z.re + ", " + nQ.z.im + ")]");
        console.log("(t-1)*Q = [(" + tQ.x.re + ", " + tQ.x.im + ") : (" + tQ.y.re + ", " + tQ.y.im + ") : (" + tQ.z.re + ", " + tQ.z.im + ")]");
        console.log("p*Q = [(" + pQ.x.re + ", " + pQ.x.im + ") : (" + pQ.y.re + ", " + pQ.y.im + ") : (" + pQ.z.re + ", " + pQ.z.im + ")]");
        console.log("frobex(Q) = [(" + fbx.x.re + ", " + fbx.x.im + ") : (" + fbx.y.re + ", " + fbx.y.im + ") : (" + fbx.z.re + ", " + fbx.z.im + ")]");
        console.log("frobex(Q) in E' = " + Et.contains(Q.frobex(1)));
        console.log("p^2*Q = [(" + p2Q.x.re + ", " + p2Q.x.im + ") : (" + p2Q.y.re + ", " + p2Q.y.im + ") : (" + p2Q.z.re + ", " + p2Q.z.im + ")]");
        console.log("frobex(Q) = [(" + fbx2.x.re + ", " + fbx2.x.im + ") : (" + fbx2.y.re + ", " + fbx2.y.im + ") : (" + fbx2.z.re + ", " + fbx2.z.im + ")]");
        console.log("frobex(Q) in E' = " + Et.contains(Q.frobex(2)));
        console.log("p^3*Q = [(" + p3Q.x.re + ", " + p3Q.x.im + ") : (" + p3Q.y.re + ", " + p3Q.y.im + ") : (" + p3Q.z.re + ", " + p3Q.z.im + ")]");
        console.log("frobex(Q) = [(" + fbx3.x.re + ", " + fbx3.x.im + ") : (" + fbx3.y.re + ", " + fbx3.y.im + ") : (" + fbx3.z.re + ", " + fbx3.z.im + ")]");
        console.log("frobex(Q) in E' = " + Et.contains(Q.frobex(3)));
        Q = Q.multiply(bn._6);
        if (!Q.multiply(bn.p).equals(Q.normalize().frobex(1))) {
            throw new Error("inconsistent Frobenius");
        }
        switch (pairType) {
            case "tate":
                var g = pair.tate(E.G, Et.Gt);
                break;
            case "ate":
                var g = pair.ate(Et.Gt, E.G);
                break;
            case "eta":
                var g = pair.eta(E.G, Et.Gt);
                break;
            case "opt":
        var g = pair.opt(Et.Gt, E.G);
                break;
            default:
                throw new Error("Invalid pairing type.");
        }
        console.log("(JSField12)((" + g.v[0].re + "," + g.v[0].im + "),(" + g.v[1].re + "," + g.v[1].im + "),(" + g.v[2].re + "," + g.v[2].im + "),(" +
                        g.v[3].re + "," + g.v[3].im + "),(" + g.v[4].re + "," + g.v[4].im + "),(" + g.v[5].re + "," + g.v[5].im + "))");
        if (g.isZero()) {
            throw new Error("degeneracy error!");
        }
        if (!g.exp(bn.n).isOne()) {
            throw new Error("G_T order error!");
        }
        switch (pairType) {
            case "tate":
                var a = pair.tate(E.G.twice(1).normalize(), Et.Gt);
                var b = pair.tate(E.G, Et.Gt.twice(1).normalize());
                break;
            case "ate":
                var a = pair.ate(Et.Gt.twice(1), E.G);
                var b = pair.ate(Et.Gt, E.G.twice(1));
                break;
            case "eta":
                var a = pair.eta(E.G.twice(1).normalize(), Et.Gt);
                var b = pair.eta(E.G, Et.Gt.twice(1).normalize());
                break;
            case "opt":
                var a = pair.opt(Et.Gt.twice(1), E.G);
                var b = pair.opt(Et.Gt, E.G.twice(1));
                break;
            default:
                throw new Error("Invalid pairing type.");
        }
        var c = g.square();
        console.log("bilinear? " + (a.equals(b) && b.equals(c)));
        if (!(a.equals(b) && b.equals(c)) || a.isOne()) {
            console.log(">>>> a = " + a);
            console.log(">>>> b = " + b);
            console.log(">>>> c = " + c);
            throw new Error("Bilinearity error!");
        }
        for (var j = 0; j < 10; j++) {
            var m = new BigInteger(this.validBitsRange[i], rng);
            switch (pairType) {
                case "tate":
                    a = pair.tate(E.G.multiply(m), Et.Gt);
                    b = pair.tate(E.G, Et.Gt.multiply(m));
                    break;
                case "ate":
                    a = pair.ate(Et.Gt.multiply(m), E.G);
                    b = pair.ate(Et.Gt, E.G.multiply(m));
                    break;
                case "eta":
                    a = pair.eta(E.G.multiply(m), Et.Gt);
                    b = pair.eta(E.G, Et.Gt.multiply(m));
                    break;
                case "opt":
                    a = pair.opt(Et.Gt.multiply(m), E.G);
                    b = pair.opt(Et.Gt, E.G.multiply(m));
                    break;
                default:
                    throw new Error("Invalid pairing type.");
            }
            c = g.exp(m);
            console.log("bilinear? " + (a.equals(b) && b.equals(c)));
            if (!(a.equals(b) && b.equals(c)) || a.isOne()) {
                console.log(">>>> a = " + a);
                console.log(">>>> b = " + b);
                console.log(">>>> c = " + c);
                throw new Error("Bilinearity error!");
            }
        }
    }
};

JSTests.test2 = function (iterations) {
    var numBits = 256; // caveat: maybe using larger values is better
    console.log("Testing E(F_p) arithmetic...");
    for (var i = 0; i < this.validBitsRange.length; i++) {
        var bn = new JSParams(this.validBitsRange[i]);
        var E = new JSCurve(bn);
        var rand = new SecureRandom();
        for (var j = 0; j < iterations; j++) {
            console.log("test #" + j);
            var x = E.G.randomize(rand);
            var y = E.G.randomize(rand);
            var z = E.G.randomize(rand);
            var ecZero = E.G.E.infinity;
            var m = new BigInteger(numBits, rand);
            var n = new BigInteger(numBits, rand);
            if (iterations === 1) {
                console.log("\nchecking cloning/comparison/pertinence");
            }
            if (!x.equals(x)) {
                throw new Error("Comparison failure");
            }
            if (!x.isOnSameCurve(x)) {
                throw new Error("Inconsistent pertinence self-comparison");
            }
            if (!x.E.contains(x)) {
                throw new Error("Inconsistent curve pertinence");
            }
            // check addition properties:
            if (iterations === 1) {
                console.log(" done.\nchecking addition properties");
            }
            if (!x.add(y).equals(y.add(x))) {
                throw new Error("x + y != y + x");
            }
            if (!x.add(ecZero).equals(x)) {
                throw new Error("x + 0 != x");
            }
            if (!x.add(x.negate()).isZero()) {
                throw new Error("x + (-x) != 0");
            }
            if (!x.add(y).add(z).equals(x.add(y.add(z)))) {
                throw new Error("(x + y) + z != x + (y + z)");
            }
            if (!x.negate().negate().equals(x)) {
                throw new Error("-(-x) != x");
            }
            // check scalar multiplication properties:
            if (iterations === 1) {
                console.log(" done.\nchecking scalar multiplication properties");
            }
            if (!x.multiply(new BigInteger("0")).equals(ecZero)) {
                throw new Error("0*x != 0");
            }
            if (!x.multiply(new BigInteger("1")).equals(x)) {
                throw new Error("1*x != x");
            }
            if (!x.multiply(new BigInteger("2")).equals(x.twice(1))) {
                throw new Error("2*x != twice x");
            }
            if (!x.multiply(new BigInteger("2")).equals(x.add(x))) {
                throw new Error("2*x != x + x");
            }
            if (!x.multiply(new BigInteger("-1")).equals(x.negate())) {
                throw new Error("(-1)*x != -x");
            }
            if (!x.multiply(m.negate()).equals(x.negate().multiply(m))) {
                throw new Error("(-m)*x != m*(-x)");
            }
            if (!x.multiply(m.negate()).equals(x.multiply(m).negate())) {
                throw new Error("(-m)*x != -(m*x)");
            }
            if (!x.multiply(m.add(n)).equals(x.multiply(m).add(x.multiply(n)))) {
                throw new Error("(m + n)*x != m*x + n*x");
            }
            var w = x.multiply(n).multiply(m);
            if (!w.equals(x.multiply(m).multiply(n))) {
                throw new Error("m*(n*x) != n*(m*x)");
            }
            if (!w.equals(x.multiply(m.multiply(n)))) {
                throw new Error("m*(n*x) != (m*n)*x");
            }
        }
    }        
};

JSTests.test3 = function (iterations) {
    var numBits = 256; // caveat: maybe using larger values is better
    console.log("Testing E'(F_{p^2}) arithmetic...");
    for (var i = 0; i < this.validBitsRange.length; i++) {
        var bn = new JSParams(this.validBitsRange[i]);
        var E = new JSCurve(bn);
        var Et = new JSCurve2(E);
        var rand = new SecureRandom();
        for (var j = 0; j < iterations; j++) {
            console.log("test #" + j);
            var x = Et.Gt.randomize(rand);
            var y = Et.Gt.randomize(rand);
            var z = Et.Gt.randomize(rand);
            var ecZero = Et.Gt.E.infinity;
            var m = new BigInteger(numBits, rand);
            var n = new BigInteger(numBits, rand);
            if (iterations === 1) {
                console.log("\nchecking cloning/comparison/pertinence");
            }
            if (!x.equals(x)) {
                throw new Error("Comparison failure");
            }
            if (!x.isOnSameCurve(x)) {
                throw new Error("Inconsistent pertinence self-comparison");
            }
            if (!x.E.contains(x)) {
                throw new Error("Inconsistent curve pertinence");
            }
            // check addition properties:
            if (iterations === 1) {
                console.log(" done.\nchecking addition properties");
            }
            if (!x.twice(1).equals(x.add(x))) {
                throw new Error("2*x != x + x");
            }
            if (!x.add(y).equals(y.add(x))) {
                throw new Error("x + y != y + x");
            }
            if (!x.add(ecZero).equals(x)) {
                throw new Error("x + 0 != x");
            }
            if (!x.add(x.negate()).isZero()) {
                throw new Error("x + (-x) != 0");
            }
            if (!x.add(y).add(z).equals(x.add(y.add(z)))) {
                throw new Error("(x + y) + z != x + (y + z)");
            }
            if (!x.negate().negate().equals(x)) {
                throw new Error("-(-x) != x");
            }
            // check scalar multiplication properties:
            if (iterations === 1) {
                console.log(" done.\nchecking scalar multiplication properties");
            }
            if (!x.multiply(new BigInteger("0")).equals(ecZero)) {
                throw new Error("0*x != 0");
            }
            if (!x.multiply(new BigInteger("1")).equals(x)) {
                throw new Error("1*x != x");
            }
            if (!x.multiply(new BigInteger("2")).equals(x.twice(1))) {
                throw new Error("2*x != twice x");
            }
            if (!x.multiply(new BigInteger("2")).equals(x.add(x))) {
                throw new Error("2*x != x + x");
            }
            if (!x.multiply(new BigInteger("-1")).equals(x.negate())) {
                throw new Error("(-1)*x != -x");
            }
            if (!x.multiply(m.negate()).equals(x.negate().multiply(m))) {
                throw new Error("(-m)*x != m*(-x)");
            }
            if (!x.multiply(m.negate()).equals(x.multiply(m).negate())) {
                throw new Error("(-m)*x != -(m*x)");
            }
            if (!x.multiply(m.add(n)).equals(x.multiply(m).add(x.multiply(n)))) {
                throw new Error("(m + n)*x != m*x + n*x");
            }
            var w = x.multiply(n).multiply(m);
            if (!w.equals(x.multiply(m).multiply(n))) {
                throw new Error("m*(n*x) != n*(m*x)");
            }
            if (!w.equals(x.multiply(m.multiply(n)))) {
                throw new Error("m*(n*x) != (m*n)*x");
            }
            if (!x.multiply(bn.p).equals(x.normalize().frobex(1))) {
                console.log("x^p    = " + x.multiply(bn.p));
                console.log("Phi(x) = " + x.normalize().frobex(1));
                throw new Error("inconsistent Frobenius");
            }
            if (!x.multiply(bn.p).multiply(bn.p).equals(x.normalize().frobex(2))) {
                throw new Error("inconsistent Frobenius");
            }
            if (!x.multiply(bn.p).multiply(bn.p).multiply(bn.p).equals(x.normalize().frobex(3))) {
                throw new Error("inconsistent Frobenius");
            }
        }
    }
};

JSTests.test4 = function (iterations) {
    var numBits = 256; // caveat: maybe using larger values is better
    console.log("Testing F_{p^12} arithmetic...");
    var before = Date.now();
    for (var i = 0; i < this.validBitsRange.length; i++) {
        var bn = new JSParams(this.validBitsRange[i]);
        var E = new JSCurve(bn);
        var rand = new SecureRandom();
        for (var j = 0; j < iterations; j++) {
            var f = E.bn.Fp12_0.randomize(rand);
            var g = E.bn.Fp12_0.randomize(rand);
            var h = E.bn.Fp12_0.randomize(rand);
            var m = new BigInteger(numBits, rand);
            var n = new BigInteger(numBits, rand);
            // addition/subtraction tests
            if (!f.add(E.bn.Fp12_0.bn.Fp12_0).equals(f)) {
                throw new Error("Inconsistent Fp12 field addition");
            }
            if (!f.add(f.negate()).isZero()) {
                throw new Error("Inconsistent Fp12 field addition");
            }
            if (!f.subtract(g).equals(f.add(g.negate()))) {
                throw new Error("Inconsistent Fp12 field addition");
            }
            if (!f.subtract(g).negate().equals(g.subtract(f))) {
                throw new Error("Inconsistent Fp12 field addition");
            }
            if (!f.add(g).add(h).equals(f.add(g.add(h)))) {
                throw new Error("Inconsistent Fp12 field addition");
            }
            if (!f.add(g).equals(g.add(f))) {
                throw new Error("Inconsistent Fp12 field addition");
            }
            if (!f.add(g).subtract(g).equals(f)) {
                throw new Error("Inconsistent Fp12 field addition");
            }
            // multiplication tests
            if (!f.multiply(g).multiply(h).equals(f.multiply(g.multiply(h)))) {
                throw new Error("Inconsistent Fp12 field multiplication");
            }
            if (!f.multiply(E.bn.Fp12_0.bn.Fp12_0).isZero()) {
                throw new Error("Inconsistent Fp12 field multiplication");
            }
            if (!f.multiply(E.bn.Fp12_0.bn.Fp12_1).equals(f)) {
                throw new Error("Inconsistent Fp12 field multiplication");
            }
            if (!f.multiply(g).equals(g.multiply(f))) {
                throw new Error("Inconsistent Fp12 field multiplication");
            }
            if (!f.multiply(f).equals(f.square())) {
                throw new Error("Inconsistent Fp12 field multiplication");
            }
            // inversion tests
            var z = f.inverse();
            if (!f.multiply(z).isOne()) {
                throw new Error("Inconsistent Fp12 field inversion");
            }
            if (!f.multiply(g.multiply(z)).equals(g)) {
                throw new Error("Inconsistent Fp12 field inversion");
            }
            // distribution tests
            if (!f.multiply(g.add(h)).equals(f.multiply(g).add(f.multiply(h)))) {
                throw new Error("Inconsistent Fp12 field distribution");
            }
            // exponentiation tests
            if (!f.exp(m).exp(n).equals(f.exp(n).exp(m))) {
                throw new Error("Inconsistent Fp12 field exponentiation");
            }
            if (!f.exp(m).exp(n).equals(f.exp(m.multiply(n)))) {
                throw new Error("Inconsistent Fp12 field exponentiation");
            }
            if (!f.exp(m).multiply(f.exp(n)).equals(f.exp(m.add(n)))) {
                throw new Error("Inconsistent Fp12 field exponentiation");
            }
            if (!f.frobenius().equals(f.exp(E.bn.Fp12_0.bn.p))) {
                console.log("frob(f) = " + f.frobenius());
                console.log("f^p     = " + f.exp(E.bn.Fp12_0.bn.p));
                throw new Error("Inconsistent Fp12 field Frobenius");
            }
            z = f;
            for (var k = 0; k < 6; k++) {
                if (!f.conjugate(k).equals(z)) {
                    console.log("f.conjugate(" + k + ") = " + f.conjugate(k));
                    console.log("f.exp((p^2)^" + k + ") = " + z);
                    throw new Error("Inconsistent Fp12 field conjugate");
                }
                z = z.exp(E.bn.Fp12_0.bn.p).exp(E.bn.Fp12_0.bn.p);
            }
        }
        var after = Date.now();
        console.log("OK; all " + iterations + " tests done in " + (after-before)/1000 + " s.");
    }
};

JSTests.prototype.multiplyFp = JSTests.multiplyFp;
JSTests.prototype.multiplyFp2 = JSTests.multiplyFp2;
JSTests.prototype.precomputedMultiplyFp = JSTests.precomputedMultiplyFp;
JSTests.prototype.addFp = JSTests.addFp;
JSTests.prototype.pairing = JSTests.pairing;
JSTests.prototype.test1 = JSTests.test1;
JSTests.prototype.test2 = JSTests.test2;
JSTests.prototype.test3 = JSTests.test3;
JSTests.prototype.test4 = JSTests.test4;

module.exports = JSTests;