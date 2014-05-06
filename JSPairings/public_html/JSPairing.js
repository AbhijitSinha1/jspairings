/*
 * JSPairing.js
 * Bilinear pairings over Barreto-Naehrig (BN) elliptic curves.
 * Copyright (C) Jorge H. C. Lucema.
 */

var MillerLoop = true;
var FinalExp = true;

var JSPairing = (function () {

    var JSPairing = function (Et) {

        /*
         * statistics
         */
        this.addcount = 0;
        this.mulcount = 0;
        this.sqrcount = 0;
        this.modcount = 0;
        this.fpmcount = 0;

        this.E2 = Et;
        this.E = Et.E;
        this.bn = this.E.bn;
        this.Fp12_0 = this.bn.Fp12_0;
        this.Fp12_1 = this.bn.Fp12_1;
//        this.gp16 = new Array(Math.round((this.bn.n.bitLength() + 3)/4));
//        for (var i = 0; i < Math.round((this.bn.n.bitLength() + 3)/4); i++) {
//            this.gp16[i] = new Array(16);
//        }
    };

    /*
     * 
     */
    JSPairing.gl = function (V, P, Q) {
        var p = this.bn.p;
        if (V.isZero() || P.isZero() || Q.isZero()) {
            return this.Fp12_1;
        }
        var Vz3 = V.z.multiply(V.z).multiply(V.z).mod(p);
        if (V.equals(P)) {
            // y = Y/Z^3 => 1/2y = Z^3/2Y
            // x = X/Z^2 => 3x^2 = 3X^2/Z^4 =>
            // => lambda = 3x^2/2y = 3X^2/(2Y*Z)
            var n = V.x.multiply(V.x).multiply(new BigInteger("3"));//.mod(p);
            var d = V.y.multiply(V.z).shiftLeft(1);//.mod(p);
        } else {
            // lambda = (P.y - V.y)/(P.x - V.x) // P.Z = 1
            // = (P.Y - V.Y/V.Z^3) / (P.X - V.X/V.Z^2)
            // = (P.Y*V.Z^3 - V.Y) / (P.X*V.Z^3 - V.X*V.Z)
            var n = P.y.multiply(Vz3).subtract(V.y);//.mod(p);
            var d = P.x.multiply(Vz3).subtract(V.x.multiply(V.z));//.mod(p);
        }
        // lambda = n/d
        var w = new Array(6);
        //n*(Qt[1]*z^2 - V.x) + d*(V.y - Qt[2]*z^3);
        //n*Q.x*z^2 - n*V.x + d*V.y - d*Q.y*z^3;
        //(d*V.y - n*V.x) + n*Q.x*z^2 - d*Q.y*z^3;
        //(d*V.Y/V.Z^3 - n*V.X/V.Z^2) + n*Q.x*z^2 - d*Q.y*z^3;
        //(d*V.Y - n*V.X*V.Z) + n*Q.x*V.Z^3*z^2 - d*Q.y*V.Z^3*z^3;
        w[0] = new JSField2(this.bn, d.multiply(V.y).subtract(n.multiply(V.x).multiply(V.z)).mod(this.bn.p));
        w[2] = Q.x.multiply(n.multiply(Vz3));
        w[3] = Q.y.multiply(p.subtract(d).multiply(Vz3));
        w[1] = w[4] = w[5] = this.E2.Fp2_0;
        return new JSField12(this.bn, w);
    };

    /*
     * 
     */
    JSPairing.reset = function () {
        this.addcount = 0;
        this.mulcount = 0;
        this.sqrcount = 0;
        this.modcount = 0;
        this.fpmcount = 0;
        JSField2.reset();
    };

    /*
     * 
     */
    JSPairing.update = function () {
        this.addcount += JSField2.getadd();
        this.mulcount += JSField2.getmul();
        this.sqrcount += JSField2.getsqr();
        this.modcount += JSField2.getmod();
        this.fpmcount += JSField2.getfpm();
    };

    /*
     * 
     */
    JSPairing.tate = function (P, Q) {
        var f = this.Fp12_1;
        P = P.normalize();
        Q = Q.normalize();
        if (!P.isZero() && !Q.isZero()) {
            var bn = this.E.bn;
            var V = P;
            for (var i = bn.n.bitLength() - 2; i >= 0; i--) {
                f = f.square().multiply(this.gl(V, V, Q));
                V = V.twice(1);
                if (bn.n.testBit(i)) {
                    f = f.multiply(this.gl(V, P, Q));
                    V = V.add(P);
                }
            }
            f = f.finExp();
        }
        return f;
    };

    /*
     * The eta (sometimes called twisted ate) pairing for points P and Q on BN curves E and E'.
     */
    JSPairing.eta = function (P, Q) {
        var f = this.Fp12_1;
        P = P.normalize();
        Q = Q.normalize();
        if (!P.isZero() && !Q.isZero()) {
            var bn = this.E.bn;
            var V = P;
            var ord = bn.rho; // the Tate pairing would have order bn.n instead of bn.rho
            for (var i = ord.bitLength() - 2; i >= 0; i--) {
                f = f.square().multiply(this.gl(V, V, Q));
                V = V.twice(1);
                if (ord.testBit(i)) {
                    f = f.multiply(this.gl(V, P, Q));
                    V = V.add(P);
                }
            }
            if (bn.u.signum() < 0) {
                // Aranha's trick:
                f = f.conjugate(3);
                //f = f.inverse(); // f = f_{6u+2,Q}
            }
            f = f.finExp();
        }
        return f;
    };

    /*
     * Improved implementation of the ate pairing
     */
    JSPairing.ate = function (Q, P) {
        var f = this.Fp12_1;
        P = P.normalize();
        Q = Q.normalize();
        if (!P.isZero() && !Q.isZero()) {
            var ord = this.bn.t.subtract(this.E.bn._1);
            var X = Q.x;
            var Y = Q.y;
            var Z = Q.z;
            var w = new Array(6);
            var start = ord.bitLength() - 2;
            for (var i = start; i >= 0; i--) {
                // Costello et al's double-and-line technique
                var A = X.square();
                var B = Y.square();
                var C = Z.square();
                if (this.bn.b === 3) {
                    var D = C.multiply(new Number(3*this.bn.b)).multiplyV();
                } else {
                    var D = C.multiply(new Number(3*this.bn.b)).divideV();
                }
                var F = Y.add(Z).square().subtract(B).subtract(C);
                if (i > 0) {
                    this.E = X.add(Y).square().subtract(A).subtract(B);
                    var G = D.multiply(new Number(3));
                    X = this.E.multiply(B.subtract(G));
                    Y = B.add(G).square().subtract(D.square().twice(2).multiply(new Number(3)));
                    Z = B.multiply(F).twice(2);
                }
                // line = L_10*x_P + L_01*y_P*z + L_00*z^3
                w[0] = F.multiply(P.y.negate()); // L_{0,1}
                w[1] = A.multiply(new Number(3)).multiply(P.x); // L_{1,0}
                w[3] = D.subtract(B); // L_{0,0}
                w[2] = w[4] = w[5] = this.E2.Fp2_0;
                var line = new JSField12(this.bn, w);
                if (i !== ord.bitLength() - 2) {
                    f = f.square().multiply(line);
                } else {
                    f = new JSField12(line);
                }
                if (ord.testBit(i)) {
                    // Costello et al's add-and-line technique
                    A = X.subtract(Z.multiply(Q.x)); B = Y.subtract(Z.multiply(Q.y));
                    //gADD = B*Q.x - A*Q.y - B*P.x + A*P.y;
                    w[0] = A.multiply(P.y); // L_{0,1}
                    w[1] = B.multiply(P.x.negate()); // L_{1,0}
                    w[3] = B.multiply(Q.x).subtract(A.multiply(Q.y)); // L_{0,0}
                    w[2] = w[4] = w[5] = this.E2.Fp2_0;
                    line = new JSField12(this.bn, w);
                    f = f.multiply(line);
                    C = A.square(); X = X.multiply(C); C = C.multiply(A);
                    D = B.square().multiply(Z).add(C).subtract(X.twice(1));
                    Y = B.multiply(X.subtract(D)).subtract(Y.multiply(C));
                    X = A.multiply(D);
                    Z = Z.multiply(C);
                }
            }
            f = f.finExp();
        }
        return f;
    };

    /*
     * Improved implementation of the optimal pairing
     */
    JSPairing.opt = function (Q, P) {
        var f = this.Fp12_1;
        P = P.normalize();
        Q = Q.normalize();
        this.reset();
        JSField2.countoff(true);
        if (!P.isZero() && !Q.isZero()) {
            var ord = this.bn.optOrd; // |6u+2|
            JSField2.counton(MillerLoop);
            var X = Q.x;
            var Y = Q.y;
            var Z = Q.z;
            var w = new Array(6);
            var start = ord.bitLength() - 2;
            for (var i = start; i >= 0; i--) {
                // Costello et al's double-and-line technique
                var A = X.square();
                var B = Y.square();
                var C = Z.square();
                if (this.bn.b === 3) {
                    var D = C.multiply(new Number(3*this.bn.b)).multiplyV();
                } else {
                    var D = C.multiply(new Number(3*this.bn.b)).divideV();
                }
                var F = Y.add(Z).square().subtract(B).subtract(C);
                if (i > 0) {
                    this.E = X.add(Y).square().subtract(A).subtract(B);
                    var G = D.multiply(new Number(3));
                    X = this.E.multiply(B.subtract(G));
                    Y = B.add(G).square().subtract(D.square().twice(2).multiply(new Number(3)));
                    Z = B.multiply(F).twice(2);
                }
                // line = L_10*x_P + L_01*y_P*z + L_00*z^3
                w[0] = F.multiply(P.y.negate());    // L_{0,1}
                w[1] = A.multiply(new Number(3)).multiply(P.x); // L_{1,0}
                w[3] = D.subtract(B);               // L_{0,0}
                w[2] = w[4] = w[5] = this.E2.Fp2_0;
                var line = new JSField12(this.bn, w);
                if (i !== ord.bitLength() - 2) {
                    f = f.square().multiply(line);
                } else {
                    f = new JSField12(line);
                }
                if (ord.testBit(i)) {
                    // Costello et al's add-and-line technique
                    A = X.subtract(Z.multiply(Q.x)); B = Y.subtract(Z.multiply(Q.y));
                    //gADD = B*Q.x - A*Q.y - B*P.x + A*P.y;
                    w[0] = A.multiply(P.y);                           // L_{0,1}
                    w[1] = B.multiply(P.x.negate());                  // L_{1,0}
                    w[3] = B.multiply(Q.x).subtract(A.multiply(Q.y)); // L_{0,0}
                    w[2] = w[4] = w[5] = this.E2.Fp2_0;
                    line = new JSField12(this.bn, w);
                    f = f.multiply(line);
                    C = A.square(); X = X.multiply(C); C = C.multiply(A);
                    D = B.square().multiply(Z).add(C).subtract(X.twice(1));
                    Y = B.multiply(X.subtract(D)).subtract(Y.multiply(C));
                    X = A.multiply(D);
                    Z = Z.multiply(C);
                    //JSField2.countoff(MillerLoop);
                }
            }
            // now T = [|6u+2|]Q and f = f_{|6u+2|,Q}
            if (this.bn.u.signum() < 0) {
                // Aranha's trick:
                f = f.conjugate(3);
                //f = f.inverse(); // f = f_{6u+2,Q}
            }
            // optimal pairing: f = f_{6u+2,Q}(P)*l_{Q3,-Q2}(P)*l_{-Q2+Q3,Q1}(P)*l_{Q1-Q2+Q3,[6u+2]Q}(P)
            var Q1 = Q.frobex(1);
            var Q2 = Q.frobex(2).negate();
            var Q3 = Q.frobex(3);
            // Costello et al's add-and-line technique
            X = Q2.x; Y = Q2.y; Z = Q2.z;
            A = X.subtract(Q3.x); B = Y.subtract(Q3.y);
            //A = X.subtract(Z.multiply(Q3.x)); B = Y.subtract(Z.multiply(Q3.y));
            //gADD = B*Q3.x - A*Q3.y - B*P.x + A*P.y;
            w[0] = A.multiply(P.y); // L_{0,1}
            w[1] = B.multiply(P.x.negate()); // L_{1,0}
            w[3] = B.multiply(Q3.x).subtract(A.multiply(Q3.y)); // L_{0,0}
            w[2] = w[4] = w[5] = this.E2.Fp2_0;
            line = new JSField12(this.bn, w);
            var line1 = new JSField12(line);

            C = A.square(); X = X.multiply(C); C = C.multiply(A);
            D = B.square().add(C).subtract(X.twice(1));
            Y = B.multiply(X.subtract(D)).subtract(Y.multiply(C));
            X = A.multiply(D);
            Z = C;

            // Costello et al's add-and-line technique
            //X = Q4.x; Y = Q4.y; Z = Q4.z;
            A = X.subtract(Z.multiply(Q1.x)); B = Y.subtract(Z.multiply(Q1.y));
            //gADD = B*Q1.x - A*Q1.y - B*P.x + A*P.y;
            w[0] = A.multiply(P.y); // L_{0,1}
            w[1] = B.multiply(P.x.negate()); // L_{1,0}
            w[3] = B.multiply(Q1.x).subtract(A.multiply(Q1.y)); // L_{0,0}
            w[2] = w[4] = w[5] = this.E2.Fp2_0;
            line = new JSField12(this.bn, w);
            var line2 = new JSField12(line);
            f = f.multiply(line1.multiply(line2));

            JSField2.countoff(MillerLoop);

            JSField2.counton(FinalExp);
            f = f.finExp();
            JSField2.countoff(FinalExp);
        }
        this.update();
        return f;
    };

    JSPairing.prototype.gl = JSPairing.gl;
    JSPairing.prototype.reset = JSPairing.reset;
    JSPairing.prototype.update = JSPairing.update;
    JSPairing.prototype.tate = JSPairing.tate;
    JSPairing.prototype.eta = JSPairing.eta;
    JSPairing.prototype.ate = JSPairing.ate;
    JSPairing.prototype.opt = JSPairing.opt;

    return JSPairing;
})();
