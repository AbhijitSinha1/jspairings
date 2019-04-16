const JSTests = require('./examples/JSTests');
const jsTests = new JSTests();

let multiplyFp = jsTests.multiplyFp();
let multiplyFp2 = jsTests.multiplyFp2();
let precomputedMultiplyFp = jsTests.precomputedMultiplyFp();
let addFp = jsTests.addFp();
let pairing = jsTests.pairing();

console.log({multiplyFp, multiplyFp2, precomputedMultiplyFp, addFp, pairing })

// jsTests.test1();
// jsTests.test2();
// jsTests.test3();
// jsTests.test4();