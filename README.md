# mtprng
Mersenne twister MT19937 pseudorandom number generator for haXe language

based on http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
and https://gist.github.com/banksean/300494  

## Install using haxelib
  
```sh
haxelib install mtprng
```

## Creating new MT instance

```haxe
var mt = new mtprng.MT( 42 ); // 42 is the PRNG seed
var mt = new mtprng.MT();			// use default seed (haxe.Timer.stamp())
```

## Methods

* ``randomUInt()`` -- returns pseudorandom ``UInt`` in range [0,2^32)
* ``randomInt()`` -- returns pseudorandom ``Int`` in range [-2^31,2^31)
* ``randomFloat()`` -- returns pseudorandom ``Float`` with 53-bits precision
* ``randomFloat32()`` -- returns pseudorandom ``Float`` with 32-bits precision

## Using static instance

```haxe
var x = mtprng.MT.instance.randomInt();
var y = mtprng.MT.instance.randomFloat();
```

## Platforms

* cpp, java, js, python, neko -- works, same results
* as3, swf -- not tested
* cs -- not compiling
* php -- wrong results
