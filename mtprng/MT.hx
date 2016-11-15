/* 
	 Haxe implementation of MT19937 pseudorandom number generator.
	 (see authors page http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html)
	 Working on cpp,js,java,python,neko. Not tested on as3,swf.
	 Not working on cs (not compiling) and php (wrong results).

	 Git repository https://github.com/iskolbin/mt
	 Written by Ilya Kolbin (iskolbin@gmail.com)
*/

/* 
	 JavaScript and Python conditional compilation based on
	 https://gist.github.com/banksean/300494 by Sean McCullough (banksean@gmail.com)
	 
	 see also npm package https://github.com/boo1ean/mersenne-twister
*/

/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.
 
   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).
 
   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          
 
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:
 
     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
 
     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
 
     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.
 
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 
   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

package mtprng;

class MT {
	public static inline var N = 624;
	public static inline var M = 397;
	public static inline var MATRIX_A: UInt   = 0x9908b0df;
	public static inline var UPPER_MASK: UInt = 0x80000000;
	public static inline var LOWER_MASK: UInt = 0x7fffffff;
	public static inline var FF_MASK: UInt    = 0xffffffff;
	public static inline var MULT: UInt       = 1812433253;
	public static inline var ZERO: UInt 			= 0;
	public static inline var DEFAULT: UInt 		= 5489;
	public static inline var ONE: UInt 				= 1;
	public static inline var TEMPER_1: UInt 	= 0x9d2c5680;
	public static inline var TEMPER_2: UInt 	= 0xefc60000;
	public static inline var AH_MASK: UInt		= 0xffff0000;
	public static inline var AL_MASK: UInt		= 0x0000ffff;

	public static var instance(default,null) = new MT();

	static var mag01 = [ZERO, MATRIX_A];

	public var mt(default,null): Array<UInt>;
	public var mti(default,null) = 0;

	public function new( ?s: UInt ) {
		mt = new Array<UInt>();
		for ( i in 0...N ) mt[i] = 0;
		init( s == null ? Std.int( haxe.Timer.stamp()) : s );
	}

	public function init( s: UInt ) {
#if (js||python||lua)
		mt[0] = s >>> 0;
#else
		mt[0] = s & FF_MASK;
#end
		for ( j in 1...N ) {
			var s = (mt[j-1] ^ (mt[j-1] >> 30));
#if (js||python||lua)
			mt[j] = ((((((s & AH_MASK) >>> 16) * MULT) << 16) + (s & AL_MASK) * MULT) + j) >>> 0;
#else
			mt[j] = (MULT * s + j) & FF_MASK;
#end
		}
		mti = N;
	}

	static inline var SAR: UInt = 19650218;
	static inline var SBR: UInt = 1664525;
	static inline var SCR: UInt = 1566083941;

	// TODO fix for js,python,lua
	public static function makeFromArray( initKey: Array<UInt> ) {
		var i = 1;
	 	var j = 0;
		var k = N > initKey.length ? N : initKey.length;
		var self = new MT( SAR );
		var mt = self.mt;
		while ( k > 0 ) {
#if (js||python||lua)
			var s = mt[i-1] ^ (mt[i-1] >> 30);
			mt[i] = (mt[i]^(((((s & AH_MASK) >>> 16) * SBR) << 16) + (s & AL_MASK) * SBR) + initKey[j] + j) >>> 0; /* non linear */
#else
			var s = mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * SBR);
			mt[i] = (s + initKey[j] + j) & FF_MASK;
#end
			i++;
			j++;
			if ( i >= N ) {
				mt[0] = mt[N-1];
				i = 1;
			}
			if (j >= initKey.length ) {
				j = 0;
			}
			k--;
		}
		var k = N-1;
		while ( k > 0 ) {
#if (js||python||lua)
			var s = mt[i-1] ^ (mt[i-1] >> 30);
			mt[i] = (mt[i]^(((((s & AH_MASK) >>> 16) * SCR) << 16) + (s & AL_MASK) * SCR) - i) >>> 0; /* non linear */
#else
			var s = mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * SCR);
			mt[i] = (s - i) & FF_MASK;
#end
			i++;
			if ( i>= N ) {
				mt[0] = mt[N-1];
				i=1;
			}
			k--;
		}

		mt[0] = UPPER_MASK; /* MSB is 1; assuring non-zero initial array */
		return self;
	}

	public function randomUInt(): UInt {
		var mt = this.mt;
		var y: UInt;
		var mag01 = MT.mag01;

		if ( mti >= N ) { 		/* generate N words at one time */
			if ( mti == N+1 )   /* if init() has not been called, */
				init( DEFAULT ); 	/* a default initial seed is used */

			for ( kk in 0...N-M ) {
				y = ( mt[kk] & UPPER_MASK ) | ( mt[kk+1] & LOWER_MASK );
				mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & ONE];
			}
			for ( kk in N-M...N-1 ) {
				y = ( mt[kk] & UPPER_MASK ) | ( mt[kk+1] & LOWER_MASK );
				mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & ONE];
			}
			y = ( mt[N-1] & UPPER_MASK ) | ( mt[0] & LOWER_MASK );
			mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & ONE];

			mti = 0;
		}
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & TEMPER_1;
    y ^= (y << 15) & TEMPER_2;
    y ^= (y >> 18);

#if python
		return y >>> 0;
#else
    return y;
#end
	}	

	public inline function random( limit: Int ) {
		return randomUInt() % limit;
	}

	public inline function randomInt(): Int {
		var x: Int = randomUInt();
		return x;
	}

	public function randomFloat(): Float {
		var a = randomUInt() >> 5;
		var b = randomUInt() >> 6;
		var a_: Float = a;
		var b_: Float = b;
		return (a_ * 67108864.0 + b_) * (1.0 / 9007199254740992.0);
	}

	public function randomFloat32(): Float {
		var a = randomUInt();
		var a_: Float = a;
		return a_ / 4294967296.0;
	}
}
