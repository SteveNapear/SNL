using System;

namespace SNL {
	#region oldfft
	//////////////////////////////////////////////////////////////////////////
	/*
		This class implements a Fast Fourier Transform.
		When the object is constructed, the single integer value denotes
		the Log base 2 of the number of elements.
		At the instantiation, the roots of 1 are generated, along with
		a bit reversed index table.
		Calls to the doFFT method use these two tables.
	
	
	public class CFFT
	{
		//................................................................
		public CFFT(uint NL2) { Initialize(NL2); }  //constructor

		//................................................................
		private static int pow2n(int i) { return (1 << i); }    //2**n
		private static int pow2n(uint u) { return (1 << (int)u); }
		//		private uint UMask(int ix){uint m=1; return(m<<ix);}

		private uint[] BR = null;               //bit reversed index table
		private uint nl2;                       //N log2
		private uint[] MT = null;                   //mask table
		private CSNplex[] Roots = null;         //complex unit circle
												//================================================================
		private void GenBR(uint NL2)            //build bit reversed index table and roots of 1
		{
			int ix, i, len; uint mask; double inc;
			MT = new uint[NL2];                 //allocate mask table
			for (ix = 0; ix < MT.Length; ix++) MT[ix] = (uint)pow2n((int)(NL2 - ix - 1));

			nl2 = NL2;
			len = 1 << (int)nl2;                //len = 2**nl2
			BR = new uint[len];                 //allocate bit reversed index table (uint array)
			Roots = new CSNplex[Len];           //allocate roots of 1 array (Complex valued array)
			inc = Math.PI * 2.0 / Len;
			//fill arrays with appropriate values
			for (ix = 0; ix < BR.Length; ix++)
			{
				mask = 0;                       //reset mask
				for (i = 0; i < nl2; i++) { if ((pow2n(i) & ix) > 0) mask |= MT[i]; }
				BR[ix] = mask;
				Roots[ix] = new CSNplex(Math.Cos(ix * inc), Math.Sin(ix * inc));    //generate roots of 1
			}
		}
		//================================================================
		private void Initialize(uint NL2) { GenBR(NL2); }   //generate roots of 1 and bit reversed index table
															//.................................................................
		public uint BRMask(int ix) { return (ix < MT.Length) ? MT[ix] : 0; }
		public int Len { get { return BR.Length; } }
		public int L2N { get { return (int)nl2; } }
		public uint BRIT(int ix) { return (ix < BR.Length) ? BR[ix] : 0; }
		public CSNplex Root(int ix) { return (ix < Roots.Length) ? Roots[ix] : new CSNplex(0, 0); }
		 
			---FFT from code by Steve Napear, 1971
		C	A is the array to be transformed (real array)
		C	C is the complex array that is returned
		C	N is the  size of the array A (Must be a power of 2)
		C	L2 is the Log2 of N

			Dimension A(N),CS(64),SN(64)
			Integer PT,GO
			Complex C(N),W,T1,T2,VAC
			Common/WCB/W(512)
			Common/LCB/PT(1024)
			Data GO /0/
		C
		C	Pickup Log(2)N
		C
			L2 = LOG2(N)
			if(L2>0) goto 8
			Print -- (Invalid array size ...)
			return
		8	continue
			Call BitR(N,L2) - generate reverse bit indexing array

		C	---Setup n roots of 1 complex table W
		C	---generating the complex array W containing the n roots of 1
			N4 = N/4
			N4P2 = N4 + 2
			N4M2 = N4 - 2
			N2P2 = N/2 + 2;
			PI = 3.14159265897932
			VAL = PI*2.0/N
			CS(1) = 1.0
			CS(N4+1) = SN(1) = 0.0
			SN(N4+1) = - 1.0
			DO 9 i=2,N4
				R = i-1
				S1 = Sin(R*VAL)
				CS(N4P2-i) = S1
				CS(N4+i) = SN(i) = SN(N2P2-i) = -S1
		9	continue

			DO 905 i=1,64
		905	W(i) = CMPLX(CS(i),SN(i))

			GO=N
		91	continue
		C
		C	--Set into bit reversed addressed order and set complex
		C
			RN=N
			Do 10 i=1,N
			J = PT(i)+1
		10	C(J) = CMPLX(A(i)/RN,0,0)
		C	
		C	---DO FFT
		C
			DO 11 ii=1,L2
			M1 = 2**(ii-1)
			M2 = M1*2
			M3 = N-M2+1
			NOM = N/M2
			Do 11 I=1,M3,M2
			DO 11 J=1,M1
			L = J*I-1
			T1 = C(L)
			T2 = C(L+M1)*W(NOM*(J-1)+1)
			C(L) = T1+T2
			C(L+M1) = T1-T2
		11	continue
			return
			end

		
	//================================================================
	//This algorithm was developed from code generated in 1971 
	//at Naval Intelligence Command by Steve Napear
	public CSNplex[] Dofft(double[] A)
		{
			if ((A == null) || (A.Length != Len)) return null;  //must match length of Root and BR tables
			CSNplex[] C = new CSNplex[Len];
			for (int ii = 0; ii < Len; ii++) C[BR[ii]] = new CSNplex(A[ii] / Len, 0);   //place in Bit reversed order
			return Dofft(C);
		}

		//================================================================
		public CSNplex[] Dofft(CSNplex[] C)
		{
			int i, ii, j, l, m1, m2, m3, nom; CSNplex t1, t2;

			//do fft					
			for (ii = 0; ii < L2N; ii++)
			{      //DO 11 II=1,L2  where L2 is the log base 2 of N
				m1 = (int)pow2n(ii);    //M1 = 2**(ii-1)
				m2 = m1 << 1;               //M2 = M1*2
				m3 = Len - m2;          //M3 = N-M2+1
				nom = Len / m2;         //NOM = N/M2
				for (i = 0; i < m3; i += m2)
				{  //DO 11 I=1,M3,M2
					for (j = 0; j < m1; j++)
					{  //DO 11 J=1,M1
						l = j + i;      //L=J+I-1
						t1 = C[l];      //T1 = C[l];
						t2 = C[l + m1] * Roots[nom * j];    //T2 = C(L+M1)*W(NOM*(J-1)+1];
						C[l] = t1 + t2; //C[L] = T1+T2
						C[l + m1] = t1 - t2;//C[L+M1] = T1-T2;
					}
				}
			}
			return C;
		}
	}
	*/

	#endregion
	//////////////////////////////////////////////////////////////////////
	/**  

015      * Performs an in-place complex FFT.  
016      *  
017      * Released under the MIT License  
019      * Copyright (c) 2010 Gerald T. Beauregard  
020      *  
021      * Permission is hereby granted, free of charge, to any person obtaining a copy  
022      * of this software and associated documentation files (the "Software"), to  
023      * deal in the Software without restriction, including without limitation the  
024      * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or  
025      * sell copies of the Software, and to permit persons to whom the Software is  
026      * furnished to do so, subject to the following conditions:  
027      *  
028      * The above copyright notice and this permission notice shall be included in  
029      * all copies or substantial portions of the Software.  
031      * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  
032      * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,  
033      * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE  
034      * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER  
035      * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING  
036      * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS  
037      * IN THE SOFTWARE.  
038      */
	////////////////////////////////////////////////////////////////////////////////////
	[Serializable] public class CFFT2
	{

		// Element for linked list in which we store the  
		// input/output data. We use a linked list because  
		// for sequential access it's faster than array index.  
		[Serializable]  class FFTElement
		{

			public double re = 0.0;     // Real component  
			public double im = 0.0;     // Imaginary component  
			public FFTElement next;     // Next element in linked list  
			public uint revTgt;         // Target position post bit-reversal  
		}

		private uint m_logN = 0;        // log2 of FFT size  
		private uint m_N = 0;           // FFT size  
		private FFTElement[] m_X;       // Vector of linked list elements  

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		public CFFT2(uint log2N) { init(log2N); }  //Constructor
		 //Initialize class to perform FFT of specified size.  
		//@param   logN    Log2 of FFT length. e.g. for 512 pt FFT, logN = 9.  
		//================================================================
		public uint nFFTsize { get { return m_N; } }
		public int Length { get => (int)m_N; }
		private void init(uint logN)    //initialize class for 2**logN elements  
		{
			m_logN = logN;
			m_N = (uint)(1 << (int)m_logN);
			// Allocate elements for linked list of complex numbers.  
			m_X = new FFTElement[m_N];
			for (uint k = 0; k < m_N; k++) m_X[k] = new FFTElement();
			// Set up "next" pointers.  
			for (uint k = 0; k < m_N - 1; k++) m_X[k].next = m_X[k + 1];
			// Specify target for bit reversal re-ordering.  
			for (uint k = 0; k < m_N; k++) m_X[k].revTgt = BitReverse(k, logN);
		}

		/**  
089          * Performs in-place complex FFT.  
090          *  
091          * @param   xRe     Real part of input/output  
092          * @param   xIm     Imaginary part of input/output  
093          * @param   inverse If true, do an inverse FFT  
094          */
		//================================================================
		public CSNplex[] run(CSNplex[] zA,bool inverse)	//returns a new CSNplex array with FFT values
		{
			int ix; CSNplex z;
			double[] rA = new double[zA.Length];
			double[] iA = new double[zA.Length];
			for(ix=0;ix<zA.Length;ix++){rA[ix] = zA[ix].r; iA[ix] = zA[ix].i;}
			if(!run(rA,iA,inverse)) return null;
			CSNplex[] zB = new CSNplex[zA.Length];
			for(ix=0;ix<zA.Length;ix++){z = new CSNplex(rA[ix],iA[ix]); zB[ix] = z;}
			return zB;
		}
		//================================================================
		public bool run(double[] xRe, double[] xIm, bool inverse)
		{

			if ((xRe.Length != m_N) || (xIm.Length != m_N)) return false;
			uint numFlies = m_N >> 1; // Number of butterflies per sub-FFT  
			uint span = m_N >> 1;     // Width of the butterfly  
			uint spacing = m_N;         // Distance between start of sub-FFTs  
			uint wIndexStep = 1;        // Increment for twiddle table index  
										// Copy data into linked complex number objects  
										// If it's an IFFT, we divide by N while we're at it  

			FFTElement x = m_X[0];
			uint k = 0;
			double scale = inverse ? 1.0 / m_N : 1.0;
			while (x != null)
			{
				x.re = scale * xRe[k];
				x.im = scale * xIm[k];
				x = x.next;
				k++;
			}

			// For each stage of the FFT  
			for (uint stage = 0; stage < m_logN; stage++)
			{
				/* 
				 Compute a multiplier factor for the "twiddle factors".  
				 The twiddle factors are complex unit vectors spaced at  
				 regular angular intervals. The angle by which the twiddle  
				 factor advances depends on the FFT stage. In many FFT  
				 implementations the twiddle factors are cached, but because  
				 array lookup is relatively slow in C#, it's just  
				 as fast to compute them on the fly.  
			   */
				double wAngleInc = wIndexStep * 2.0 * Math.PI / m_N;
				if (inverse) wAngleInc *= -1;

				double wMulRe = Math.Cos(wAngleInc);
				double wMulIm = Math.Sin(wAngleInc);
				for (uint start = 0; start < m_N; start += spacing)
				{
					FFTElement xTop = m_X[start];
					FFTElement xBot = m_X[start + span];
					double wRe = 1.0;
					double wIm = 0.0;
					// For each butterfly in this stage  
					for (uint flyCount = 0; flyCount < numFlies; ++flyCount)
					{
						// Get the top & bottom values  
						double xTopRe = xTop.re;
						double xTopIm = xTop.im;
						double xBotRe = xBot.re;
						double xBotIm = xBot.im;
						// Top branch of butterfly has addition  
						xTop.re = xTopRe + xBotRe;
						xTop.im = xTopIm + xBotIm;
						// Bottom branch of butterly has subtraction,  
						// followed by multiplication by twiddle factor  
						xBotRe = xTopRe - xBotRe;
						xBotIm = xTopIm - xBotIm;
						xBot.re = xBotRe * wRe - xBotIm * wIm;
						xBot.im = xBotRe * wIm + xBotIm * wRe;

						// Advance butterfly to next top & bottom positions  
						xTop = xTop.next;
						xBot = xBot.next;
						// Update the twiddle factor, via complex multiply  
						// by unit vector with the appropriate angle  
						// (wRe + j wIm) = (wRe + j wIm) x (wMulRe + j wMulIm)  

						double tRe = wRe;
						wRe = wRe * wMulRe - wIm * wMulIm;
						wIm = tRe * wMulIm + wIm * wMulRe;
					}
				}
				numFlies >>= 1;   // Divide by 2 by right shift  
				span >>= 1;
				spacing >>= 1;
				wIndexStep <<= 1;     // Multiply by 2 by left shift  
			}

			// The algorithm leaves the result in a scrambled order.  
			// Unscramble while copying values from the complex  
			// linked list elements back to the input/output vectors.  

			x = m_X[0];
			while (x != null)
			{
				uint target = x.revTgt;
				xRe[target] = x.re;
				xIm[target] = x.im;
				x = x.next;
			}
			return true;
		}

		/**  
		 * Do bit reversal of specified number of places of an int  
		* For example, 1101 bit-reversed is 1011  
		 * @param   x       Number to be bit-reverse.  
		 * @param   numBits Number of bits in the number.  
		 */

		private uint BitReverse(uint x, uint numBits)
		{
			uint y = 0;
			for (uint i = 0; i < numBits; i++)
			{
				y <<= 1;
				y |= x & 0x0001;
				x >>= 1;
			}
			return y;
		}
	}
	//////////////////////////////////////////////////////////////////////
	[Serializable]
	public class CSNplex
	{
		private static readonly double π = Math.PI;
		private static char degsign = (char)176;
		private static int CSNPLEX_VERSION = 2;
		private static string sPrec = "F6";
		private double real;
		private double imag;
		//constructors
		public CSNplex() { real = 0; imag = 0; }
		public CSNplex(double r) { real = r; imag = 0; }
		public CSNplex(double rv, double iv) { real = rv; imag = iv; }
		public CSNplex(CSNplex rhs) { real = rhs.real; imag = rhs.imag; }   //copy constructor
																			//=================================================================																	//boolean operations
		public static bool operator ==(CSNplex lhs, CSNplex rhs)
		{   //comparision operator
			if (lhs.real != rhs.real) return false;
			if (lhs.imag != rhs.imag) return false;
			return true;
		}
		//================================================================
		public override int GetHashCode() { return base.GetHashCode(); }
		//================================================================
		public override bool Equals(object obj)
		{
			if (!(obj is CSNplex)) return false;
			return this == (CSNplex)obj;
		}
		public static bool operator !=(CSNplex lhs, CSNplex rhs) { return !(lhs == rhs); }  //comparision operator
																							//conversion operators
		public static implicit operator CSNplex(int v) { return new CSNplex((double)v, 0.0); }
		public static implicit operator CSNplex(float v) { return new CSNplex((double)v, 0.0); }
		public static implicit operator CSNplex(double v) { return new CSNplex(v, 0.0); }
		//arithmetic operations
		//public static CSNplex operator+(CSNplex lhs, double d) { return new CSNplex(lhs.r+d, lhs.i); }
		//	public static CSNplex operator-(CSNplex lhs, double d) { return new CSNplex(lhs.r-d, lhs.i); }
		public static CSNplex operator -(CSNplex rhs) { return new CSNplex(-rhs.r, -rhs.i); }
		public static CSNplex operator +(CSNplex lhs, CSNplex rhs)
		{
			return new CSNplex(lhs.real + rhs.real, lhs.imag + rhs.imag);
		}
		public static CSNplex operator -(CSNplex lhs, CSNplex rhs)
		{
			return new CSNplex(lhs.real - rhs.real, lhs.imag - rhs.imag);
		}
		public static CSNplex operator *(CSNplex lhs, double m) { return new CSNplex(lhs.real * m, lhs.imag * m); }
		public static CSNplex operator *(CSNplex lhs, CSNplex rhs)
		{   //multiplication by a complex number
			return new CSNplex(lhs.real * rhs.real - lhs.imag * rhs.imag, lhs.imag * rhs.real + lhs.real * rhs.imag);
		}
		public static CSNplex operator /(CSNplex num, double dem) { return (dem != 0.0) ? (new CSNplex(num.real / dem, num.imag / dem)) : CSNplex.Invalid; }
		public static CSNplex operator /(CSNplex lhs, CSNplex rhs)
		{   //division by a complex number
			double dem = rhs.real * rhs.real + rhs.imag * rhs.imag;
			double rnum = lhs.real * rhs.real + lhs.imag * rhs.imag;
			double inum = lhs.imag * rhs.real - lhs.real * rhs.imag;
			return (dem != 0.0) ? (new CSNplex(rnum / dem, inum / dem)) : CSNplex.Invalid;
		}
		//standard functions
		public double r { get { return real; } set { real = value; } }
		public double i { get { return imag; } set { imag = value; } }
		public double Real { get { return real; } set { real = value; } }
		public double Imag { get { return imag; } set { imag = value; } }
		public double Amp { get { return Math.Sqrt(real * real + imag * imag); } }
		public double Phase { get { return Math.Atan2(imag, real); } }
		public string sPhase { get { return (180.0 * Phase / π).ToString("f0") + degsign; } }
		public double Arg { get { return Phase; } }
		public bool isValid { get { if (double.IsNaN(real) || double.IsNaN(imag)) return false; return true; } }
		//====================================================================
		public override string ToString()
		{
			if (!isValid) return "Invalid";
			string s = '(' + real.ToString(sPrec) + ", " + imag.ToString(sPrec) + "i)" +
			", [" + Amp.ToString(sPrec) + ", " + (Phase * 180 / Math.PI).ToString("f2") + degsign + ']';
			return s;
		}
		//statics
		public static CSNplex Invalid { get { return new CSNplex(double.NaN, 0); } }
		public static int Version { get { return CSNPLEX_VERSION; } }
		//.....................................................................
		public static int Log2N(int N)
		{
			int log = 0;
			do { N >>= 1; if (N == 0) break; log += 1; } while (true);
			return log;
		}
		//.....................................................................
		private static CSNplex ln(CSNplex z)
		{
			double rad = Math.Log(z.Amp);
			double theta = Math.Atan2(z.i, z.r);
			return new CSNplex(rad, theta);
		}
		//.....................................................................
		public static CSNplex Cln(CSNplex z)    //Complex ln
		{
			double amp, theta, rad, r = z.r, im = z.i;
			if (z.r == 0)
			{                                       //either pi or -pi
				theta = (z.i < 0) ? -Math.PI : Math.PI;
				return new CSNplex(0, theta);
			}
			amp = z.Amp;
			if (amp == 0) return CSNplex.Invalid;   //Can't take the log of 0
			rad = Math.Log(z.Amp);
			theta = Math.Atan2(z.i, z.r);
			return new CSNplex(rad, theta);
		}
		//.....................................................................
		public static CSNplex Cexp(CSNplex z) //Complex e**z
		{
			double c1 = Math.Exp(z.r), c2 = z.i;
			return new CSNplex(c1 * Math.Cos(c2), c1 * Math.Sin(c2));
		}
		//.....................................................................
		public static CSNplex Cpower(CSNplex z, int n) //Complex z**n
		{
			CSNplex ans = 1, v = z;
			do
			{
				if ((n & 1) > 0) ans *= v;
				v *= v;
				n >>= 1;
			} while (n > 0);
			return ans;

		}
		//.....................................................................
		private static CSNplex CSin1(CSNplex z) //test module
		{
			int nc = 0, n = 0, n2, sign;
			CSNplex sum = 0, term, z1;
			double small = 0.0000001;
			do
			{
				if (!z.isValid) return CSNplex.Invalid;
				sign = ((n & 1) > 0) ? -1 : 1;
				n2 = 2 * n + 1;
				z1 = CSNplex.Cpower(z, n2);
				double fac = Factorial(n2);
				term = z1 / fac;
				sum += term;
				n += 1;
				nc += 1;
			} while (term.Amp > small);
			return sum;
		}
		//.....................................................................
		public static double Factorial(int n) //n!
		{
			if (n <= 0) return double.NaN;
			double v = 1.0;
			do
			{
				v *= (double)n;
				n -= 1;
			} while (n > 0);
			return v;
		}
		//.....................................................................
		public static CSNplex CSin(CSNplex z) //Complex sin(z)
		{
			CSNplex ans, num, zi = new CSNplex(-z.i, z.r), zmi = new CSNplex(z.i, -z.r), den = new CSNplex(0, 2);
			num = CSNplex.Cexp(zi) - CSNplex.Cexp(zmi);
			ans = num / den;
			return ans;
		}
		//.....................................................................
		public static CSNplex CCos(CSNplex z) //Complex cos(z)
		{
			CSNplex ans, num, zi = new CSNplex(-z.i, z.r), zmi = new CSNplex(z.i, -z.r), den = new CSNplex(0, 2);
			num = CSNplex.Cexp(zi) + CSNplex.Cexp(zmi);
			ans = num / 2;
			return ans;
		}
		//.....................................................................
		public static CSNplex CTan(CSNplex z)
		{
			CSNplex sin = CSin(z), cos = CCos(z);
			return (cos.isValid) ? sin / cos : CSNplex.Invalid;
		}
		//special functions
	};
    // ///////////////////////////////////////////////////////////////////
    public class CKalman
    {
        private double _oldest;
        private double _est;
        private double _olderror;
        private double _kg;

        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        public CKalman(double InitialEst, double InitialErr)
        {
            _oldest = InitialEst;
            _olderror = InitialErr;
        }
        //================================================================
        private double KG => _oldest / (_oldest + _olderror);
        public double Mea(double newmea)
        {
            _kg = KG;
            _est = _oldest + KG * (newmea - _oldest);
            _olderror = (1.0 - _kg) * _oldest;
            _oldest = _est;
            double delta = _est - newmea;
            return _est;
        }
    }

}