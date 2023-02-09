
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using System.Windows;
using System.Windows.Forms;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace SNL {

	//////////////////////////////////////////////////////////////////////////
	[Serializable]
	public partial class CSNC {
		// //////////////////////////////////////////////////////////////////
		[Serializable]
		public class CRGB/*Color in RGB*/{
			private byte a;
			private byte r;
			private byte g;
			private byte b;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			public CRGB(byte R, byte G, byte B) { a = 255; r = R; g = G; b = B; }
			public CRGB(System.Windows.Media.Color C) { a = C.A; r = C.R; g = C.G; b = C.B; }
			public byte A { get => a; set => a = value; }
			public byte R { get => r; set => r = value; }
			public byte G { get => g; set => g = value; }
			public byte B { get => b; set => b = value; }
			public System.Windows.Media.Color Color {
				get {
					System.Windows.Media.Color c = new Color();
					c.A = 255; c.R = r; c.G = g; c.B = b;
					return c;
				}
			}
			public static double operator -(CRGB c1, CRGB c2) { return CSNC.ColorDistance(c1.Color, c2.Color); }
			public static bool operator ==(CRGB O1, CRGB O2) { return (O1.a == O2.a && O1.r == O2.r && O1.g == O2.g && O1.b == O2.b); }
			public static bool operator !=(CRGB O1, CRGB O2) { return (O1.a != O2.a || O1.r != O2.r || O1.g != O2.g || O1.b != O2.b); }
			public override int GetHashCode() => base.GetHashCode();
			public override bool Equals(object obj) {
				if (obj == null || GetType() != obj.GetType()) return false;
				return (this == (CRGB)obj);
			}
			public override string ToString() {
				return '[' + r.ToString("d0") + ',' + g.ToString("d0") + ',' + b.ToString("d0") + ']';
			}
			public CLab LabD50 => CRGB2CLab(R, G, B, CXYZ.D50);
			public CLab LabD65 => CRGB2CLab(R, G, B, CXYZ.D65);
			public static string sColor(System.Drawing.Color c) => '[' + c.R.ToString("###") + ',' + c.G.ToString("###") + ',' + c.B.ToString("###") + ']';
			public static string sColor(System.Windows.Media.Color c) => '[' + c.R.ToString("###") + ',' + c.G.ToString("###") + ',' + c.B.ToString("###") + ']';
		}
		// /////////////////////////////////////////////////////////////
		[Serializable]
		public class CLab {
			//use structure
			private double l;
			private double a;
			private double b;
			private System.Windows.Media.Color rgb;
			private System.Windows.Media.Color newrgb;

			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			public CLab() {
				l = 0.0; a = 0.0; b = 0.0;
				rgb = System.Windows.Media.Colors.Black;
				newrgb = System.Windows.Media.Colors.White;
			}
			public CLab(double L, double A, double B) {
				l = L; a = A; b = B;
				rgb = System.Windows.Media.Colors.Black;
				newrgb = System.Windows.Media.Colors.White;
			}
			public CLab(byte red, byte green, byte blue, CXYZ icc) {
				CLab temp = CRGB2CLab(red, green, blue, icc);
				l = temp.L; a = temp.A; b = temp.B;
				rgb = System.Windows.Media.Colors.Black;
				rgb.R = red; rgb.G = green; rgb.B = blue;
				newrgb = System.Windows.Media.Colors.White;
			}
			//============================================================
			public System.Windows.Media.Color Color() { return Color(CXYZ.D50); }
			public System.Windows.Media.Color Color(bool bD50) { return Color(bD50 ? CXYZ.D50 : CXYZ.D65); }
			public System.Windows.Media.Color Color(CXYZ icc) {
				CXYZ xyz = CLab2CXYZ(l, a, b, icc);
				return xyz.Color;
			}
			//============================================================
			public static bool operator ==(CLab L1, CLab L2) { return (L1.L == L2.L && L1.A == L2.A && L1.B == L2.B); }
			public static bool operator !=(CLab L1, CLab L2) { return (L1.L != L2.L || L1.A != L2.A || L1.B != L2.B); }
			public double L { get => l; set => l = value; }
			public double A { get => a; set => a = value; }
			public double B { get => b; set => b = value; }
			public System.Windows.Media.Color RGB { get => rgb; set => rgb = value; }
			public System.Windows.Media.Color NewRGB { get => newrgb; set => newrgb = value; }
			public string sRGB { get => sColor(rgb); }
			//============================================================
			public override string ToString() {
				return "L: " + l.ToString("f1") + ", a: " + a.ToString("f1") + ", b: " + b.ToString("f1");
			}
			//============================================================
			public override int GetHashCode() => base.GetHashCode();
			public override bool Equals(Object obj) {
				if (obj == null || GetType() != obj.GetType()) return false;
				return (this == (CLab)obj);
			}
			//=============================================================
			//			public CRGB RGB50=>Lab2RGB(this,CXYZ.D50); 
		}
		// /////////////////////////////////////////////////////////////
		[Serializable]
		public class ClabList : List<CLab>/*List<CLab>*/ {
			private string spath = " ";
			private int nrows, ncols;
			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			public ClabList(int Capacity, int nR, int nC) : base(Capacity) {
				nrows = nR; ncols = nC;
			}
			//============================================================
			public string sPath { get => spath; set => spath = value; }
			public int nRows { get => nrows; }
			public int nCols { get => ncols; }
		}
		// /////////////////////////////////////////////////////////////
		[Serializable]
		public class CXYZ /*XYZ values*/{
			public static readonly CXYZ D65 = new CXYZ(0.9505, 1.0, 1.0890);
			public static readonly CXYZ D50 = new CXYZ(0.9642, 1.0000, 0.8251);
			private double x;
			private double y;
			private double z;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			public CXYZ() { x = 0.0; y = 0.0; z = 0.0; }
			public static bool operator ==(CXYZ item1, CXYZ item2) {
				return (item1.X == item2.X && item1.Y == item2.Y && item1.Z == item2.Z);
			}
			public static bool operator !=(CXYZ item1, CXYZ item2) {
				return (item1.X != item2.X || item1.Y != item2.Y || item1.Z != item2.Z);
			}
			public double X { get => x; set { this.x = (value > 0.9505) ? 0.9505 : ((value < 0) ? 0 : value); } }
			public double Y { get => y; set { y = (value > 1.0) ? 1.0 : (value < 0.0) ? 0.0 : value; } }
			public double Z { get => z; set { z = (value > 1.089) ? 1.089 : (value < 0.0) ? 0.0 : value; } }
			public override int GetHashCode() => base.GetHashCode();
			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			public CXYZ(double X, double Y, double Z) {
				x = (X > 0.9505) ? 0.9505 : ((X < 0) ? 0 : X);
				y = (Y > 1.0) ? 1.0 : ((Y < 0) ? 0 : Y);
				z = (Z > 1.089) ? 1.089 : ((Z < 0) ? 0 : Z);
			}
			//============================================================
			public override string ToString() {
				return "(" + x.ToString("f3") + ',' + y.ToString("f3") + ',' + z.ToString("f3") + ')';
			}
			//============================================================
			public override bool Equals(Object obj) {
				if (obj == null || GetType() != obj.GetType()) return false;
				return (this == (CXYZ)obj);
			}
			//================================================================
			public System.Windows.Media.Color Color { get { return GetColor(x, Y, Z); } }
			//................................................................
			public static System.Windows.Media.Color GetColor(double X, double Y, double Z) {
				double[] Cn = new double[3];
				Cn[0] = X * 3.2410 - Y * 1.5374 - Z * 0.4986; // red
				Cn[1] = -X * 0.9692 + Y * 1.8760 + Z * 0.0416; // green
				Cn[2] = X * 0.0556 - Y * 0.2040 + Z * 1.0570; // blue

				for (int i = 0; i < 3; i++) {
					Cn[i] = (Cn[i] <= 0.0031308) ? 12.92 * Cn[i] : (1 + 0.055) * Math.Pow(Cn[i], (1.0 / 2.4)) - 0.055;
				}
				System.Windows.Media.Color c = new System.Windows.Media.Color();
				c.R = (byte)(Cn[0] * 255.0 + 0.5);
				c.G = (byte)(Cn[1] * 255.0 + 0.5);
				c.B = (byte)(Cn[2] * 255.0 + 0.5);
				c.A = 255;
				return c;
			}
		}
		//.................................................................
		public static CLab CRGB2CLab(System.Windows.Media.Color C, CXYZ ic) { return CRGB2CLab(C.R, C.G, C.B, ic); }
		//................................................................
		public static CLab CRGB2CLab(int red, int green, int blue, CXYZ ic) {
			CXYZ xyz = CRGB2CXYZ(red, green, blue);
			return CXYZ2CLab(xyz.X, xyz.Y, xyz.Z, ic);
		}
		//..................................................................
		public static CLab CXYZ2CLab(double x, double y, double z, CXYZ icc) {
			CLab lab = new CLab();
			lab.L = 116.0 * Fxyz(y / icc.Y) - 16;
			lab.A = 500.0 * (Fxyz(x / icc.X) - Fxyz(y / icc.Y));
			lab.B = 200.0 * (Fxyz(y / icc.Y) - Fxyz(z / icc.Z));
			return lab;
		}
		//.................................................................
		public static CXYZ CRGB2CXYZ(int red, int green, int blue) {
			// normalize red, green, blue values
			double rL = (double)red / 255.0;
			double gL = (double)green / 255.0;
			double bL = (double)blue / 255.0;

			// convert to a sRGB form
			double r = Byte2R(red);
			double g = Byte2R(green);
			double b = Byte2R(blue);

			//			double r = (rL > 0.04045) ? Math.Pow((rL + 0.055) / (1 + 0.055), 2.2) : (rL / 12.92);
			//			double g = (gL > 0.04045) ? Math.Pow((gL + 0.055) / (1 + 0.055), 2.2) : (gL / 12.92);
			//			double b = (bL > 0.04045) ? Math.Pow((bL + 0.055) / (1 + 0.055), 2.2) : (bL / 12.92);
			// converts
			return new CXYZ(
				(r * 0.4124 + g * 0.3576 + b * 0.1805),
				(r * 0.2126 + g * 0.7152 + b * 0.0722),
				(r * 0.0193 + g * 0.1192 + b * 0.9505)
				);
		}
		//.................................................................
		private static double Byte2R(int v) {
			double vpc = ((double)v) / byte.MaxValue;
			double ans = (vpc > 0.04045) ? Math.Pow((vpc + 0.055) / 1.055, 2.2) : (vpc / 12.92);
			return ans;
		}
		//.................................................................
		public static CXYZ CLab2CXYZ(double l, double a, double b, CXYZ icc) {   //go back from Lab to XYZ.  icc should be either CXYZ.D50 or CXYZ.D65
			const double delta = 0.206896552, deltasq = delta * delta, rq = 0.137931034;
			double fy = (l + 16) / 116.0;
			double fx = fy + (a / 500.0);
			double fz = fy - (b / 200.0);
			return new CXYZ(
				(fx > delta) ? icc.X * (fx * fx * fx) : (fx - rq) * 3 * deltasq * icc.X,
				(fy > delta) ? icc.Y * (fy * fy * fy) : (fy - rq) * 3 * deltasq * icc.Y,
				(fz > delta) ? icc.Z * (fz * fz * fz) : (fz - rq) * 3 * deltasq * icc.Z
				);
		}
		//.................................................................
		public static CIEXYZ RGB2XYZ(int red, int green, int blue) {
			// normalize red, green, blue values
			double rL = (double)red / 255.0;
			double gL = (double)green / 255.0;
			double bL = (double)blue / 255.0;

			// convert to a sRGB form
			double r = (rL > 0.04045) ? Math.Pow((rL + 0.055) / (1 + 0.055), 2.2) : (rL / 12.92);
			double g = (gL > 0.04045) ? Math.Pow((gL + 0.055) / (1 + 0.055), 2.2) : (gL / 12.92);
			double b = (bL > 0.04045) ? Math.Pow((bL + 0.055) / (1 + 0.055), 2.2) : (bL / 12.92);
			// converts
			return new CIEXYZ(
				(r * 0.4124 + g * 0.3576 + b * 0.1805),
				(r * 0.2126 + g * 0.7152 + b * 0.0722),
				(r * 0.0193 + g * 0.1192 + b * 0.9505)
				);

		}
		//.................................................................
		public static CIELab RGB2Lab(System.Windows.Media.Color C, CIEXYZ icc) { return RGB2Lab(C.R, C.G, C.B, icc); }
		//................................................................
		public struct CIERGB {
			byte a;
			byte r;
			byte g;
			byte b;
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			public CIERGB(byte R, byte G, byte B) { a = 255; r = R; g = G; b = B; }
			public CIERGB(System.Windows.Media.Color C) { a = C.A; r = C.R; g = C.G; b = C.B; }
			public byte A { get => a; set => a = value; }
			public byte R { get => r; set => r = value; }
			public byte G { get => g; set => g = value; }
			public byte B { get => b; set => b = value; }
			public System.Windows.Media.Color Color {
				get {
					System.Windows.Media.Color c = new Color();
					c.A = 255; c.R = r; c.G = g; c.B = b;
					return c;
				}
			}
			public static double operator -(CIERGB c1, CIERGB c2) { return CSNC.ColorDistance(c1.Color, c2.Color); }
			public static bool operator ==(CIERGB O1, CIERGB O2) { return (O1.a == O2.a && O1.r == O2.r && O1.g == O2.g && O1.b == O2.b); }
			public static bool operator !=(CIERGB O1, CIERGB O2) { return (O1.a != O2.a || O1.r != O2.r || O1.g != O2.g || O1.b != O2.b); }
			public override int GetHashCode() => base.GetHashCode();
			public override bool Equals(object obj) {
				if (obj == null || GetType() != obj.GetType()) return false;
				return (this == (CIERGB)obj);
			}
			public override string ToString() {
				return '[' + r.ToString("d0") + ',' + g.ToString("d0") + ',' + b.ToString("d0") + ']';
			}

		}
		///////////////////////////////////////////////////////////////
		public struct CIELab {
			//use structure
			private double l;
			private double a;
			private double b;
			private System.Windows.Media.Color rgb;
			private System.Windows.Media.Color newrgb;

			public static bool operator ==(CIELab L1, CIELab L2) { return (L1.L == L2.L && L1.A == L2.A && L1.B == L2.B); }
			public static bool operator !=(CIELab L1, CIELab L2) { return (L1.L != L2.L || L1.A != L2.A || L1.B != L2.B); }
			public double L { get => l; set => l = value; }
			public double A { get => a; set => a = value; }
			public double B { get => b; set => b = value; }
			public System.Windows.Media.Color RGB { get => rgb; set => rgb = value; }
			public System.Windows.Media.Color NewRGB { get => newrgb; set => newrgb = value; }
			public string sRGB { get => sColor(rgb); }
			//============================================================
			public override string ToString() {
				return "L: " + l.ToString("f1") + ", a: " + a.ToString("f1") + ", b: " + b.ToString("f1");
			}
			//============================================================
			public override int GetHashCode() => base.GetHashCode();
			//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			public CIELab(double L, double A, double B) {
				l = L; a = A; b = B;
				rgb = System.Windows.Media.Colors.Black;
				newrgb = System.Windows.Media.Colors.White;
			}
			public CIELab(byte red, byte green, byte blue, CIEXYZ icc) {
				CIELab temp = RGB2Lab(red, green, blue, icc);
				l = temp.L; a = temp.A; b = temp.B;
				rgb = System.Windows.Media.Colors.Black;
				rgb.R = red; rgb.G = green; rgb.B = blue;
				newrgb = System.Windows.Media.Colors.White;
			}

			public override bool Equals(Object obj) {
				if (obj == null || GetType() != obj.GetType()) return false;
				return (this == (CIELab)obj);
			}
		}
		//////////////////////////////////////////////////////////////////
		public struct CIEXYZ {  //use structure
			public static readonly CIEXYZ D65 = new CIEXYZ(0.9505, 1.0, 1.0890);
			public static readonly CIEXYZ D50 = new CIEXYZ(0.9642, 1.0000, 0.8251);
			private double x;
			private double y;
			private double z;

			public static bool operator ==(CIEXYZ item1, CIEXYZ item2) {
				return (item1.X == item2.X && item1.Y == item2.Y && item1.Z == item2.Z);
			}

			public static bool operator !=(CIEXYZ item1, CIEXYZ item2) {
				return (item1.X != item2.X || item1.Y != item2.Y || item1.Z != item2.Z);
			}
			public double X { get => x; set { this.x = (value > 0.9505) ? 0.9505 : ((value < 0) ? 0 : value); } }
			public double Y { get => y; set { y = (value > 1.0) ? 1.0 : (value < 0.0) ? 0.0 : value; } }
			public double Z { get => z; set { z = (value > 1.089) ? 1.089 : (value < 0.0) ? 0.0 : value; } }
			public override int GetHashCode() => base.GetHashCode();
			//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
			public CIEXYZ(double X, double Y, double Z) {
				x = (X > 0.9505) ? 0.9505 : ((X < 0) ? 0 : X);
				y = (Y > 1.0) ? 1.0 : ((Y < 0) ? 0 : Y);
				z = (Z > 1.089) ? 1.089 : ((Z < 0) ? 0 : Z);
			}
			//============================================================
			public override string ToString() {
				return "(" + x.ToString("f3") + ',' + y.ToString("f3") + ',' + z.ToString("f3") + ')';
			}
			//============================================================
			public override bool Equals(Object obj) {
				if (obj == null || GetType() != obj.GetType()) return false;
				return (this == (CIEXYZ)obj);
			}
		}
		//................................................................
		public static CIELab RGB2Lab(int red, int green, int blue, CIEXYZ icc) {
			CIEXYZ xyz = RGB2XYZ(red, green, blue);
			return XYZ2Lab(xyz.X, xyz.Y, xyz.Z, icc);
		}
		/// Converts CIEXYZ to RGB structure.
		//................................................................
		public static System.Windows.Media.Color XYZ2RGB(CIEXYZ xyz) { return XYZ2RGB(xyz.X, xyz.Y, xyz.Z); }
		//................................................................
		public static System.Windows.Media.Color XYZ2RGB(double X, double Y, double Z) {
			double[] Cn = new double[3];
			Cn[0] = X * 3.2410 - Y * 1.5374 - Z * 0.4986; // red
			Cn[1] = -X * 0.9692 + Y * 1.8760 + Z * 0.0416; // green
			Cn[2] = X * 0.0556 - Y * 0.2040 + Z * 1.0570; // blue

			for (int i = 0; i < 3; i++) {
				Cn[i] = (Cn[i] <= 0.0031308) ? 12.92 * Cn[i] : (1 + 0.055) * Math.Pow(Cn[i], (1.0 / 2.4)) - 0.055;
			}
			System.Windows.Media.Color c = new System.Windows.Media.Color();
			c.R = (byte)(Cn[0] * 255.0 + 0.5);
			c.G = (byte)(Cn[1] * 255.0 + 0.5);
			c.B = (byte)(Cn[2] * 255.0 + 0.5);
			c.A = 255;
			return c;
		}
		//.............................................................
		public static Color Lab2RGB(CLab Lab, bool isD50) {
			return new Color();
		}
		//.............................................................
		public static System.Windows.Media.Color Lab2RGB(CIELab L) { return XYZ2RGB(Lab2XYZ(L)); }
		//.................................................................
		/// XYZ to L*a*b* transformation function.
		private static double Fxyz(double t) {  //internal routine
			return ((t > 0.008856) ? Math.Pow(t, (1.0 / 3.0)) : (7.787 * t + 16.0 / 116.0));
		}
		//..................................................................
		public static CIELab XYZ2Lab(double x, double y, double z, CIEXYZ icc) {
			CIELab lab = new CIELab();
			lab.L = 116.0 * Fxyz(y / icc.Y) - 16;
			lab.A = 500.0 * (Fxyz(x / icc.X) - Fxyz(y / icc.Y));
			lab.B = 200.0 * (Fxyz(y / icc.Y) - Fxyz(z / icc.Z));
			return lab;
		}
		//................................................................
		public static CIEXYZ Lab2XYZ(CIELab lab) { return Lab2XYZ(lab.L, lab.A, lab.B, CIEXYZ.D65); }
		public static CIEXYZ Lab2XYZ(double l, double a, double b, CIEXYZ icc) {   //go back from Lab to XYZ
			const double delta = 0.206896552, deltasq = delta * delta, rq = 0.137931034;
			double fy = (l + 16) / 116.0;
			double fx = fy + (a / 500.0);
			double fz = fy - (b / 200.0);
			return new CIEXYZ(
				(fx > delta) ? icc.X * (fx * fx * fx) : (fx - rq) * 3 * deltasq * icc.X,
				(fy > delta) ? icc.Y * (fy * fy * fy) : (fy - rq) * 3 * deltasq * icc.Y,
				(fz > delta) ? icc.Z * (fz * fz * fz) : (fz - rq) * 3 * deltasq * icc.Z
				);
		}
		//............................................................
		public static double ColorDistance(System.Windows.Media.Color c1, System.Windows.Media.Color c2) {
			double d = 0;
			double dR = c1.R - c2.R;
			double dG = c1.G - c2.G;
			double dB = c1.B - c2.B;
			d = dR * dR + dG * dG + dB * dB;
			return Math.Sqrt(d);
		}
		//...........................................................
		public static string sColour(Color c) => '(' + SC(c.A) + ',' + SC(c.R) + ',' + SC(c.G) + ',' + SC(c.B) + ')';
		public static string sColor(Color c) => "R: " + c.R.ToString() + ", G: " + c.G.ToString() + ", B: " + c.B.ToString();
		public static string xSWMColor(System.Windows.Media.Color c) => "#" + XC(c.A) + XC(c.R) + XC(c.G) + XC(c.B);
		public static string xSDColor(System.Drawing.Color c) => '#' + XC(c.A) + XC(c.R) + XC(c.G) + XC(c.B);
		private static string SC(byte b) => b.ToString("000");
		//............................................................
		private static string XC(byte B) {
			string xHex = "0123456789ABCDEF";
			int j, mask = 0xf, ix = (int)B;
			StringBuilder sb = new StringBuilder(3);
			j = (B >> 4) & mask;
			sb.Append(xHex[j]);
			j = B & mask;
			sb.Append(xHex[j]);
			return sb.ToString();
		}
		//.................................................................
		public static Color Byte2SWMediaColor(byte[] b, int ix) {
			Color c = new Color();
			c.B = b[ix];
			c.G = b[ix + 1];
			c.R = b[ix + 2];
			c.A = 255;
			return c;
		}
		//.................................................................
		public static System.Drawing.Color Byte2SDColor(byte[] b, int ix) {
			System.Drawing.Color sdc = System.Drawing.Color.FromArgb(b[ix+3],b[ix+2],b[ix+1],b[ix+0]);
			return sdc;
		}
		//...........................................................
		public static Color OppositeColor(System.Windows.Media.Color color) {
			System.Windows.Media.Color oc = new Color();
			oc.A = color.A;
			oc.R = (byte)(255 - color.R);
			oc.G = (byte)(255 - color.G);
			oc.B = (byte)(255 - color.B);
			return oc;
		}
		//...........................................................
		public static Color Drawing2Media(System.Drawing.Color color) {
			return Color.FromArgb(color.A, color.R, color.G, color.B);
		}
		public static System.Drawing.Color Media2Drawing(System.Windows.Media.Color color) {
			return System.Drawing.Color.FromArgb(color.A, color.R, color.G, color.B);

		}

	}
}
