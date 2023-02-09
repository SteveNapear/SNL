//////////////////////////////////////////////////////////////////////////
using System;
using System.Collections;
using System.Collections.Generic;
using System.Drawing;
using System.Globalization;
using System.IO;
using System.Management;
using System.Runtime.InteropServices;
using System.Security;
using System.Text;
using System.Windows;

namespace SNL {
	public static partial class CSNUtil {
		//...................................................................
		public static double Dist2Pts(System.Windows.Point p1, System.Windows.Point p2) {
			double d2 = (p2.X - p1.X) * (p2.X - p1.X) + (p2.Y - p1.Y) * (p2.Y - p1.Y);
			return Math.Sqrt(d2);
		}
		//...................................................................
		public static double Dist2Pts(System.Drawing.Point p1, System.Drawing.Point p2) {
			double d2 = (p2.X - p1.X) * (p2.X - p1.X) + (p2.Y - p1.Y) * (p2.Y - p1.Y);
			return Math.Sqrt(d2);
		}
		//.....................................................................
		public static System.Windows.Size GetScreenWA {
			get {
				System.Windows.Forms.Screen S = System.Windows.Forms.Screen.PrimaryScreen;
				//				Rectangle R = S.Bounds;
				//				int bpp = S.BitsPerPixel;
				Rectangle WA = System.Windows.Forms.Screen.GetWorkingArea(new System.Drawing.Point(0, 0));
				System.Windows.Size sz = new System.Windows.Size(WA.Width, WA.Height);
				return sz;
			}
		}
		//..........................................................
		/// <summary>
		/// Reads a byte array of specified length from a starting position in the file
		/// </summary>
		/// <param name="sPath" full path of file></param>
		/// <param name="pos" position in file></param>
		/// <param name="len" # of bytes></param>
		/// <returns byte[] ></returns>
		public static byte[] ReadBlock(string sPath, int pos, int len) {
			FileStream fs = null; BinaryReader br = null; byte[] B = null; string se;
			if (sPath == null || !File.Exists(sPath) || len < 1 || pos < 0) return null;

			try {
				fs = new FileStream(sPath, FileMode.Open, FileAccess.Read);
				br = new BinaryReader(fs);
				br.BaseStream.Position = pos;
				B = new byte[len];
				int n = br.Read(B, 0, len);
				while (n < len) B[n++] = 0;
				if (n < 1) throw new Exception("Read error, no bytes read.");
				br.Close();
				fs.Close();
			} catch (Exception exp) {
				se = "\nError while reading file\n:" + sPath
					+ "\n Error code: " + exp.Message
					+ "\nStack: " + exp.StackTrace;
				B = null;
			} finally {
				if (br != null) br.Dispose();
				if (fs != null) fs.Dispose();
			}
			return B;
		}
		//...............................................................
		/// <summary>
		/// 
		/// </summary>
		/// <param name="sClassName"></param>
		/// <param name="sa"></param>
		/// <param name="sOutPath"></param>
		/// <returns></returns>
		public static bool ClassGenerator(string sClassName, string[] sa, string sOutPath) {
			if (sa == null || sa.Length < 2 || sOutPath == null) return false;
			int j, ix = 0, nterms = 0, nbad, nempty = 0; string sv; bool bOK = false;
			List<string> list = new List<string>(sa.Length + 1);
			string st, se, sdir, sfns, spath;


			foreach (string s in sa) {
				sv = s;
				nbad = 0;
				if (sv.Length == 0) { nempty += 1; sv = "Empty" + nempty.ToString("d"); } else if (!char.IsLetter(sv[0])) sv = "Alpha";
				sv = Convert2LegalChar(sv);
				sv = sv.Replace(' ', '_');
				do {
					j = list.FindIndex(x => x.CompareTo(sv) == 0);
					if (j < 0) break;
					nbad += 1;
					if (j >= 0) sv += nbad.ToString("0000");
				} while (true);
				list.Add(sv);
			}
			nterms = list.Count;
			string sCN = sClassName;
			//generate output file
			sdir = System.IO.Path.GetDirectoryName(sOutPath);
			sfns = System.IO.Path.GetFileNameWithoutExtension(sOutPath) + ".cs";
			spath = sdir + '\\' + sfns;

			//generate class cs
			FileStream fs = null; StreamWriter sw = null;
			try {
				fs = new FileStream(spath, FileMode.Create, FileAccess.Write, FileShare.None);
				sw = new StreamWriter(fs);

				//write out class definition
				st = "[serializable] public class " + sCN + " {" + "//" + nterms.ToString("n0") + " variables." +
				   "\t//list of fields for class " + sCN;
				sw.WriteLine(st);
				for (ix = 0; ix < list.Count; ix++) {
					st = "\tprivate string " + list[ix].ToLower() + ';';
					sw.WriteLine(st);
				}

				//write CQBclass instantiator
				sw.WriteLine("\t//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
				st = "\t" + "public " + sCN + "(string[] sa) {";
				sw.WriteLine(st);
				//write instantiator internal routine
				for (ix = 0; ix < list.Count; ix++) {
					string s1 = list[ix].ToLower();
					string s2 = " = sa[" + ix.ToString("d") + "];";
					sw.WriteLine("\t\t" + s1 + s2);
				}
				sw.WriteLine("\t\t}");
				//write public data access routines
				sw.WriteLine("\t//========================================================================");
				for (ix = 0; ix < list.Count; ix++) {
					st = "\tpublic string s" + CSNUtil.ULower(list[ix]) + "{get => " + list[ix].ToLower() + ";}";
					sw.WriteLine(st);
				}
				//write TryParse
				sw.WriteLine("\t\t//...............................................................");
				st = "\tpublic static bool TryParse(string s, out " + sCN + " obj){";
				sw.WriteLine(st);
				sw.WriteLine("\t\t\tint ix;");
				st = "\t\t\tobj = null;"; sw.WriteLine(st);
				st = "\t\t\tif(s == null) return false;"; sw.WriteLine(st);
				st = "\t\t\tstring[]sa = CDASlib.CSVSplit(s, '\t');"; sw.WriteLine(st);
				st = "\t\t\tif(sa.Length != " + nterms.ToString() + ") return false;"; sw.WriteLine(st);
				string sname = sCN;

				st = "\t\t\tobj = new " + sCN + "(sa);";
				sw.WriteLine(st);

				sw.WriteLine("\t\t\treturn (obj!=null); \n\t}");


				sw.WriteLine("}"); //class end

				//do ClassList
				string sCL = sCN + "List";
				string sCLHeader = sCL + ": List<" + sCN + ">";
				sw.WriteLine(@"//////////////////////////////////////////////////////////////////////////");
				sw.WriteLine("[serializeable] public class " + sCLHeader + " {");

				sw.WriteLine("\tprivate string spath = \"\";");
				sw.WriteLine("\t//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
				sw.WriteLine("\tpublic " + sCL + "(int capacity) : base(capacity) {}");
				sw.WriteLine("\t//=====================================================================");
				sw.WriteLine("\tpublic string sPath { get => spath; }");
				sw.WriteLine("\t//......................................................................");
				st = "\tpublic static " + sCL + " ReadCSV(string sP) {";
				sw.WriteLine(st);
				sw.WriteLine("\t\tif (!System.IO.File.exists(sP)) return null;");

				sw.WriteLine("\t\tSystem.IO.FileStream fs = null; System.IO.StreamReader sr = null;");
				st = "\t\t" + sCL + " list = null;";
				sw.WriteLine(st);
				sw.WriteLine("\t\t string sl, se; CSKU sku = null;");

				sw.WriteLine("\t\ttry {");
				sw.WriteLine("\t\t\tfs = new System.IO.FileStream(sP, System.IO.FileMode.Open, System.IO.FileAccess.Read, System.IO.FileShare.Read);");
				sw.WriteLine("\t\t\tif (fs == null) return null;");
				sw.WriteLine("\t\t\tsr = new System.IO.StreamReader(fs);");
				sw.WriteLine("\t\t\tsl = sr.ReadLine();         //ignore header line");
				sw.WriteLine("\t\t\tlist = new CSKUList(1000);");
				sw.WriteLine("\t\t\twhile ((sl = sr.ReadLine()) != null) {");
				sw.WriteLine("\t\t\t\tif (!CSKU.TryParse(sl, out sku)) continue;");
				sw.WriteLine("\t\t\t\tlist.Add(sku);");
				sw.WriteLine("\t\t\t}");
				sw.WriteLine("\t\t\tsrt.close();");
				sw.WriteLine("\t\t\tfs.Close();");

				sw.WriteLine("\t\t} catch (Exception exp) {");
				sw.WriteLine("\t\t\tse = \"Error while reading CSV file:\\n\" + sP");
				sw.WriteLine("\t\t\t+ \"\\nError: \" + exp.Message");
				sw.WriteLine("\t\t\t+ \"\\nStack: \" + exp.StackTrace;");
				sw.WriteLine("\t\t\tMessageBox.show(se);");
				sw.WriteLine("\t\t\tlist = null;");
				sw.WriteLine("\t\t} finally {");
				sw.WriteLine("\t\t\tif (sr != null) sr.Close();");
				sw.WriteLine("\t\t\tif (fs != null) fs.Close();");
				sw.WriteLine("\t\t}");
				sw.WriteLine("\t\tlist.TrimExcess();");
				sw.WriteLine("\t\treturn list;");
				sw.WriteLine("\t}");

				sw.WriteLine("}");
				sw.WriteLine();
				sw.WriteLine();

				sw.Flush();
				sw.Close();
				fs.Close();
				bOK = true;
			} catch (Exception exp) {
				se = "Error while writing file:\n" + spath +
					"\nError code: " + exp.Message +
					"\nStack trace: " + exp.StackTrace;
				bOK = false;
			} finally {
				if (sw != null) sw.Dispose();
				if (fs != null) fs.Dispose();
			}
			return bOK;
		}

		//................................................................
		/// <summary>
		/// For a byte[] array trys to find the position within the array containing string s
		/// </summary>
		/// <param name="A" byte[] must be at least 1 char long></param>
		/// <param name="s" string will be converted to a byte arrary></param>
		/// <returns></returns>
		public static int Findindex(byte[] A, string s) {
			if (s == null || s.Length == 0) return -1;
			char[] C = s.ToCharArray();
			byte[] B = new byte[C.Length];
			for (int ix = 0; ix < C.Length; ix++) B[ix] = (byte)C[ix];
			return Findindex(A, B);
		}
		//................................................................
		/// <summary>
		/// For a byte[] array trys to find the position within the array containing byte array B
		/// </summary>
		/// <param name="A"></param>
		/// <param name="B"></param>
		/// <returns></returns>
		public static int Findindex(byte[] A, byte[] B) {
			int ix, j;

			if (B.Length > A.Length) return -1;
			ix = 0; j = 0;
			do {
				j = 1;
				while (A[ix++] != B[0] && ix < A.Length) ;
				if (ix >= A.Length) return -1;
				//ASSERT(A[ix] == B[0]);
				while (A[ix++] == B[j++] && j < B.Length) ;
				if (j == B.Length) return ix - j;

			} while (ix < A.Length && j < B.Length);
			return -1;
		}
		//................................................................
		/// <summary>
		/// Converts string to letter or digit values.  Non => '_'
		/// </summary>
		/// <param name="st" string to be converted></param>
		/// <returns></returns>
		private static string Convert2LegalChar(string st) {
			StringBuilder sb = new StringBuilder(st.Length);
			string s = "";
			for (int ix = 0; ix < st.Length; ix++) sb.Append(char.IsLetterOrDigit(st[ix]) ? st[ix] : '_');
			s = sb.ToString();
			return s;
		}
		//................................................................
		private const char degsign = (char)0xba;  //degree sign
												  //general utility routines
												  //.....................................................................
												  //	1000 kB kilobyte 
												  //10002 MB megabyte 
												  //10003 GB gigabyte 
												  //10004 TB terabyte 
												  //10005 PB petabyte 
												  //10006 EB exabyte 
												  //10007 ZB zettabyte 
												  //10008 YB yottabyte 
		private static char[] AO = { 'B', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y' };
		//.....................................................................
		public static string sO(long nbytes)    //returns compressed Order of Magnitude string
		{
			int i = 0; ulong n = (ulong)nbytes, div = 1000;
			while (n >= div) { i += 1; n /= div; }
			return (i < AO.Length) ? n.ToString("n0") + AO[i] : nbytes.ToString("n0") + 'B';
		}
		//.....................................................................
		public static bool cmp2s(string s, string s2)   //compares beginning of s1 to s2
		{
			int ix = 0;
			if (s2.Length > s.Length) return false;
			while (ix < s2.Length) { if (s[ix] > s2[ix]) return false; else if (s[ix] < s2[ix]) return false; ix += 1; }
			return true;
		}
		//...................................................................
		private static char u2c(ulong ul) { return (char)(((ulong)'0') + (ul % 10)); }
		//.....................................................................
		private static string u2s(ulong u, string st) {
			char c = u2c(u);
			u /= 10;
			if (u > 0) st += u2s(u, st);
			st += c;
			return st;
		}
		//................................................................
		public static double Removedollarsign(string s) /*remove '$', return value*/{
			if (s is null || s.Length == 0) return 0;
			string ss = s.Trim();
			if (ss[0] == '$') ss = ss.Substring(1);
			if (!double.TryParse(ss, out double v)) return 0;
			return v;
		}

		//................................................................
		public static string sExponent(char X, int expenent)   //returns x superscript n
		{
			char[] SuperDigits = { '\u00ba', '\u00b9', '\u00b2', '\u00b3', '\u2074', '\u2075', '\u2076', '\u2077', '\u2078', '\u2079' };
			char e, minus = '\u02c9';
			int[] A = { 0, 0, 0, 0 };
			if (expenent > 9999) return "?????";  //maximum of four digits for exponent
			if (expenent == 0) return " ";      //no exponent string when exponent is 1
			int i, ix, v = expenent < 0 ? -expenent : expenent;      //if n neg, use positive and remember sign
			string ss = " "; ss += X;
			ix = 0; while (v > 0) { A[ix++] = v % 10; v /= 10; }
			if (expenent < 0) ss += minus;
			do {
				i = A[--ix];
				e = SuperDigits[i];
				ss += e;
			} while (ix > 0);
			return ss;
		}
		//................................................................
		public static string sRect(System.Windows.Rect R) {
			string f = "f1";
			return '[' + R.X.ToString(f) + ',' + R.Y.ToString(f) + ',' + R.Width.ToString(f) + ',' + R.Height.ToString(f) + ']';
		}
		//...............................................................
		public static string sRect(System.Drawing.Rectangle R) => '[' + R.X.ToString("d0") + ',' + R.Y.ToString("d0") + ',' + R.Width.ToString("d0") + ',' + R.Height.ToString("d0") + ']';
		public static string sRect(System.Drawing.RectangleF R) => '[' + R.X.ToString("f1") + ',' + R.Y.ToString("f1") + ',' + R.Width.ToString("f1") + ',' + R.Height.ToString("f1") + ']';
		//................................................................
		public static string sSize(System.Windows.Size sz) => '(' + sz.Width.ToString("f2") + ',' + sz.Height.ToString("f2") + ')';
		//................................................................
		public static string sSize(System.Drawing.Size sz) => '(' + sz.Width.ToString("d0") + ',' + sz.Height.ToString("d0") + ')';
		//...............................................................
#pragma warning disable IDE1006 // Naming Styles
		public static string sPoint(System.Windows.Point p) => '(' + p.X.ToString("000.00") + ',' + p.Y.ToString("000.00") + ')';
#pragma warning restore IDE1006 // Naming Styles
		//................................................................
		public static string sColor(System.Windows.Media.Color c) => '(' + c.A.ToString("000") + ',' + c.R.ToString("000") + ',' + c.G.ToString("000") + ',' + c.B.ToString("000") + ')';
		//................................................................
		public static string sPoint(System.Drawing.Point p) => '(' + p.X.ToString("000") + ',' + p.Y.ToString("000") + ')';
		public static char sDegree => '\u00ba';
		//................................................................
		public static string sI2SRgtAdj(int value, int fieldlength)     //int to string, right adjusted w fixed field length
		{
			bool bSign = (value < 0); int ix, v = (bSign) ? -value : value;
			if (fieldlength < 1) return "???";
			if (fieldlength > 255) fieldlength = 255;
			char[] ln = new char[fieldlength];
			for (ix = 0; ix < ln.Length - 1; ix++) ln[ix] = ' '; ln[ix] = '\0';
			do { ln[--fieldlength] = (char)('0' + v % 10); v /= 10; } while ((fieldlength > 1) && (v > 0));
			if (bSign) ln[--fieldlength] = '-';
			StringBuilder sb = new StringBuilder(ln.Length + 1);
			for (ix = 0; ix < ln.Length; ix++) sb.Append(ln[ix]);
			return sb.ToString();
		}
		//...................................................................
		public static string ULower(string s) //Set 1st char to upper, remaining to lower case
		{
			if ((s == null) || (s.Length == 0)) return " ";
			if (s.Length == 1) return s.ToUpper();
			return Char.ToUpper(s[0]) + s.Substring(1).ToLower();
		}
		//...................................................................
		public static string f2c(float v)   //float to currency $99.99
		{
			CultureInfo ci = new CultureInfo("en-us");
			return v.ToString("C", ci);
		}
		//....................................................................
		public static bool Currency2float(string s, out float v)    /*converts -$9.99 to float*/	{
			v = 0.0f;
			if (s is null || s.Length == 0) return false;
			string sm, st = s.Trim();
			if (st.Length == 0) return false;

			int i = 0, len = s.Length; bool b;
			if (i >= len) return false;

			b = (st[i] == '-'); if (b) i += 1;
			if (i >= len) return false;
			if (st[i] == '$') i += 1;

			sm = st.Substring(i);
			bool bOK = float.TryParse(sm, out v);
			if (b) v = -v;
			return bOK;
		}
		//.......................................................................
		public static bool Currency2Double(string s, out double v)/*converts -$8.88 to double*/{
			v = 0.0;
			if (s == null || s.Length == 0) return false;
			string sm, st = s.Trim();
			if (st.Length == 0) return false;

			int i = 0, len = s.Length; bool b;
			if (i >= len) return false;

			b = (st[i] == '-'); if (b) i += 1;
			if (i >= len) return false;
			if (st[i] == '$') i += 1;

			sm = st.Substring(i);
			bool bOK = double.TryParse(sm, out v);
			if (b) v = -v;
			return bOK;

		}
		//..................................................................
		// CSVHeader2IA returns an integer array containing indecies of the sEnum string array.
		//bUpper, if true, will make compare Upper case strings in both the header and the sEnum list
		public static int[] CSVHeader2IA(string sHeader, string[] sEnum, bool bUpper) //
		{
			if (sHeader == null || sEnum == null || sEnum.Length < 2) return null;
			string[] sa = CSVSplit(sHeader, ',');
			if (sa.Length < 2) return null;

			int i, ix, mx = Math.Max(sa.Length, sEnum.Length);

			int[] IX = new int[mx]; string sn, st;
			for (ix = 0; ix < IX.Length; ix++) {
				IX[ix] = -1;
				if (ix >= sEnum.Length) continue;  //ignore remaining entries
				sn = (bUpper ? sEnum[ix].ToUpper() : sEnum[ix]);
				for (i = 0; i < sa.Length; i++) {
					st = (bUpper ? sa[i].ToUpper() : sa[i]);
					if (st == sn) { IX[ix] = i; break; }
				}
			}
			return IX;
		}
		//.....................................................................
		public static string WrapinQuotes(string s) { return '"' + s + '"'; }
		//.....................................................................
		public static string Encap(string ss, char c)//wraps quotes if special char present
		{
			string st = "";
			if (ss == null) return st;
			if (ss.IndexOf(c) < 0) return ss;
			st = '\"' + ss + '\"';
			return st;
		}
		//.....................................................................
		public static string Q(string s) => (s == null || s.Length == 0) ? string.Empty : '\"' + s + '\"';
		//.....................................................................
		public static string ExpandTabs(string s, int tabsize) {
			if (s == null || s.Length == 0) return null;
			StringBuilder sb = new StringBuilder(s.Length * 2);

			int ix = 0; char c;
			for (ix = 0; ix < tabsize; ix++) sb.Append(' ');
			string sblank = sb.ToString();
			sb.Clear();
			ix = 0;
			while (ix < s.Length) {
				c = s[ix++];
				if (c == '\t') sb.Append(sblank);
				else sb.Append(c);
			}
			return sb.ToString();
		}
		//.....................................................................
		public static string Decap(string ss)   /*removes leading and/or trailing quotes*/{
			if (ss == null) return null;
			string delim = "\"";
			return ss.Trim(delim.ToCharArray());
		}
		//.....................................................................
		public static string Pack(string ss)        //removes multiple blanks
		{
			int i, j, len = ss.Length; char[] sa = new char[len];
			try {
				i = 0; j = 0;
				do {
					while ((i < len) && (ss[i] == ' ')) i += 1;
					while ((i < len) && (ss[i] != ' ')) sa[j++] = ss[i++]; sa[j++] = ' ';
				} while (i < len);
				string s = new string(sa, 0, j);
				return s;
			} catch { return null; }
		}
		//.....................................................................
		public static DateTime GetYYMMDD(string ss) {
			ss = ss.Trim(); int yr, mn, dy;

			string[] sa = ss.Split('.');
			yr = 2000 + 10 * (sa[0][0] - '0') + (sa[0][1] - '0');
			mn = 10 * (sa[1][0] - '0') + (sa[1][1] - '0');
			dy = 10 * (sa[2][0] - '0') + (sa[2][1] - '0');
			DateTime dt = new DateTime(yr, mn, dy);
			return dt;
		}
		//................................................................
		public static int Getbitstream(byte[] A, int start, int len)//builds Little Endian stream
		{
			int v = 0, d, r;
			d = start / 8; r = start % 8;
			if (d >= A.Length) return v; else v = A[d] << 24;
			if (d + 1 >= A.Length) return v; else v += A[d + 1] << 16;
			if (d + 2 >= A.Length) return v; else v += A[d + 2] << 8;
			if (d + 3 >= A.Length) return v; else v += A[d + 3];
			v <<= r;
			v >>= (32 - len);
			v &= GetMask(len);
			return v;
		}
		//................................................................
		public static int GetMask(int len) {
			int m = 0;
			while (len-- > 0) { m <<= 1; m += 1; }
			return m;
		}
		//.................................................................
		public static short Byte2ShortBigEndian(byte[] B, int offset) {
			short u = 0; int v1, v2, v3;
			if (offset + 2 > B.Length) return u;
			v1 = B[offset + 1];
			v2 = (B[offset] << 8);
			v3 = v2 + v1;
			u = (short)v3;
			return u;
		}
		//.................................................................
		public static short Byte2ShortLittleEndian(byte[] B, int offset) {
			short u = 0; int v1, v2, v3;
			if (offset + 2 > B.Length) return u;
			v1 = B[offset];
			v2 = (B[offset + 1] << 8);
			v3 = v2 + v1;
			u = (short)v3;
			return u;
		}
		//...................................................................
		public static short Byte2ShortEndian(byte[] B, int offset, bool bBigEndian) {
			if (bBigEndian) return Byte2ShortBigEndian(B, offset);
			return Byte2ShortLittleEndian(B, offset);
		}
		//...................................................................
		public static int Byte2IntLittleEndian(byte[] B, int offset) {
			int v = 0;
			if (offset + 4 > B.Length) return 0;
			v = B[offset];
			v += (B[offset + 1] << 8);
			v += (B[offset + 2] << 16);
			v += (B[offset + 3] << 24);
			return v;
		}
		//....................................................................
		public static System.Drawing.Color Byte2SDColor(byte[] B,int offset) /*Little Endian order B[0]:B, B[1]:G ...*/{
			System.Drawing.Color c = Color.FromArgb(B[offset+3],B[offset+2],B[offset+1],B[offset+0]);
			return c;
		}
		//....................................................................
		public static System.Drawing.Color Byte2SDColor(byte[] B,int offset,byte Alpha) /*Little Endian order B[0]:B, B[1]:G ...*/{
			System.Drawing.Color c = Color.FromArgb(Alpha,B[offset+2],B[offset+1],B[offset+0]);
			return c;
		}
		//....................................................................
		public static int Byte2IntBigEndian(byte[] B, int offset) {
			int v = 0;
			if (offset + 4 > B.Length) return 0;
			v = B[offset + 3];
			v += (B[offset + 2] << 8);
			v += (B[offset + 1] << 16);
			v += (B[offset] << 24);
			return v;
		}
		//.....................................................................
		public static int Byte2IntEndian(byte[] B, int offset, bool bBigEndian) {
			if (bBigEndian) return Byte2IntBigEndian(B, offset);
			return Byte2IntLittleEndian(B, offset);
		}
		//.....................................................................
		public static long Byte2LongEndian(Byte[] B, int offset, bool bBigEndian) {
			if (B == null || B.Length < 8 + offset) return 0L;
			long L = 0; int ix = 0;
			if (bBigEndian) for (ix = 0; ix < 8; ix++) { L <<= 8; L += B[ix + offset]; } else for (ix = 7; ix >= 0; ix--) { L <<= 8; L += B[ix + offset]; }
			return L;
		}
		//................................................................
		public static int SwitchEndian(int v) {
			int i, w = 0; int[] V = new int[4];
			for (i = 0; i < V.Length; i++) { V[i] = (v & 0xff); v >>= 8; }
			for (i = 0; i < V.Length; i++) { w <<= 8; w |= V[i]; }
			return w;
		}
		//................................................................
		public static uint SwitchEndian(uint v) {
			uint i, w = 0; uint[] V = new uint[4];
			for (i = 0; i < V.Length; i++) { V[i] = (v & 0xff); v >>= 8; }
			for (i = 0; i < V.Length; i++) { w <<= 8; w |= V[i]; }
			return w;
		}
		//................................................................
		public static short SwitchEndian(short v) {
			byte b0, b1;
			b0 = (byte)(v >> 8);
			b1 = (byte)(v & 0xff);
			int iv = (b1 << 8) + b0;
			return (short)iv;
		}
		//................................................................
		public static ushort SwitchEndian(ushort v) {
			byte b0, b1;
			b0 = (byte)(v >> 8);
			b1 = (byte)(v & 0xff);
			int iv = (b1 << 8) + b0;
			return (ushort)iv;
		}
		//................................................................
		public static void SetInt2BytesLE(int v, byte[] B, int offset) {
			B[offset + 0] = (byte)(v & 0xff);
			B[offset + 1] = (byte)((v >> 8) & 0xff);
			B[offset + 2] = (byte)((v >> 16) & 0xff);
			B[offset + 3] = (byte)((v >> 24) & 0xff);

		}
		//..............................................................
		public static void SetInt2BytesBE(int v, byte[] B, int offset) {
			B[offset + 3] = (byte)(v & 0xff);
			B[offset + 2] = (byte)((v >> 8) & 0xff);
			B[offset + 1] = (byte)((v >> 16) & 0xff);
			B[offset + 0] = (byte)((v >> 24) & 0xff);

		}
		//................................................................
		public static void SetShort2BytesLE(int v, byte[] B, int offset) {
			B[offset + 0] = (byte)(v & 0xff);
			B[offset + 1] = (byte)((v >> 8) & 0xff);
		}
		//................................................................
		public static void SetShort2BytesBE(int v, byte[] B, int offset) {
			B[offset + 1] = (byte)(v & 0xff);
			B[offset + 0] = (byte)((v >> 8) & 0xff);
		}
		//................................................................
		public static void Int2ByteLE(int v, byte[] B, int offset) {
			B[offset + 0] = (byte)(v & 0xff);
			B[offset + 1] = (byte)((v >> 8) & 0xff);
			B[offset + 2] = (byte)((v >> 16) & 0xff);
			B[offset + 3] = (byte)((v >> 24) & 0xff);
		}
		//................................................................
		public static void Int2ByteBE(int v, byte[] B, int offset) {
			B[offset + 3] = (byte)(v & 0xff);
			B[offset + 2] = (byte)((v >> 8) & 0xff);
			B[offset + 1] = (byte)((v >> 16) & 0xff);
			B[offset + 0] = (byte)((v >> 24) & 0xff);
		}
		//................................................................
		public static char Byte2Unicode(byte[] B, int offset)/*Assume BigEndian*/ {
			if (offset + 1 > B.Length) return '\0';
			int ians, w0 = B[offset + 1], w1 = B[offset];
			ians = (w1 << 8) + w0;
			return (char)ians;
		}
		//................................................................
		public static char Byte2Unicode(byte[] B, int offset, bool bBigEndian) {
			if (offset + 1 >= B.Length) return '\0';
			int ians, w0 = B[offset + 1], w1 = B[offset];
			if (bBigEndian) ians = (w1 << 8) + w0;
			else ians = (w0 << 8) + w1;
			return (char)ians;
		}
		//===============================================================
		public static string Byte2UniString(byte[] B, int off, int len) {
			int blen = B.Length; char c;
			if (off + len >= blen) return string.Empty;
			StringBuilder sb = new StringBuilder(len + 1);
			for (int ix = 0; ix < len; ix += 2) { c = Byte2Unicode(B, off + ix, true); sb.Append(c); }
			string st = sb.ToString();
			return st;
		}

		//................................................................
		public static string ByteArray2Hex(byte[] B, int offset, int length) {
			StringBuilder sb = new StringBuilder(B.Length);
			int ix = 0; char[] ca = null;
			sb.Append("x");

			while (ix < length && ix < B.Length) {
				ca = Byte2Hex(B[ix++]);
				sb.Append(ca[0]);
				sb.Append(ca[1]);
			}
			return sb.ToString();
		}
		//.............................................................
		public static char[] Hex = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E', 'F' };
		public static char[] Byte2Hex(byte B) {
			char[] ca = new char[2];
			int v = B, mask = 0xf;

			ca[0] = (char)(v & mask);
			ca[1] = (char)((v >> 4) & mask);
			return ca;
		}
		//................................................................
		public static string NoControlChars(string st) {
			StringBuilder sb = new StringBuilder(st.Length);
			for (int ix = 0; ix < st.Length; ix++) {
				if (char.IsControl(st[ix])) continue;
				sb.Append(st[ix]);
			}
			return sb.ToString();
		}
		//................................................................
		public static string mmddyy(DateTime dt)    //returns mm/dd/yy
		{
			char[] ca = new char[8];
			int i = 0, v = dt.Month;
			ca[i++] = (char)('0' + v / 10);
			ca[i++] = (char)('0' + v % 10);
			ca[i++] = '/';
			v = dt.Day;
			ca[i++] = (char)('0' + v / 10);
			ca[i++] = (char)('0' + v % 10);
			ca[i++] = '/';
			v = (dt.Year % 100);
			ca[i++] = (char)('0' + v / 10);
			ca[i++] = (char)('0' + v % 10);
			return new string(ca);
		}
		//....................................................................
		public static string yymmdd(DateTime dt) //returns yy.mm.dd
		{
			string st;

			if (dt.Year < 1900) return "00.00.00";
			st = dt.Year.ToString("D2"); st += '.';
			st = st.Remove(0, 2);
			st += dt.Month.ToString("D2"); st += '.';
			st += dt.Day.ToString("D2");
			return st;
		}
		//.....................................................................
		public static string yyyymmdd(DateTime dt) //returns yyyy.mm.dd
		{
			int i, v, year = dt.Year, month = dt.Month, day = dt.Day;
			char[] ar = new char[10];
			for (i = 0; i < 4; i++) { v = year % 10; ar[3 - i] = (char)(v + '0'); year /= 10; }
			ar[4] = '.';
			ar[6] = (char)(month % 10 + '0');
			ar[5] = (char)(month / 10 + '0');
			ar[7] = '.';
			ar[8] = (char)(day / 10 + '0');
			ar[9] = (char)(day % 10 + '0');
			StringBuilder sb = new StringBuilder();
			sb.Append(ar);
			return sb.ToString();
		}
		//....................................................................
		public static char[] I2CA(int v)    //return 12 char right justified in char array
		{
			int i = 0; bool b = (v < 0) ? true : false;
			if (b) v = -v;
			char[] ca = new char[12];
			for (i = 0; i < ca.Length - 1; i++) ca[i] = ' '; ca[i] = '\0';
			do {
				ca[--i] = (char)('0' + v % 10);
				v /= 10;
			} while (v > 0);
			ca[--i] = (b ? '-' : ' ');
			return ca;
		}
		//....................................................................
		private static char[] i2ca2(int v)  //converts v to two char array
		{
			char[] A = new char[2];
			A[0] = (char)('0' + v % 10); A[1] = (char)('0' + v / 10);
			return A;
		}
		//....................................................................
		public static string yymmddhhmmss(DateTime dt)//returns "yy.mm.dd hh:mm:ss"
		{
			StringBuilder sb = new StringBuilder(30);
			sb.Append(yymmdd(dt));
			sb.Append(' ');
			sb.Append(hhmmss(dt));
			return sb.ToString();
		}
		//.....................................................................
		public static string hhmmss(DateTime dt)    //returns "hh:mm:ss" for given datetime object
		{
			StringBuilder sb = new StringBuilder(30);
			char[] A;
			A = i2ca2(dt.Hour); sb.Append(A[1]); sb.Append(A[0]); sb.Append(':');
			A = i2ca2(dt.Minute); sb.Append(A[1]); sb.Append(A[0]); sb.Append(':');
			A = i2ca2(dt.Second); sb.Append(A[1]); sb.Append(A[0]);
			return sb.ToString();
		}
		//.....................................................................
		public static string StripQuotes(string s)//strips quotes if string is in a quote
		{
			s = s.Trim();
			if ((s.Length < 2) || (s[0] != '"')) return s;
			return s.Substring(1, s.Length - 2);
		}
		//.....................................................................
		//CSVSplit splits phrases into substrings ignoring delim char inside quotes 
		public static string[] CSVSplit(string s, Char delim)   // split for csv files 
		{
			bool bQuote = false; int i = 0, j = 0, ndelim = 1; string st; Char c;

			Char[] ca = s.ToCharArray();

			if (ca.Length == 0) return null;
			for (i = 0; i < s.Length; i++) if (ca[i] == delim) ndelim += 1;
			string[] sa = new string[ndelim];
			i = 0; do {
				st = "";
				do {
					c = ca[i++];
					if (c == '"') bQuote = !bQuote;
					else if (!bQuote && (c == delim)) break;
					else st += c;
				} while (i < s.Length);
				sa[j++] = st;
			} while (i < ca.Length);
			if (j == sa.Length) return sa;
			string[] sb = new string[j];
			for (i = 0; i < j; i++) sb[i] = sa[i];
			return sb;
		}
		//...................................................................
		public static void DrawChart(Graphics g,    //draws a chart with (optional) scale on left side
		Rectangle r,        //full rectangle for chart, including scale area
		float wdScale,      //width of scale area. If 0, then no scale is drawn
		int Ylow,           //low Y value
		int Yhigh,          //high Y value
		int Yinc,           //Y increment
		int nVertical,      //# of vertical lines (r.Width/nVertical is x increment
		System.Drawing.Color graphcolor,   //color of chart lines
		System.Drawing.Color background,   //color of background
		bool bSigned        //if set, 0 is in middle, else  bottom
		) {
			float em, x, y, yy, wd, ht, dy, dhy, htinc, deltaV; int ix, iy;
			bool bScale = (wdScale == 0.0f) ? false : true;

			em = bScale ? wdScale / 5.0f : 10.0f;
			Font lfont = new Font("Courier New", em, System.Drawing.FontStyle.Bold, GraphicsUnit.Pixel);
			SolidBrush br = new SolidBrush(graphcolor);
			Pen pen = new Pen(br);
			StringFormat sfr = new StringFormat();
			sfr.Alignment = StringAlignment.Far;
			sfr.LineAlignment = StringAlignment.Center;

			deltaV = (Yinc == 0) ? Yhigh - Ylow : (Yhigh - Ylow) / Yinc;

			g.Clear(background);
			ht = (bSigned ? r.Height / 2.0f : r.Height);
			htinc = ht / deltaV; dhy = htinc / 10;
			//draw horizontal lines
			r.X += (int)wdScale; r.Width -= (int)wdScale;   //adjust for scale ara

			y = ht; dy = 0;
			ix = Ylow;
			while ((y + dy) <= r.Height) {
				if (bScale) g.DrawString(ix.ToString("D"), lfont, br, r.X, y - dy, sfr);
				pen.Width = 2.0f; g.DrawLine(pen, r.Left, y - dy, r.Right, y - dy);
				pen.Width = 0.1f;
				for (iy = 1; iy < 10; iy++) {
					yy = y - dy - iy * dhy;
					g.DrawLine(pen, r.Left, yy, r.Right, yy);
				}
				//do mirror 
				if (bSigned && ((y + dy) != ht)) {
					if (bScale) g.DrawString((-ix).ToString("D"), lfont, br, r.X, y + dy, sfr);
					pen.Width = 2.0f; g.DrawLine(pen, r.Left, y + dy, r.Right, y + dy);
					pen.Width = 0.1f;
					for (iy = 1; iy < 10; iy++) {
						yy = y + dy - iy * dhy;
						g.DrawLine(pen, r.Left, yy, r.Right, yy);
					}
				}
				ix += Yinc;
				dy += htinc;
			}
			//draw vertical lines
			pen.Width = 0.1f;
			wd = r.Width / (1 + nVertical);
			x = r.X; while (x <= r.Right) { g.DrawLine(pen, x, r.Top, x, r.Bottom); x += wd; }

			sfr.Dispose();
			pen.Dispose();
			br.Dispose();
			lfont.Dispose();

		}
		//...................................................................
		private static float[] diy = { 0.0f, 0.1375f, 0.2630f, 0.3785f, 0.4854f, 0.5850f, 0.6781f, 0.7655f, 0.8480f, 0.9260f };
		//...................................................................
		public static void DrawLogChart(Graphics g, //draws a chart with (optional) scale on left side
			Rectangle r,        //full rectangle for chart, including scale area
			float wdScale,      //width of scale area. If 0, then no scale is drawn
			int low,            //low value of log
			int high,           //high log value
			int nVertical,      //# of vertical lines (r.Width/nVertical is x increment
			System.Drawing.Color graphcolor,   //color of chart lines
			System.Drawing.Color background,   //color of background
			bool bSigned        //if set, 0 is in middle, else  bottom
		) {
			float em, x, y, yy, wd, ht, dy, dhy, deltalog = high - low; int ix, iy;
			bool bScale = (wdScale == 0) ? false : true;

			em = bScale ? wdScale / 5.0f : 10.0f;
			Font lfont = new Font("Courier New", em, System.Drawing.FontStyle.Bold, GraphicsUnit.Pixel);
			SolidBrush br = new SolidBrush(graphcolor);
			Pen pen = new Pen(br);
			StringFormat sfr = new StringFormat();
			sfr.Alignment = StringAlignment.Far;
			sfr.LineAlignment = StringAlignment.Center;

			g.Clear(background);
			dy = r.Height / deltalog;

			ht = (bSigned ? r.Height / 2.0f : r.Height);
			//draw lines
			r.X += (int)wdScale; r.Width -= (int)wdScale;   //adjust for scale ara

			y = r.Y + r.Height;                 //start at bottom of chart
			dhy = dy;
			ix = low;
			while (y >= r.Y) {
				if (bScale) g.DrawString(ix++.ToString("D"), lfont, br, r.X, y, sfr);
				pen.Width = 2.0f; g.DrawLine(pen, r.Left, y, r.Right, y);
				pen.Width = 0.1f;
				for (iy = 1; iy < 10; iy++) {   //draw log sublines
					yy = y - dhy * diy[iy];
					g.DrawLine(pen, r.Left, yy, r.Right, yy);
				}
				y -= dy;
			}

			//draw vertical lines
			pen.Width = 0.1f;
			wd = r.Width / (1 + nVertical);
			x = r.X; while (x <= r.Right) { g.DrawLine(pen, x, r.Top, x, r.Bottom); x += wd; }

			sfr.Dispose();
			pen.Dispose();
			br.Dispose();
			lfont.Dispose();

		}
		//................................................................
		public static double[] LeastSquares(double[] X, double[] Y, out double Slope, out double YIntercept) {
			Slope = 0.0; YIntercept = 0.0;
			if (X.Length != Y.Length || X.Length <= 1) return null;
			double avgx = 0.0, avgy = 0.0, n = X.Length, nm1 = (1.0 / (n - 1.0)); int ix;
			double S12, Sxy, sx = 0.0, sx2 = 0.0, sxy = 0.0, sy = 0.0;

			avgx = Avg(X); avgy = Avg(Y); sx = Sum(X); sy = Sum(Y);
			for (ix = 0; ix < n; ix++) { sx2 += X[ix] * X[ix]; sxy += X[ix] * Y[ix]; }

			S12 = (sx2 - sx * sx / n);
			Sxy = (sxy - sx * sy / n);
			Slope = Sxy / S12;
			YIntercept = avgy / (Slope * avgx);

			double[] ans = new double[X.Length];
			for (ix = 0; ix < ans.Length; ix++) ans[ix] = (X[ix] - avgx) * Slope + avgy;
			return ans;
		}
		//....................................................................\
		public static double Avg(double[] A) {
			if (A == null || A.Length == 0) return 0.0;
			if (A.Length == 1) return A[0];
			double d = Sum(A);
			return d / A.Length;
		}
		public static double Sum(double[] A) {
			double sum = 0.0;
			if (A == null || A.Length == 0) return 0.0;
			for (int ix = 0; ix < A.Length; ix++) sum += A[ix];
			return sum;
		}
		//....................................................................
		public static double Gaussian(double x, double μ, double σ) {
			double ee = (x - μ);
			double ee2 = ee * ee;
			double exp = ee2 / (2 * σ * σ);
			double v = Math.Exp(-exp);
			double den = Math.Sqrt(2 * Math.PI * σ * σ);
			double a = 1.0 / den;
			return a * v;
		}
		//.....................................................................
		public static double Exp(double x, double μ, double σ)  //adjusted so max == 1.0
		{
			return Math.Sqrt(2 * Math.PI * σ * σ) * (Gaussian(x, μ, σ));
		}
		//....................................................................
		public static double Getμ(double[] A) {
			if (A == null || A.Length < 1) return 0.0;
			double sum = 0.0;
			for (int ix = 0; ix < A.Length; ix++) sum += A[ix];
			return sum / A.Length;
		}
		//....................................................................
		public static double Getσ(double μ, double[] A) {
			if (A == null || A.Length < 1) return 0;
			double v, variance = 0.0;
			for (int ix = 0; ix < A.Length; ix++) { v = A[ix] - μ; variance += (v * v); }
			return Math.Sqrt(variance / A.Length);
		}
		//....................................................................
		public static PointF MapGaussian(double v, double μ, double σ, double MaxX, RectangleF R)     //returns point within rectangle for gaussian curve
		{
			double yy = Gaussian(v, μ, σ);
			double maxH = 1.0 / (SQRTΠ2 * σ);
			double dy = R.Height * yy / maxH;
			double v1 = (v - μ);
			double v2 = v1 / σ;
			double x = R.Left + R.Width * v / MaxX;
			double y = R.Bottom - dy;
			PointF pt = new PointF((float)x, (float)y);
			return pt;
		}
		//......................................................................
		/*
	   public static System.Windows.Point MapGaussian(double v, double μ, double σ, double Maxσ, System.Windows.Rect R)     //returns point within rectangle for gaussian curve
	   {
		   var σLength = R.Width / (2.0 * Maxσ);   //distance per 1 SD
		   var yy = CDASlib.Gaussian(v, μ, σ);
		   var maxH = 1.0 / (CDASlib.SQRTΠ2 * σ);
		   var dy = R.Height * yy / maxH;
		   var midX = R.Width / 2.0;
		   var v1 = (v - μ);
		   var v2 = v1 / σ;
		   var x = R.Left + midX + (((v - μ) / σ) * σLength);
		   var y = R.Bottom - dy;
		   return new System.Windows.Point(x, y);
	   }
		*/
		//.....................................................................
		public static double Getμg(double[] A)  //returns mean of log(A[i])
		{
			if (A == null || A.Length < 1) return 0.0;
			double sum = 0.0;
			for (int ix = 0; ix < A.Length; ix++) {
				if (A[ix] <= 0.0) continue;
				double v = Math.Log(A[ix]);
				sum += v;
			}
			sum /= A.Length;
			return sum; //ln of log-Normal mean
		}
		//......................................................................
		public static double Getσg(double[] A, double μg)   //returns standard deviation of log(A[i])
		{
			double sum = 0.0;
			if (A == null || A.Length < 1) return 0.0;
			for (int ix = 0; ix < A.Length; ix++) {
				if (A[ix] <= 0) continue;
				double v = (Math.Log(A[ix]) - μg);
				sum += (v * v);
			}
			sum = sum / A.Length;
			return Math.Sqrt(sum);
		}
		//....................................................................
		public static double LogNormal(double x, double μg, double σg)   //retrurns the PDF for LogNormal distribution, where μg is the ln(Geometric mean)...
		{
			if (x <= 0.0) return 0.0;
			double K = x * σg * SQRTΠ2;
			double v = Math.Log(x) - μg;
			double vv = v * v;
			double d1 = 2.0 * σg * σg;
			double d2 = Math.Exp(-vv / d1);
			double ans = d2 / K;
			return ans;
		}
		//......................................................................
		public static double GetGeometricMean(double[] A)  //returns mean of Σ( ln(x)	 )/N
		{
			if (A == null || A.Length < 1) return 0.0;
			double sum = 0.0, ans; int count = 0, ix;
			for (ix = 0; ix < A.Length; ix++) {
				if (A[ix] <= 0) continue;
				count += 1;
				sum += Math.Log(A[ix]);
			}
			ans = (count > 0) ? Math.Exp(sum / count) : 0.0;
			return ans;
		}
		//.....................................................................
		public static double GetGeometricSD(double[] A, double GM) {
			if (A == null || A.Length < 1) return 0.0;
			double logsum = 0.0, ans = 0.0, d0, v;
			int ix;
			for (ix = 0; ix < A.Length; ix++) {
				v = Math.Log(A[ix] / GM);
				logsum += (v * v);
			}
			d0 = Math.Sqrt(logsum / A.Length);
			ans = Math.Exp(d0);
			return ans;
		}

		//....................................................................
		public static PointF MapLogNormal(double v, double μg, double σg, double MaxX, double MaxY, RectangleF R)     //returns point within rectangle for gaussian curve
		{
			double x, y;
			x = R.Left + R.Width * v / MaxX;

			double yy = LogNormal(v, μg, σg);
			double dy = R.Height * yy / MaxY;   //get percentage of height from max
			y = R.Bottom - dy;
			if (y < R.Top) y = R.Top;

			PointF pt = new PointF((float)x, (float)y);
			return pt;
		}

		//....................................................................
		public static readonly string sGreek = "ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩαβγδεζηθικλμνξοπρςστυφχψω";
		//....................................................................
		public static double Π { get { return Math.PI; } }
		public static double Π2 { get { return Math.PI * 2; } }
		public static double SQRTΠ2 { get { return Math.Sqrt(Π2); } }
		//....................................................................
		public static bool GetShareList(string smachine, out string[] SA) {
			SA = null; string se;
			try {
				string sm = smachine.Trim('\\');
				string path = string.Format(@"\\{0}\root\cimv2", sm);
				var query = "select * from win32_share";
				var worker = new ManagementObjectSearcher(path, query);
				var shares = worker.Get();
				if (shares.Count == 0) return true;
				SA = new string[shares.Count];
				int i = 0;
				foreach (ManagementObject share in shares) {
					SA[i++] = (string)share["Name"];
				}
			} catch (Exception exp) {
				SA = null;
				se = "Error: " + exp.Message;
				return false;
			}
			return true;
		}
		//.....................................................................
		public static Font OptimumFont(Graphics g, String sText, ref SizeF sz, string sFamily, System.Drawing.FontStyle Style) //returns font that fits text
		{
			if (sz.Width * sz.Height <= 1.0) return null;
			if (sText.Length < 1) return null;

			float mxwidth = sz.Width, mxhgt = sz.Height, em = mxhgt; Font font = null;

			do {
				font = new Font(sFamily, em, Style, GraphicsUnit.Pixel);
				if (font == null) return null;
				sz = g.MeasureString(sText, font);
				if ((sz.Height <= mxhgt) && (sz.Width < mxwidth)) return font;
				font.Dispose();
				em *= 0.95f;
			} while (em > 1.0f);

			return null;

		}
	}

	///////////////////////////////////////////////////////////////////////
	public class TSort  //implements a tree sort algorithm
	{
		public delegate int QCompare(int ix, int iy);
		public delegate void QSwap(int ix, int iy);

		public event QCompare Compare;  //instance of delegate
		public event QSwap Swap;            //instance of delegate
											//=================================================================
		private void Short(int len)//short sort
		{
			int i, j;
			for (i = 0; i < len - 1; i++) {
				for (j = i + 1; j < len; j++) if (Compare(i, j) > 0) Swap(i, j);
			}
		}
		//=================================================================
		public void Sort(int len) //Implements a Tree Sort
		{
			if (len < 2) return;
			if (len < 6) { Short(len); return; }
			//grow tree 
			int m = (len - 1) / 2;
			while (m >= 0) {
				Downswap(m, len);
				m -= 1;
			}
			//cut tree down from root
			while (len > 0) {
				Swap(0, --len); //swap 1st, last.
				Downswap(0, len);
			}
		}
		//==================================================================
		private void Downswap(int m, int len) {
			int n = m * 2 + 1, p = n + 1;
			if (n >= len) return;
			if ((p < len) && (Compare(n, n + 1) < 0)) n += 1;
			if (Compare(m, n) > 0) return;
			Swap(m, n);
			Downswap(n, len);
		}
		//..................................................................
		public static string sWL(int ix, int[] A) {
			return '[' + ix.ToString("n0") + "]=" + A[ix].ToString("n0");
		}
		public static void WL(int ix, int[] A) { Console.WriteLine(sWL(ix, A)); }
		public static string sSWL(int ix, int iy, int[] A) {
			return '\t' + sWL(ix, A) + "<--->" + sWL(iy, A);
		}
		public static void SWL(int ix, int iy, int[] A) { Console.WriteLine(sSWL(ix, iy, A)); }
		//..................................................................
		public static bool isTree(int[] T) {
			int len = T.Length, i, m;
			for (i = 0; i < T.Length; i++) {
				m = i * 2 + 1;
				if (m >= T.Length) return true;
				if (T[i] < T[m++]) return false;
				if (m < T.Length && T[i] < T[m]) return false;
			}
			return true;
		}
		//...................................................................
		public bool IsSorted(int len) {
			for (int ix = 0; ix < len - 2; ix++) if (Compare(ix, ix + 1) > 0) return false;
			return true;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	#region NetworkBrowser CLASS
	/// <summary>
	/// Provides a mechanism for supplying
	// a list of all PC names in the local network.
	/// This collection of PC names is used in the form 
	/// 


	/// This class makes use of a DllImport instruction.
	/// The purpose of which is as follows:
	/// When a DllImport declaration is made
	/// in managed code (C#) it is a call to a legacy
	/// unmanaged code module, normally
	/// a C++ Dynamic Link Library. These C++ Dll's are
	/// usually part of the operating system API,
	/// or some other vendors API, and must be 
	/// used to carry out operations that are not
	/// native within the managed code C# framework. 
	/// This is fairly normal within the windows world.
	/// The only thing that needs careful consideration
	/// is the construction of the correct type of STRUCTS,
	/// object pointers, and attribute markers,
	/// which all contribute to making the link
	/// between managed (C#) and unmanaged code (C++)
	/// more seamless
	/// 

	/// This class makes use of the following Dll calls
	/// <list type="bullet">
	/// <item>
	/// <description> Netapi32.dll : NetServerEnum,
	/// The NetServerEnum function lists all servers
	/// of the specified type that are visible in
	/// a domain. For example, an application can call 
	/// NetServerEnum to list all domain controllers
	/// only or all SQL servers only.
	/// You can combine bit masks to list several
	/// types. For example, a value of 0x00000003 
	/// combines the bit masks for SV_TYPE_WORKSTATION
	/// (0x00000001) and SV_TYPE_SERVER (0x00000002).
	/// </description>
	/// </item>
	/// <item>
	/// <description> Netapi32.dll : NetApiBufferFree,
	/// The NetApiBufferFree function frees 
	/// the memory that the NetApiBufferAllocate
	/// function allocates. Call NetApiBufferFree 
	/// to free the memory that other network
	/// management functions return.</description>
	/// </item>
	/// </list>
	/// </summary>
	public sealed class NetworkBrowser {
		#region Dll Imports

		//declare the Netapi32 : NetServerEnum method import
		[DllImport("Netapi32", CharSet = CharSet.Auto,
		SetLastError = true),
		SuppressUnmanagedCodeSecurityAttribute]

		/// <summary>
		/// Netapi32.dll : The NetServerEnum function lists all servers
		/// of the specified type that are
		/// visible in a domain. For example, an 
		/// application can call NetServerEnum
		/// to list all domain controllers only
		/// or all SQL servers only.
		/// You can combine bit masks to list
		/// several types. For example, a value 
		/// of 0x00000003  combines the bit
		/// masks for SV_TYPE_WORKSTATION 
		/// (0x00000001) and SV_TYPE_SERVER (0x00000002)
		/// </summary>
		public static extern int NetServerEnum(
			 string ServerNane, // must be null
			 int dwLevel,
			 ref IntPtr pBuf,
			 int dwPrefMaxLen,
			 out int dwEntriesRead,
			 out int dwTotalEntries,
			 int dwServerType,
			 string domain, // null for login domain
			 out int dwResumeHandle
			 );

		//declare the Netapi32 : NetApiBufferFree method import
		[DllImport("Netapi32", SetLastError = true),
		SuppressUnmanagedCodeSecurityAttribute]

		/// <summary>
		/// Netapi32.dll : The NetApiBufferFree function frees 
		/// the memory that the NetApiBufferAllocate function allocates. 
		/// Call NetApiBufferFree to free
		/// the memory that other network 
		/// management functions return.
		/// </summary>
		public static extern int NetApiBufferFree(
			 IntPtr pBuf);

		//create a _SERVER_INFO_100 STRUCTURE
		[StructLayout(LayoutKind.Sequential)]
		public struct _SERVER_INFO_100 {
			internal int sv100_platform_id;
			[MarshalAs(UnmanagedType.LPWStr)]
			internal string sv100_name;
		}
		#endregion
		#region Public Constructor
		/// <SUMMARY>
		/// Constructor, simply creates a new NetworkBrowser object
		/// </SUMMARY>
		public NetworkBrowser() {

		}
		#endregion
		#region Public Methods
		/// <summary>
		/// Uses the DllImport : NetServerEnum
		/// with all its required parameters
		/// (see http://msdn.microsoft.com/library/default.asp?
		///      url=/library/en-us/netmgmt/netmgmt/netserverenum.asp
		/// for full details or method signature) to
		/// retrieve a list of domain SV_TYPE_WORKSTATION
		/// and SV_TYPE_SERVER PC's
		/// </summary>
		/// <returns>Arraylist that represents
		/// all the SV_TYPE_WORKSTATION and SV_TYPE_SERVER
		/// PC's in the Domain</returns>
		public ArrayList getNetworkComputers() {
			//local fields
			ArrayList networkComputers = new ArrayList();

			const int MAX_PREFERRED_LENGTH = -1;
			int SV_TYPE_WORKSTATION = 1;
			int SV_TYPE_SERVER = 2;
			IntPtr buffer = IntPtr.Zero;
			IntPtr tmpBuffer = IntPtr.Zero;
			int entriesRead = 0;
			int totalEntries = 0;
			int resHandle = 0;
			int sizeofINFO = Marshal.SizeOf(typeof(_SERVER_INFO_100));


			try {
				//call the DllImport : NetServerEnum 
				//with all its required parameters
				//see http://msdn.microsoft.com/library/
				//default.asp?url=/library/en-us/netmgmt/netmgmt/netserverenum.asp
				//for full details of method signature

				int ret = NetServerEnum(null, 100, ref buffer,
					 MAX_PREFERRED_LENGTH,
					 out entriesRead,
					 out totalEntries,
					 SV_TYPE_WORKSTATION | SV_TYPE_SERVER,
					 null,
					 out resHandle);

				//if the returned with a NERR_Success 
				//(C++ term), =0 for C#
				if (ret == 0) {
					//loop through all SV_TYPE_WORKSTATION 
					//and SV_TYPE_SERVER PC's
					for (int i = 0; i < totalEntries; i++) {
						//get pointer to, Pointer to the 
						//buffer that received the data from
						//the call to NetServerEnum. 
						//Must ensure to use correct size of 
						//STRUCTURE to ensure correct 
						//location in memory is pointed to
						tmpBuffer = new IntPtr((int)buffer + (i * sizeofINFO));
						//Have now got a pointer to the list 
						//of SV_TYPE_WORKSTATION and 
						//SV_TYPE_SERVER PC's, which is unmanaged memory
						//Needs to Marshal data from an 
						//unmanaged block of memory to a 
						//managed object, again using 
						//STRUCTURE to ensure the correct data
						//is marshalled 
						_SERVER_INFO_100 svrInfo = (_SERVER_INFO_100)
							 Marshal.PtrToStructure(tmpBuffer,
										typeof(_SERVER_INFO_100));

						//add the PC names to the ArrayList
						networkComputers.Add(svrInfo.sv100_name);
					}
				}
			} catch (Exception ex) {
				MessageBox.Show("Problem with acessing " +
						"network computers in NetworkBrowser " +
						"\r\n\r\n\r\n" + ex.Message,
					 "Error");
				return null;
			} finally {
				//The NetApiBufferFree function frees 
				//the memory that the 
				//NetApiBufferAllocate function allocates
				NetApiBufferFree(buffer);
			}
			//return entries found
			return networkComputers;

		}
		#endregion
	}
	#endregion

}
