using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SNL {
	// ////////////////////////////////////////////////////////////////////
	[Serializable]
	public class CVector {
		private string sname;
		private double[] v = null;
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		public CVector(String Name, int N) {
			if (N <= 0) N = 1;
			sname = Name;
			v = new double[N];
			Clear();
		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		public CVector(params double[] Y) {
			v = new double[Y.Length];
			for (int i = 0; i < Y.Length; i++) v[i] = Y[i];
			sname = "V";
		}
		//================================================================
		public bool isValid {
			get {
				for (int ix = 0; ix < v.Length; ix++) if (double.IsInfinity(v[ix]) || double.IsNaN(v[ix])) return false;
				return true;
			}
		}
		//================================================================
		public bool Fill(params double[] Y) {
			if (Y == null) return false;
			if (v.Length < Y.Length) v = new double[Y.Length];
			for (int i = 0; i < Y.Length; i++) v[i] = Y[i];
			return true;
		}
		//================================================================
		public string sName { get { return sname; } set { sname = value; } }
		public void Clear() { for (int i = 0; i < v.Length; i++) v[i] = 0.0; }
		public double Sum => v.Sum();
		public double Min => v.Min();
		public double Max => v.Max();
		public int Length => v.Length;
		public double[] DoubleArray => v;
		public double this[int i] { get { return (i < v.Length ? v[i] : 0.0); } set { v[i] = value; } }
		public string sV(string sf) { return sV(v, sf); }
		//===================================================================
		public double Mag {
			get {
				int i, n = v.Length;
				double mag = 0.0;
				for (i = 0; i < n; i++) mag += v[i] * v[i];
				return Math.Sqrt(mag);
			}
		}
		//===================================================================
		public string sVall {
			get {
				int i, n = v.Length;
				StringBuilder sb = new StringBuilder(100);
				sb.Clear();
				sb.Append("[");
				for (i = 0; i < n - 1; i++) sb.Append(v[i].ToString("f1") + ", ");
				sb.Append(v[i] + ']');
				return sb.ToString();
			}
		}
		public override string ToString() => sname + '[' + Length.ToString("d0") + ']' + " => " + sV("f3");
		//.................................................................
		public static CVector operator +(CVector V1, CVector V2) {
			if (V1 is null || V2 == null) return null;
			if (V1.Length != V2.Length) return null;
			CVector VPlus = new CVector("VPlus", V1.Length);
			for (int ix = 0; ix < VPlus.Length; ix++) VPlus[ix] = V1[ix] + V2[ix];
			return VPlus;
		}
		//.................................................................
		public static CVector operator -(CVector V1, CVector V2) {
			if (V1 is null || V2 == null) return null;
			if (V1.Length != V2.Length) return null;
			CVector VMinus = new CVector("VMinus", V1.Length);
			for (int ix = 0; ix < VMinus.Length; ix++) VMinus[ix] = V1[ix] - V2[ix];
			return VMinus;
		}
		//.................................................................
		public static bool operator ==(CVector V1, CVector V2) {
			if (V1 is null) return (V2 is null ? true : false);
			if (V2 is null) return false;
			if (V1.Length != V2.Length) return false;

			for (int ix = 0; ix < V1.Length; ix++) if (V1[ix] != V2[ix]) return false;
			return true;
		}
		public static bool operator !=(CVector V1, CVector V2) => !(V1 == V2);
		//.................................................................
		public static string sV(double[] V, string sf) {
			if (V == null) return " null ";
			string s; StringBuilder sb = new StringBuilder(1000);
			int ix;
			sb.Append("[");
			for (ix = 0; ix < V.Length; ix++) {
				s = (ix > 0 ? ", " : "  ") + V[ix].ToString(sf);
				sb.Append(s);
			}
			sb.Append("]");
			s = sb.ToString();
			return s;
		}
		//................................................................
		public static double DotProduct(CVector v1, CVector v2) {
			if (v1.Length != v2.Length) return double.NaN;
			double sum = 0;
			double den = v1.Mag * v2.Mag;
			for (int ix = 0; ix < v1.Length; ix++) sum += (v1[ix] * v2[ix]);
			return sum / den;
		}

		//..................................................................
		public static CVector operator *(double v, CVector v2) {
			if (v2 == null || v2.Length == 0) return null;
			CVector v1 = new CVector(v2.sName + "A", v2.Length);
			for (int ix = 0; ix < v2.Length; ix++) v1[ix] = v2[ix] * v;
			return v1;
		}
		//....................................................................
		public static CVector operator *(CVector v1, CVector v2) {
			if (v1 == null || v2 == null || v1.Length != v2.Length) return null;
			CVector ans = new CVector(v1.sname + v2.sname, v1.Length);
			for (int ix = 0; ix < v1.Length; ix++) ans[ix] = v1[ix] * v2[ix];
			return ans;
		}
		//....................................................................
		public static double ACos(CVector v1, CVector v2) => Math.Acos(CVector.DotProduct(v1, v2));

		//==================================================================
		public override bool Equals(object obj) => this == (CVector)obj;
		public override int GetHashCode() => base.GetHashCode();
	}
	// ////////////////////////////////////////////////////////////////////
	[Serializable]
	public class CMatrix {
		private const char SQ = (char)178;

		private string sname = "";
		private int nrow;
		private int ncol;

		double[] matrix = null;
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		public CMatrix(string Name, int Row, int Col)/*  M(row,col) = matrix(row*ncol + col) */ {
			sname = Name;
			nrow = Row;
			ncol = Col;
			matrix = new double[nrow * ncol];
			Clear();
		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		public CMatrix(string Name, CMatrix M)/*copy matrix M*/{
			sname = Name;
			nrow = M.nrow;
			ncol = M.ncol;
			matrix = new double[nrow * ncol];
			for (int ir = 0; ir < nrow; ir++) {
				for (int ic = 0; ic < ncol; ic++) this[ir, ic] = M[ir, ic];
			}
		}
		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		public CMatrix(string Name, double[,] D) {
			sname = Name;
			nrow = D.GetLength(0);
			ncol = D.GetLength(1);
			matrix = new double[nrow * ncol];
			for (int i = 0; i < nrow; i++) {
				for (int j = 0; j < ncol; j++) matrix[ncol * i + j] = D[i, j];
			}
		}
		//=================================================================
		public override string ToString() => sname + "[r:" + nrow.ToString("n0") + ", c:" + ncol.ToString("n0") + ']';
		//=================================================================
		public bool isSquare => nrow == ncol;
		//=================================================================
		public bool isSingular(double ε) /*compare each row against all others*/{
		
				int i, j, n = nRows; double dot = 0.0;
				for (i = 0; i < n; i++) {
					CVector Vi = Row(i, "vi");
					for (j = i + 1; j < n; j++) {
						CVector Vj = Row(j, "vj");
						dot = CVector.DotProduct(Vi, Vj);
						if (dot < ε) return true;
					}
				}
				return false;
		}
		//==================================================================
		public int[] RemoveRow(int row,int[] I) {
			if (row >= nrow) return I;
			int i, j, n = (nrow-1) * ncol;
			i = row * ncol; //start
			while (i < n) {
				j = i + ncol;
				matrix[i] = matrix[j];
				matrix[j] = 0;
				i += 1;
			}
			j = row;
			do { I[j] = I[j + 1]; j += 1; } while (j < I.Length - 1);
			nrow -= 1;
			return I;
		}
		//==================================================================
		public int[] MakeNonSingular(double ε) /*returns and index */{
			double dot = 0.0;
			int i,row;
			int[] I = new int[nRows];
			for (i = 0; i < nRows; i++) I[i] = i;
			for (row = 0; row < nrow - 1; row++) {
				if (I[row] < 0) continue;
				CVector VRow = new CVector(Row(row));
				for (int r = row + 1; r < nrow; r++) {
					if (I[r] < 0) continue;
					CVector Vr = new CVector(Row(r));
					dot = CVector.DotProduct(VRow, Vr);
					if (dot < ε) I = RemoveRow(r,I);
				}
			}
			int[] II = new int[nRows];
			for (i = 0; i < II.Length; i++) II[i] = I[i];
			return II;
		}
		//==================================================================
		public CVector Row(int row, string sN) {
			if (row < 0 || row >= nRows) return null;
			CVector V = new CVector(sN, nCols);
			V.Fill(Row(row));
			return V;
		}
		//=================================================================
		public void Clear() { for (int i = 0; i < matrix.Length; i++) matrix[i] = 0; }
		public double Sum => matrix.Sum();
		public double Min => matrix.Min();
		public double Max => matrix.Max();
		//=================================================================
		public double this[int r, int c] { get { return matrix[r * ncol + c]; } set { matrix[r * ncol + c] = value; } }
		//=================================================================
		public string sName => sname;
		public bool isEmpty => (nrow * ncol == 0);
		//===================================================================
		public bool SetRow(int Row, params double[] V) {
			if (V.Length != ncol || Row >= nRows) return false;
			int rx = Row * ncol;
			for (int ix = 0; ix < ncol; ix++) matrix[rx + ix] = V[ix];
			return true;
		}
		//====================================================================
		public bool SetRow(int Row, CVector V) {
			if (V.Length != ncol || Row >= nRows) return false;
			int ix, rx = Row * ncol;
			for (ix = 0; ix < ncol; ix++) matrix[rx + ix] = V[ix];
			return true;
		}
		//==================================================================
		public bool SetCol(int Col, params double[] V) {
			if (V.Length != nrow || Col >= nCols) return false;
			int cx = Col * nrow, ix;
			for (ix = 0; ix < nrow; ix++) matrix[cx + ix] = V[ix];
			return true;
		}
		//================================================================
		public double[] Row(int row) {
			if (row >= nrow) return null;
			double[] V = new double[ncol];
			for (int ix = 0; ix < ncol; ix++) V[ix] = this[row, ix];
			return V;
		}
		//=================================================================
		public double[] Col(int col) {
			if (col >= ncol) return null;
			double[] V = new double[nrow];
			for (int ix = 0; ix < nrow; ix++) V[ix] = this[ix, col];
			return V;
		}
		//================================================================
		public int nRows => nrow;
		public int nCols => ncol;
		//================================================================
		public CMatrix Transpose {
			get {
				CMatrix M = new CMatrix(sname + "\u036d", ncol, nrow);
				int ir, ic;
				for (ir = 0; ir < nrow; ir++) {
					for (ic = 0; ic < ncol; ic++) {
						M[ic, ir] = this[ir, ic];
					}
				}
				return M;
			}
		}
		//=================================================================
		public string sRow(int row, string sf) {
			double[] V = Row(row);
			return sV(V, "0.00");
		}
		//=================================================================
		public double Determinant {
			get {
				if (nRows != nCols) return double.NaN;
				double v1, v2, d1 = 0.0, d2 = 0.0, ans = 0.0;
				int p, n = nCols, i, ic, ir, imc;
				for (p = 0; p < n; p++) {
					v1 = 1.0;
					v2 = 1.0;
					for (ir = 0, ic = p, imc = inn(n - p, n), i = 0; i < n; i++) {
						v1 *= this[ir, ic];
						v2 *= this[ir, imc];
						ic = imn(ic + 1, n);
						ir = imn(ir + 1, n);
						imc = inn(imc, n);
					}
					d1 += v1;
					d2 += v2;
				}
				ans = d1 - d2;
				return ans;
			}
		}
		//=================================================================
		public CMatrix Inverse {
			get {
				if (nrow != ncol) return null;
				CMatrix I = new CMatrix(sName + "\u0304", this);
				int p, row, col, n = nrow; double v, d = 1.0;

				for (p = 0; p < n; p++) {
					v = I[p, p];
					if (v == 0) return null;
					d *= v;
					for (row = 0; row < n; row++)/*pivot row*/ {
						if (row == p) continue;
						I[row, p] /= v;
					}
					for (col = 0; col < n; col++) /*pivot col*/{
						if (col == p) continue;
						I[p, col] /= (-v);
					}
					for (row = 0; row < n; row++) {
						if (row == p) continue;
						for (col = 0; col < n; col++) {
							if (col == p) continue;
							I[row, col] = I[row, col] + I[p, col] * I[row, p];
						}
					}
					I[p, p] = 1.0 / I[p, p];
				}
				return I;
			}
		}
		//=================================================================
		public CMatrix InV {
			get {
				if (nrow != ncol) return null;
				CMatrix I = new CMatrix(sName + "-I", this);
				double[] b = new double[ncol];

				int ic, ir, k; double ratio = 0.0;
				for (k = 0; k < b.Length; k++) b[k] = 0.0;
				for (ic = 0; ic < ncol - 1; ic++) {
					for (ir = ic; ir < ncol; ir++) {
						ratio = I[ir, ic] / I[ic, ic];
						for (k = ic; k < ncol; k++) {
							I[ir, k] -= (ratio * I[ic, k]);
							b[ir] -= (ratio * b[ic]);
						}
					}
				}

				return I;
			}
		}
		//.................................................................
		// CPP Program to decompose a matrix into 
		// lower and upper traingular matrix 
		public static bool LU(CMatrix M, out CMatrix Lower, out CMatrix Upper) {
			Lower = null;
			Upper = null;
			if (M == null || M.nCols != M.nRows) return false;
			int i, j, k, N = M.ncol;
			double sum;
			Lower = new CMatrix(M.sname + "-Lower", N, N);
			Upper = new CMatrix(M.sname + "-Upper", N, N);
			Lower.Clear();
			Upper.Clear();


			// Decomposing matrix into Upper and Lower triangular matrix
			for (i = 0; i < N; i++) {
				//upper Triangular
				for (k = i; k < N; k++) {
					sum = 0.0;
					for (j = 0; j < i; j++) sum += Lower[i, j] * Upper[j, k];
					Upper[i, k] = M[i, k] - sum;
				}
				//Lower Triangular
				for (k = i; k < N; k++) {
					if (i == k) Lower[i, i] = 1.0;     //diagonal set to 1
					else {
						sum = 0.0;
						for (j = 0; j < i; j++) sum += Lower[k, j] * Upper[j, i];
						Lower[k, i] = (M[k, i] - sum) / Upper[i, i];
					}
				}
			}
			return true;
		}
		//..........................................................
		public static CVector SolveLower(CMatrix Lower, CVector B, string sName) {
			if (Lower == null || B == null || !Lower.isSquare || Lower.nrow != B.Length) return null;
			CVector Y = new CVector(sName, B.Length);
			int col, row; double sum;
			double v = B[0];
			Y[0] = v;
			for (row = 1; row < Lower.nRows; row++) {
				sum = 0.0;
				for (col = 0; col < row; col++) sum += (Lower[row, col] * Y[col]);
				Y[row] = B[row] - sum;
			}
			return Y;
		}
		//..........................................................
		public static CVector SolveUpper(CMatrix Upper, CVector Z, string sName) {
			if (Upper == null || !Upper.isSquare || Upper.nrow != Z.Length) return null;
			CVector C = new CVector(sName, Z.Length);
			int row, col, n = Z.Length - 1; double sum, v;
			v = Z[n];
			for (row = n; row >= 0; row--) {
				sum = 0.0;
				for (col = row + 1; col <= n; col++) sum += Upper[row, col] * C[col];
				v = Z[row] - sum;
				C[row] = v / Upper[row, row];
			}

			return C;
		}
		//..........................................................
		private static double[] SolveLU(CMatrix M, double[] rightPart, out CMatrix LU) {
			// decomposition of matrix
			LU = null;
			if (M == null || M.isEmpty) return null;
			int n = M.nCols;

			double[,] lu = new double[n, n];
			double sum = 0;
			for (int i = 0; i < n; i++) {
				for (int j = i; j < n; j++) {
					sum = 0;
					for (int k = 0; k < i; k++)
						sum += lu[i, k] * lu[k, j];
					lu[i, j] = M[i, j] - sum;
				}
				for (int j = i + 1; j < n; j++) {
					sum = 0;
					for (int k = 0; k < i; k++)
						sum += lu[j, k] * lu[k, i];
					lu[j, i] = (1 / lu[i, i]) * (M[j, i] - sum);
				}
			}
			LU = new CMatrix("LU", lu);

			// lu = L+U-I
			// find solution of Ly = b
			double[] y = new double[n];
			for (int i = 0; i < n; i++) {
				sum = 0;
				for (int k = 0; k < i; k++)
					sum += lu[i, k] * y[k];
				y[i] = rightPart[i] - sum;
			}
			// find solution of Ux = y
			double[] x = new double[n];
			for (int i = n - 1; i >= 0; i--) {
				sum = 0;
				for (int k = i + 1; k < n; k++)
					sum += lu[i, k] * x[k];
				x[i] = (1 / lu[i, i]) * (y[i] - sum);
			}
			return x;
		}
		//...................................................................
		public static CVector LeastSquareSolution(CMatrix M, CVector Y)/*Multivariable LSQ*/{
			if (M == null || M.nCols * M.nRows == 0) return null;
			int mrow, row, col, N = M.nCols + 1;
			double[] dv = null;
			CVector C = new CVector("C", N);    //internal
			CMatrix SQ = new CMatrix("SQ", N, N);

			//generate internal vector
			C[0] = Y.Sum;
			for (col = 1; col < C.Length; col++) {
				dv = M.Col(col - 1);
				for (row = 0; row < M.nRows; row++) C[col] += dv[row] * Y[row];
			}
			//generate internal least squares matrix
			//do  column 0, row 0
			SQ[0, 0] = M.nRows;
			for (col = 0, row = 1; row < SQ.nRows; row++) {
				dv = M.Col(col++);
				SQ[col, 0] = dv.Sum();
				SQ[0, col] = SQ[col, 0];
			}
			for (row = 1; row < SQ.nCols; row++) {
				for (col = 1; col < SQ.nCols; col++) {
					SQ[col, row] = 0;
					for (mrow = 0; mrow < M.nRows; mrow++) {
						double dd = SQ[row, col];
						double v1 = M[mrow, row - 1];
						double v2 = M[mrow, col - 1];
						SQ[col, row] += v1 * v2;
					}
				}
			}
			CMatrix Lower = null, Upper = null;
			if (!CMatrix.LU(SQ, out Lower, out Upper)) return null;

			CVector Z = CMatrix.SolveLower(Lower, C, "Z");
			CVector B = CMatrix.SolveUpper(Upper, Z, "β");
			return B;
		}
		//.................................................................
		public static CMatrix Identity(string Name, int n) {
			if (n == 0) return null;
			CMatrix M = new CMatrix(Name, n, n);
			for (int ir = 0; ir < n; ir++) {
				for (int ic = 0; ic < n; ic++) if (ic == ir) M[ic, ir] = 1.0;
			}
			return M;
		}
		//...................................................................
		public static CMatrix operator *(CMatrix A, CMatrix B) {
			if (A.nCols != B.nRows) return null;
			double[] V1 = null, V2 = null;
			CMatrix C = new CMatrix(A.sname + "_*_" + B.sname, A.nRows, B.nCols);
			for (int irow = 0; irow < A.nRows; irow++) {
				for (int icol = 0; icol < B.nCols; icol++) {
					V1 = A.Row(irow);
					V2 = B.Col(icol);
					C[irow, icol] = VdotV(V1, V2);
				}
			}
			return C;
		}
		//................................................................
		public static CVector operator *(CMatrix A, CVector V) {
			if (A.nCols != V.Length) return null;
			CVector B = new CVector(A.sName + "_*_" + V.sName, A.nRows);
			int row, col; double v;
			for (row = 0; row < A.nRows; row++) {
				v = 0;
				for (col = 0; col < A.nCols; col++) {
					v += A[row, col] * V[col];
				}
				B[row] = v;
			}
			return B;
		}
		//................................................................
		private int imn(int iv, int md) => (iv % md);
		private int inn(int iv, int md) => (iv <= 0 ? md - 1 : iv - 1);
		//................................................................
		public static CMatrix Multiply(CMatrix A, CMatrix B) => A * B;
		//.................................................................
		public static bool operator ==(CMatrix M1, CMatrix M2) {
			if (M1 is null) {
				if (M2 is null) return true;
				return false;
			}
			if (M2 is null) return false;
			if (M1.nRows != M2.nRows) return false;
			if (M1.nCols != M2.nCols) return false;
			for (int row = 0; row < M1.nRows; row++) {
				for (int col = 0; col < M1.nCols; col++) if (M1[row, col] != M2[row, col]) return false;
			}
			return true;
		}
		public static bool operator !=(CMatrix M1, CMatrix M2) => !(M1 == M2);
		public override bool Equals(object obj) => this == (CMatrix)obj;
		public override int GetHashCode() => base.GetHashCode();
		//.................................................................
		public static double VdotV(double[] V1, double[] V2) {
			double ans = 0.0;
			if (V1.Length != V2.Length) return 0;
			for (int ix = 0; ix < V1.Length; ix++) ans += V1[ix] * V2[ix];
			return ans;
		}
		//.................................................................
		public static string sV(double[] V, string sf) {
			if (V == null) return " null ";
			string s; StringBuilder sb = new StringBuilder(1000);
			int ix;
			sb.Append("[");
			for (ix = 0; ix < V.Length; ix++) {
				s = (ix > 0 ? ", " : "  ") + V[ix].ToString(sf);
				sb.Append(s);
			}
			sb.Append("]");
			s = sb.ToString();
			return s;
		}
	}
}
