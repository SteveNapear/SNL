using System;
using System.Runtime.Remoting.Metadata.W3cXsd2001;
using System.Windows;
using System.Windows.Media;
using System.Windows.Navigation;
using Color = System.Windows.Media.Color;

namespace SNL {
	//////////////////////////////////////////////////////////////////////
	[Serializable]
	public class CTriangle {
		private string name = string.Empty;
		protected Point[] A = new Point[3];
		private static Color white = Colors.White;

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		public CTriangle(Point L, Point R, Point up) {
			A[0] = L; A[1] = R; A[2] = up;
		}
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		public CTriangle() { A[0] = new Point(0, 0); A[1] = new Point(0, 0); A[2] = new Point(0, 0); }
		//================================================================
		public string Name { get { return name; } set { name = value; } }
		//================================================================
		public Point this[int ix] => (ix < A.Length ? A[ix] : new Point(0, 0));
		public Point[] PointArray => A;
		public System.Drawing.Point[] DrawingArray {
			get {
				System.Drawing.Point[] DA = new System.Drawing.Point[A.Length];
				for (int ix = 0; ix < DA.Length; ix++) { DA[ix].X = (int)A[ix].X; DA[ix].Y = (int)A[ix].Y; }
				return DA;
			}
		}
		public Point GetPoint(int ix) => (ix < A.Length) ? A[ix] : new Point(0, 0);
		public Point ClonePoint(int ix) => (ix < A.Length) ? new Point(A[ix].X, A[ix].Y) : new Point(0, 0);
		//================================================================
		public double Circumference => DistPoints(A[0], A[1]) + DistPoints(A[0], A[2]) + DistPoints(A[1], A[2]);
		//================================================================
		public Rect BoundingRect {
			get {
				double lft, rgt, top, bot;
				lft = A[0].X; if (lft > A[1].X) lft = A[1].X; if (lft > A[2].X) lft = A[2].X;
				rgt = A[0].X; if (rgt < A[1].X) rgt = A[1].X; if (rgt < A[2].X) rgt = A[2].X;
				bot = A[0].Y; if (bot > A[1].Y) bot = A[1].Y; if (bot > A[2].Y) bot = A[2].Y;
				top = A[0].Y; if (top < A[1].Y) top = A[1].Y; if (top < A[2].Y) top = A[2].Y;
				Rect R = new Rect(lft, bot, rgt - lft, top - bot);
				return R;
			}
		}
		//================================================================
		public bool Contains(Point p) {
			double a =
				((A[1].Y - A[2].Y) * (p.X - A[2].X) + (A[2].X - A[1].X) * (p.Y - A[2].Y)) /
				((A[1].Y - A[2].Y) * (A[0].X - A[2].X) + (A[2].X - A[1].X) * (A[0].Y - A[2].Y));
			double b =
				((A[2].Y - A[0].Y) * (p.X - A[2].X) + (A[0].X - A[2].X) * (p.Y - A[2].Y)) /
				((A[1].Y - A[2].Y) * (A[0].X - A[2].X) + (A[2].X - A[1].X) * (A[0].Y - A[2].Y));
			double c = 1.0 - a - b;
			return (a >= 0 && a <= 1 && b >= 0 && b <= 1 && c >= 0 && c <= 1);

		}
		public bool IsInside(Point p) => Contains(p);
		//================================================================
		public override string ToString() => Name + ": [" + CSNUtil.sPoint(A[0]) + ',' + CSNUtil.sPoint(A[1]) + ',' + CSNUtil.sPoint(A[2]) + ']';
		//................................................................
		private double det(Vector v1, Vector v2) { return v1.X * v2.Y - v1.Y * v2.X; }
		//................................................................
		private double det(Point p1, Point p2) { return p1.X * p2.Y - p1.Y * p2.X; }
		//................................................................
		public static double DistPoints(Point fr, Point to) {
			double dx = to.X - fr.X, dy = to.Y - fr.Y;
			double d = Math.Sqrt(dx * dx + dy * dy);
			return d;
		}
		//................................................................
		public static double Deg2Rad(double deg) => Math.PI * deg / 180.0;
		public static double Rad2Deg(double rad) => rad * 180.0 / Math.PI;
		//................................................................
		public static Point Rotate(double radians, Point pt) {
			if (pt.X == 0 && pt.Y == 0) return pt;
			double cos = Math.Cos(radians);
			double sin = Math.Sin(radians);
			Point p = new Point();
			p.X = pt.X * cos - pt.Y * sin;
			p.Y = pt.Y * cos + pt.X * sin;
			return p;
		}
		//===============================================================
		public System.Windows.Media.Color ComputeMaxwellColor(Point pt) {
			if (!Contains(pt)) return Colors.Black;
			double side = DistPoints(A[0], A[1]);
			Rect B = BoundingRect;
			double dx = pt.X / side;
			double dy = pt.Y / side;
			Color c = new Color();
			c.A = 255;
			c.B = (byte)(255.0 * DistPoints(A[2], pt) / side); ;
			c.R = (byte)(255.0 * DistPoints(A[0], pt) / side); ;
			c.G = (byte)(255.0 * dy);
			return c;
		}
		//......................................................................
		public static Point RotatePoint(Point p, Point origin, double cos, double sin) {
			double x = p.X - origin.X, y = p.Y - origin.Y;
			Point newpoint = new Point(x * cos - y * sin, x * sin + y * cos);
			newpoint.X += origin.X; newpoint.Y += origin.Y;
			return newpoint;
		}
	}
	// ///////////////////////////////////////////////////////////////////////
	public class CIsosceles : CTriangle {
		private double side;
		private double ht;
		private static double C30 = 0.866025404, S30 = 0.5;
		private static double C120 = -S30, S120 = C30;
		private static double C240 = -S30, S240 = -C30;
		private Point cp;
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		private CIsosceles() { Initialize(); }

		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		public CIsosceles(double Sid, Point Top) {
			Initialize();
			side = Sid;
			double d = Math.Tan(Deg2Rad(30)) * 0.5;
			ht = side * C30;
			A[0] = new Point(Top.X - side / 2, Top.Y + ht);
			A[1] = Top;
			A[2] = new Point(Top.X + side / 2, A[0].Y);
			double Gy = A[1].Y;
			double t30 = Math.Tan(Math.PI / 6.0);
			double dy = A[0].Y - t30 * side / 2.0;
			double y = A[1].Y + dy;
			cp = new Point(Top.X, y);
		}
		//====================================================================
		private void Initialize() {
			C30 = Math.Cos(Math.PI / 6.0);
			S30 = 0.5;
			C120 = -S30;
			S120 = C30;
			C240 = -S30;
			S240 = -C30;
		}
		//===================================================================
		public double Side => side;
		public double Height => ht;
		public Point Center => cp;
		public double Area => side * ht / 2.0;
		//===================================================================
		public CIsosceles Rotate(double Radians) /*Rotate this CIsosceles, return rotated new object*/{
			double cos = Math.Cos(Radians), sin = Math.Sin(Radians); int ix;
			CIsosceles RT = new CIsosceles();
			for (ix = 0; ix < A.Length; ix++) RT.A[ix] = RotatePoint(A[ix], cp, cos, sin);
			Rect rB = RT.BoundingRect;
			RT.cp = new Point(rB.Left + rB.Width / 2, rB.Top + rB.Height / 2);
			double dx = RT.cp.X - cp.X, dy = RT.cp.Y - cp.Y;
			for (ix = 0; ix < A.Length; ix++) { RT.A[ix].X -= dx; RT.A[ix].Y -= dy; }
			RT.cp = cp;
			RT.side = side;
			RT.ht = ht;
			return RT;
		}
		//===================================================================
		public byte Red(Point pt) {
			if (!Contains(pt)) return 0;
			double A0Y = pt.X * S120 + pt.Y * C120;
			double A1Y = A[1].X * S120 + A[1].Y * C120;
			double v = A1Y - A0Y;
			double dv = v / ht;
			//			Point nrp = RelativePoint(np);
			byte r = (byte)(255.0 * dv);
			return r;
		}
		public byte Green(Point pt) {
			if (!Contains(pt)) return 0;
			double dg = deltaY(A[1].Y, pt.Y);
			byte green = (byte)(255 * (1.0 - dg));
			return green;
		}
		public byte Blue(Point pt) {
			if (!Contains(pt)) return 0;
			double A0Y = pt.X * S240 + pt.Y * C240;
			double A1Y = A[1].X * S240 + A[1].Y * C240;
			double v = A1Y - A0Y;
			double dv = v / ht;
			byte r = (byte)(255.0 * dv);
			return r;
		}
		//===============================================================
		public Color color(Point pt) => Color.FromArgb(255, Red(pt), Green(pt), Blue(pt));
		//===============================================================
		//Convert Canvas space to unit space, where 0,0 is in middle of Isosceles triange
		//up id po
		private double deltaY(double GreenY, double PY) => (PY - GreenY) / ht;
		//==============================================================
		public Point RelativePoint(Point pt)/*place 0,0 in middle of rectangle*/ {
			Rect rb = BoundingRect;

			double x = -0.5 + (pt.X - rb.Left) / side;
			double y = (rb.Height / 2 - pt.Y - rb.Top) / side;
			return new Point(x, y);
		}
		//===============================================================
		public Color IsoscelesColor(Point pt) {
			if (!IsInside(pt)) return Colors.Black;
			Color c = new Color();
			c.A = 255;
			c.R = Red(pt);
			c.G = Green(pt);
			c.B = Blue(pt);
			return c;
		}
		//=========================================================================
		public string[] Details() {
			string[] sa = new string[4];
			sa[0] = '\n' + Name + ", origin:" + CSNUtil.sPoint(cp) + ", side = " + side.ToString("f2") + ", ht = " + ht.ToString("f2");
			sa[1] = sDetail("\nA0:", A[0], this);
			sa[2] = sDetail("\nA1:", A[1], this);
			sa[3] = sDetail("\nA2:", A[2], this);
			return sa;
		}
		//...........................................................................
		private static string sDetail(string sTitle, Point p, CIsosceles T) {
			string st = sTitle + CSNUtil.sPoint(p) + " => " + CSNUtil.sPoint(T.RelativePoint(p))
			+ ", color:" + CSNUtil.sColor(T.IsoscelesColor(p))
			+ ", 120" + CSNUtil.sDegree + ": " + CSNUtil.sPoint(RotatePoint(p, T.Center, -0.5, 0.866025404))
			+ ", 240" + CSNUtil.sDegree + ": " + CSNUtil.sPoint(RotatePoint(p, T.Center, -0.5, -0.866025404));
			return st;
		}
		//=========================================================================
		public override string ToString() => Name + ":["
		+ " cp:" + CSNUtil.sPoint(cp)
		+ " - R:" + CSNUtil.sPoint(A[0])
		+ ", G:" + CSNUtil.sPoint(A[1])
		+ ", B:" + CSNUtil.sPoint(A[2]) + ']'
		+ ", side:" + side.ToString("f2")
		+ ", ht:" + ht.ToString("f2");

	}
}
