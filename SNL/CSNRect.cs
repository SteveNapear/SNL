using System;

namespace SNL {
	//////////////////////////////////////////////////////////////////////
	[Serializable]
	public class CSNRect    //Complex plane rectangle.  Note: top is > than bottom
	{
		private static string sPrec = "F6";
		private CSNplex ll;
		private double wdt;
		private double hgt;
		//================================================================
		public CSNRect() { ll = new CSNplex(0, 0); wdt = 0; hgt = 0; }
		public CSNRect(CSNplex LL) { ll = LL; wdt = 0; hgt = 0; }
		public CSNRect(CSNplex LL, double W, double H) { ll = LL; wdt = W; hgt = H; }
		public CSNRect(double left, double bottom, double W, double H) { ll = new CSNplex(left, bottom); wdt = W; hgt = H; }
		//================================================================
		public double Width { get { return wdt; } set { wdt = value; } }
		public double Height { get { return hgt; } set { hgt = value; } }
		public double Left { get { return ll.r; } }
		public double Top { get { return ll.i + hgt; } }
		public double Right { get { return ll.r + wdt; } }
		public double Bottom { get { return ll.i; } }
		public CSNplex LowerLeft { get { return ll; } set { ll = value; } }
		public CSNplex TopRight { get { return new CSNplex(ll.r + wdt, ll.i + hgt); } }
		public double Area { get { return wdt * hgt; } }
		//===============================================================
		public override int GetHashCode() => base.GetHashCode();
		//================================================================
		public override string ToString()
		{
			return "[(" + ll.r.ToString(sPrec) + ',' + ll.i.ToString(sPrec) + ')' +
			',' + wdt.ToString(sPrec) + ',' + hgt.ToString(sPrec) + ']';
		}
		public CSNplex Center { get { return new CSNplex(ll.r + wdt / 2.0, ll.i + hgt / 2.0); } }
		//................................................................
		public static bool operator ==(CSNRect lhs, CSNRect rhs)    //comparison operator
		{
			if (lhs.Left != rhs.Left) return false;
			if (lhs.Bottom != rhs.Bottom) return false;
			if (lhs.Width != rhs.Width) return false;
			if (lhs.Height != rhs.Height) return false;
			return true;
		}
		public static bool operator !=(CSNRect lhs, CSNRect rhs) { return !(lhs == rhs); }
		public override bool Equals(Object obj)
		{
			if (!(obj is CSNRect)) return false;
			return this == (CSNRect)obj;
		}
	}
}