//************************************************************************
//This package holds classes that can be used for Geo navigation
//This includes classes for sailboat racing
//Copyrgiht(C) Equin Technology, Inc. 2011
//************************************************************************

using System;
using System.IO;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Windows.Forms;
using System.Runtime.Serialization;
using System.Runtime.Serialization.Formatters.Binary;


namespace SNL
{
	//////////////////////////////////////////////////////////////////////
	[Serializable]
	public class CGeo
	{   //A geo point
		//all angles in radians, all distances in kilometers
		public static readonly double d2r = Math.PI / 180.0;
		public static readonly double r2d = 180.0 / Math.PI;
		public static readonly double km2nm = 1.0 / 1.852;  //kilometers to nautical miles (1.852 km/nm def)
		public static readonly double km2mi = 0.621504;     //kilometers to statue miles
		public static readonly double cir = 40030.173592;   //circumference based on mean radius in km
		public static readonly double radius = 6371.0;  //mean radius km
														//..................................................................
		private static char[] dms = { (char)176, (char)39, (char)34 };

		//..................................................................
		private double lat;
		private double lon;
		//================================================================
		public CGeo(double DegLat, double DegLon) { lat = DegLat; lon = DegLon; }
		public CGeo() { lat = 0.0; lon = 0.0; }
		public CGeo(CGeo geo) { lat = geo.lat; lon = geo.lon; }
		public double Lat { get { return lat; } set { lat = value; } }
		public double Lon { get { return lon; } set { lon = value; } }
		public double LatRad { get { return lat * d2r; } }
		public double LonRad { get { return lon * d2r; } }
		public double EarthCircumferenceKm { get { return cir; } }
		public double EarthRadiusKm { get { return radius; } }
		//................................................................
		public static double Deg2Rad(double deg) { return deg * d2r; }
		public static double Rad2Deg(double rad) { return rad * r2d; }
		//................................................................
		public static double[] LatLon2XYZ(CGeo gpt) { return LatLon2XYZ(gpt.Lat, gpt.Lon); }
		public static double[] LatLon2XYZ(double Lat, double Lon)
		{
			double z = Math.Sin(d2r * Lat);
			double t = Math.Cos(d2r * Lat);
			double x = Math.Sin(d2r * Lon) * t;
			double y = Math.Cos(d2r * Lon) * t;
			double[] D = new double[3];
			D[0] = x; D[1] = y; D[2] = z;
			return D;
		}
		//====================================================================
		public bool Equals(CGeo geo) { if ((lat != geo.lat) || (lon != geo.lon)) return false; return true; }
		//====================================================================
		public double[] XYZ { get { return LatLon2XYZ(this); } }
		//====================================================================
		public double DotDist(CGeo geo) { return DotDist(this, geo); }
		//====================================================================
		public double GCArcAngleRads(CGeo geo)
		{   //returns angle in radians of great circle arc
			double ACOS = Math.Acos(this * geo);
			return ACOS;
		}
		//====================================================================
		public double Heading(CGeo to) { return Heading(this, to); }
		public double HeadingDeg(CGeo to) { return r2d * Heading(this, to); }
		public string sDegLat
		{
			get
			{
				char c; double v;
				if (lat < 0.0) { c = 'S'; v = -lat; } else { c = 'N'; v = lat; }
				return c + v.ToString("f5") + dms[0];
			}
		}
		public string sDegLon
		{
			get
			{
				char c; double v;
				if (lon < 0.0) { c = 'W'; v = -lon; } else { c = 'E'; v = lon; }
				return c + v.ToString("f5") + dms[0];
			}
		}
		public string sDeg { get { return sDegLat + ", " + sDegLon; } }
		public string sDMM { get { return CGeo.DMM(true, lat) + ", " + CGeo.DMM(false, lon); } }
		public string sDMS { get { return CGeo.DMS(true, lat) + ", " + CGeo.DMS(false, lon); } }
		//====================================================================
		public override string ToString()
		{
			return '[' + lat.ToString("f5") +
			 dms[0] + ',' + lon.ToString("f5") +
			 dms[0] + ']';
		}
		//................................................................
		public static double Heading(CGeo fr, CGeo to)
		{
			double dLon = CGeo.d2r * (to.Lon - fr.Lon);
			double cosLat2 = Math.Cos(CGeo.d2r * to.Lat);
			double t1 = Math.Sin(dLon) * cosLat2;
			double x1 = Math.Cos(CGeo.d2r * fr.Lat) * Math.Sin(CGeo.d2r * to.Lat);
			double x2 = Math.Sin(CGeo.d2r * fr.Lat) * cosLat2 * Math.Cos(dLon);
			double angle = Math.Atan2(t1, x1 - x2);
			if (angle < 0.0) angle += (2.0 * Math.PI);
			return angle;
		}
		//.....................................................................
		public static double operator *(CGeo lhs, CGeo rhs)
		{       //dot product
			double[] lft = LatLon2XYZ(lhs);
			double[] rgt = LatLon2XYZ(rhs);
			double dot = lft[0] * rgt[0] + lft[1] * rgt[1] + lft[2] * rgt[2];
			return dot;
		}
		//....................................................................
		public static double operator -(CGeo lhs, CGeo rhs) { return WGS84km(lhs, rhs); }
		//....................................................................
		public static double DotDist(CGeo lhs, CGeo rhs)
		{
			double dot = lhs * rhs;
			double ACOS = Math.Acos(dot);
			double acos = lhs.GCArcAngleRads(rhs);      //debug
			return cir * ACOS / (2.0 * Math.PI);

		}
		//used for tuning against WGS84 algorithm
		public static double DotDistRad(CGeo lhs,CGeo rhs,double circum){
			double dot = lhs * rhs;
			double ACOS = Math.Acos(dot);
			return circum * ACOS / (2.0 * Math.PI);
		}
		//....................................................................
		public static string DMM(bool bLat, double degrees)
		{
			int i, j, ideg, imin; float mn, sec; StringBuilder sb = new StringBuilder(15);

			bool b = (degrees < 0) ? true : false;
			if (b) degrees = -degrees;
			if (!Deg2DMS(degrees, out ideg, out imin, out sec)) return "??";

			char[] ca = new char[6];
			i = 0;
			ca[i++] = (bLat ? (b ? 'S' : 'N') : (b ? 'W' : 'E'));
			if (!bLat) ca[i++] = (char)('0' + ideg / 100);
			i2s2(ideg, ca, i); i += 2;
			ca[i] = dms[0];
			for (j = 0; j < i; j++) sb.Append(ca[j]);
			sb.Append(dms[0]);

			mn = imin + sec / 60.0f;
			if (mn < 10) sb.Append('0');
			sb.Append(mn.ToString("f3"));
			sb.Append(dms[1]);
			return sb.ToString();
		}
		//....................................................................
		private static void i2s2(int v, char[] ca, int ix) { ca[ix + 1] = (char)('0' + v % 10); v /= 10; ca[ix] = (char)('0' + v % 10); }
		//....................................................................
		public static string DMS(bool bLat, double degrees)
		{
			int i, j, ideg, imin; float sec; StringBuilder sb = new StringBuilder(15);

			bool b = (degrees < 0) ? true : false;
			if (b) degrees = -degrees;
			if (!Deg2DMS(degrees, out ideg, out imin, out sec)) return "??";
			char[] ca = new char[6];
			//do degrees
			i = 0;
			ca[i++] = (bLat ? (b ? 'S' : 'N') : (b ? 'W' : 'E'));
			if (!bLat) ca[i++] = (char)('0' + ideg / 100);
			i2s2(ideg, ca, i); i += 2;
			ca[i] = dms[0];
			for (j = 0; j <= i; j++) sb.Append(ca[j]);

			//do minutes
			i = 0; //ca[i++] = ' ';
			i2s2(imin, ca, i); i += 2;
			ca[i++] = dms[1];
			for (j = 0; j < i; j++) sb.Append(ca[j]);

			//do seconds
			i = 0; //ca[i++] = ' ';
			i2s2((int)sec, ca, i); i += 2;
			ca[i++] = dms[2];
			for (j = 0; j < i; j++) sb.Append(ca[j]);
			return sb.ToString();
		}
		//....................................................................
		public static bool Deg2DMS(float fdeg, out int ideg, out int min, out float sec)
		{
			return Deg2DMS((double)fdeg, out ideg, out min, out sec);
		}
		//....................................................................
		public static bool Deg2DMS(double deg, out int ideg, out int min, out float sec)
		{
			double d = (deg < 0.0) ? -deg : deg; //make postive
			ideg = min = 0; long s10, mn, rem; sec = 0.0f;

			s10 = (long)(deg * 36000.0);
			s10 = (long)Math.Round(deg * 36000.0, 1);
			mn = Math.DivRem(s10, 600, out rem);
			sec = (float)rem / 10.0f;
			ideg = (int)Math.DivRem(mn, 60, out rem);
			min = (int)rem;
			return true;
		}
		//====================================================================

		public double WGS84km(CGeo to) { return WGS84km(this, to); }
		public double WGS84nm(CGeo to) { return km2nm * WGS84km(this, to); }
		public double WGS84mi(CGeo to) { return km2mi * WGS84km(this, to); }
		//....................................................................
		private const double minoraxis = 6356752.314245;//earth minor axis
		private const double majoraxis = 6378137;       //earth major axis
		private const double fratio = 1.0 / 298.257223563;
		private const double small = 1e-12;
		//....................................................................
		public static double WGS84km(CGeo fr, CGeo to)
		{   //Vincenty's formulae 
			if (fr.Equals(to)) return 0.0;
			double a = 6378137.0;
			double b = 6356752.3142;
			double f = 1 / 298.257223563;

			double l = fr.LonRad - to.LonRad;

			double u1 = Math.Atan((1 - f) * Math.Tan(fr.LatRad));
			double u2 = Math.Atan((1 - f) * Math.Tan(to.LatRad));
			double sin_u1 = Math.Sin(u1);
			double cos_u1 = Math.Cos(u1);
			double sin_u2 = Math.Sin(u2);
			double cos_u2 = Math.Cos(u2);

			double lambda = l;
			double lambda_pi = 2 * Math.PI;
			int iter_limit = 20;

			double cos_sq_alpha = 0.0;
			double sin_sigma = 0.0;
			double cos2sigma_m = 0.0;
			double cos_sigma = 0.0;
			double sigma = 0.0;

			while (Math.Abs(lambda - lambda_pi) > 1e-12 && --iter_limit > 0)
			{
				double sin_lambda = Math.Sin(lambda);
				double cos_lambda = Math.Cos(lambda);

				sin_sigma = Math.Sqrt((cos_u2 * sin_lambda) * (cos_u2 * sin_lambda) +
					(cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda) *
					(cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda));

				cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lambda;
				sigma = Math.Atan2(sin_sigma, cos_sigma);

				double alpha = Math.Asin(cos_u1 * cos_u2 * sin_lambda / sin_sigma);
				var t1 = Math.Cos(alpha);
				cos_sq_alpha = t1 * t1;
				cos2sigma_m = cos_sigma - 2.0 * sin_u1 * sin_u2 / cos_sq_alpha;

				double cc = f / 16.0 * cos_sq_alpha * (4.0 + f * (4.0 - 3.0 * cos_sq_alpha));
				lambda_pi = lambda;
				lambda = l + (1.0 - cc) * f * Math.Sin(alpha) *
					(sigma + cc * sin_sigma * (cos2sigma_m + cc * cos_sigma * (-1.0 + 2.0 * cos2sigma_m * cos2sigma_m)));
			}

			double usq = cos_sq_alpha * (a * a - b * b) / (b * b);
			double aa = 1.0 + usq / 16384.0 * (4096.0 + usq * (-768.0 + usq * (320.0 - 175.0 * usq)));
			double bb = usq / 1024.0 * (256.0 + usq * (-128.0 + usq * (74.0 - 47.0 * usq)));
			double delta_sigma = bb * sin_sigma * (cos2sigma_m + bb / 4.0 * (cos_sigma * (-1.0 + 2.0 * cos2sigma_m * cos2sigma_m) -
				bb / 6.0 * cos2sigma_m * (-3.0 + 4.0 * sin_sigma * sin_sigma) * (-3.0 + 4.0 * cos2sigma_m * cos2sigma_m)));
			double c = b * aa * (sigma - delta_sigma);

			return c / 1000.0;  //return km
		}
		public CGeo PTBD2Geo(double distkm, double beartruedeg) { return PTBD2Geo(this, distkm, beartruedeg); }
		//=======================================================================
		public static CGeo PTBD2Geo(CGeo fr, double distkm, double beartruedeg)
		{
			if (distkm == 0.0) return new CGeo(fr); //no distance, same point!

			double a = 6378137.0;
			double b = 6356752.3142;
			double f = 1 / 298.257223563;
			double s = distkm * 1000.0;             //km-->meters
			double alpha1 = beartruedeg * CGeo.d2r; //var alpha1 = brng.toRad();
			double sinAlpha1 = Math.Sin(alpha1);    //var sinAlpha1 = Math.sin(alpha1);
			double cosAlpha1 = Math.Cos(alpha1);    //var cosAlpha1 = Math.cos(alpha1);
			double tanU1 = (1 - f) * Math.Tan(fr.LatRad);//var tanU1 = (1-f) * Math.tan(lat1.toRad());
			double cosU1 = 1 / Math.Sqrt((1 + tanU1 * tanU1));//var cosU1 = 1 / Math.sqrt((1 + tanU1*tanU1))
			double sinU1 = tanU1 * cosU1;               //, sinU1 = tanU1*cosU1;
			double sigma1 = Math.Atan2(tanU1, cosAlpha1);   //var sigma1 = Math.atan2(tanU1, cosAlpha1);
			double sinAlpha = cosU1 * sinAlpha1;        //var sinAlpha = cosU1 * sinAlpha1;
			double cosSqAlpha = 1 - sinAlpha * sinAlpha;    //var cosSqAlpha = 1 - sinAlpha*sinAlpha;  
			double uSq = cosSqAlpha * (a * a - b * b) / (b * b);    //var uSq = cosSqAlpha * (a*a - b*b) / (b*b);
			double A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq))); //var A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
			double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));//var B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));
			double sigma = s / (b * A);             //var sigma = s / (b*A), 
			double sigmaP = 2 * Math.PI;                //sigmaP = 2*Math.PI;
			double sinSigma = 0, cosSigma = 0, cos2SigmaM = 0;
			while (Math.Abs(sigma - sigmaP) > 1e-12)
			{//while (Math.abs(sigma-sigmaP) > 1e-12) {    
				cos2SigmaM = Math.Cos(2 * sigma1 + sigma);  //var cos2SigmaM = Math.cos(2*sigma1 + sigma);    
				sinSigma = Math.Sin(sigma); //var sinSigma = Math.sin(sigma);    
				cosSigma = Math.Cos(sigma); //var cosSigma = Math.cos(sigma); 
				double deltaSigma = B * sinSigma * (cos2SigmaM + B / 4 *
									(cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) -
									B / 6 * cos2SigmaM * (-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
				sigmaP = sigma;
				sigma = s / (b * A) + deltaSigma;
			}
			double tmp = sinU1 * sinSigma - cosU1 * cosSigma * cosAlpha1; //var tmp = sinU1*sinSigma - cosU1*cosSigma*cosAlpha1;
			double lat2 = Math.Atan2(           //var lat2 = Math.atan2(sinU1*cosSigma + 
						sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1,
						(1 - f) * Math.Sqrt(sinAlpha * sinAlpha + tmp * tmp)//	cosU1*sinSigma*cosAlpha1, (1-f)*Math.sqrt(sinAlpha*sinAlpha + tmp*tmp));  
					);
			double lambda = Math.Atan2(sinSigma * sinAlpha1, cosU1 * cosSigma - sinU1 * sinSigma * cosAlpha1);
			double C = f / 16 * cosSqAlpha * (4 + f * (4 - 3 * cosSqAlpha));
			double L = lambda - (1 - C) * f * sinAlpha *
					(sigma + C * sinSigma * (cos2SigmaM + C * cosSigma *
					(-1 + 2 * cos2SigmaM * cos2SigmaM)));

			double lon2 = (fr.LonRad + L + 3 * Math.PI) % (2 * Math.PI) - Math.PI;  // normalise to -180...+180  
			double revAz = Math.Atan2(sinAlpha, -tmp);  // final bearing, if required 
			CGeo geo = new CGeo(lat2 * CGeo.r2d, lon2 * CGeo.r2d);
			return geo;
		}
		//................................................................
		public static double Lat2Mercator(double Latdeg)	{
			const double π = Math.PI, pd4 = π / 4.0, d2r = π / 180.0;
			const double maxdeg = 85.0;

			if (Latdeg > maxdeg || Latdeg < -maxdeg) return double.NaN;
			double lv, d = 0.0; bool N = (Latdeg < 0 ? false : true);
			double θ = (N ? Latdeg : -Latdeg);
			double rad = d2r * θ / 2.0;
			double angle = pd4 + rad;
			d = Math.Tan(angle);
			lv = Math.Log(d);
			if (!N) lv = -lv;
			return lv;
		}
		//................................................................
		public static double Mercator2Lat(double v)		{
			const double π = Math.PI, pd4 = π / 4.0, d2r = π / 180.0;
			bool bn = (v < 0);
			if (bn) v = -v;
			double t1 = Math.Exp(v);
			double t2 = Math.Atan(t1);
			double t3 = t2 - pd4;
			Double t4 = 2.0*t3 / d2r;
			return (bn ? -t4 : t4);
		}
		//................................................................
		public static double Lat2GallPeters(double deglat)		{
			const double maxdeg = 85.0;
			if (deglat > maxdeg || deglat < -maxdeg) return double.NaN;
			double ans = Math.Sin(Math.PI * deglat / 180.0) * Math.Sqrt(2.0);
			return ans;
		}
		//................................................................
		public static double GallPeters2Lat(double v)
		{
			double t = Math.Asin(v / Math.Sqrt(2.0));
			return 180.0 * t / Math.PI;
		}
	}

	//////////////////////////////////////////////////////////////////////
	[Serializable]
	public class CWayPoint : CGeo
	{   //A CGeo with a name and a flood tide heading (true)
		private string sname;
		private double floodhdg;                    //true heading for flood tide at this location
													//================================================================
		private void Init() { sname = ""; floodhdg = 0.0; }
		//================================================================
		public CWayPoint() : base(0, 0) { Init(); }
		public CWayPoint(string sName) : base(0, 0) { sname = sName; }
		public CWayPoint(string sName, double Lat, double Lon) : base(Lat, Lon) { sname = sName; }
		public CWayPoint(CWayPoint wp) : base(wp.Lat, wp.Lon) { sname = wp.sWP; }
		public CWayPoint(string sName, double Lat, double Lon, double FloodHdg)
			: base(Lat, Lon)
		{
			sname = sName; floodhdg = FloodHdg;
		}
		//................................................................
		//===============================================================
		public string sWP { get { return sname; } set { sname = value; } }
		public double FloodHdgDeg { get { return floodhdg; } set { floodhdg = value; } }
		public double FloodHdgRad { get { return floodhdg * d2r; } }

		//================================================================
		public void CopyFrom(CWayPoint fr) { sWP = fr.sWP; Lat = fr.Lat; Lon = fr.Lon; }
		//=================================================================
		public double HeadingTo(CWayPoint wp)
		{       //returns angle, in Radians, from true north, to wp.
			double[] here = LatLon2XYZ(this);
			double[] to = LatLon2XYZ(wp);
			return 0.0;
		}
		//=================================================================
		//CBG1  32N 43’ 09" 43.150 117W 13’ 02" 13.033
		//=================================================================
		public static bool TryParse(string sline, out CWayPoint wp)
		{
			string sname, st, sm; double lat = 0, lon = 0, d; int ix; char c;

			wp = null;
			sline.Trim();
			string[] sa = sline.Split(' ');
			if (sa.Length < 8) return false;
			sname = sa[0];
			st = sa[1];
			ix = 0; while (char.IsDigit(st[ix])) ix += 1;
			c = (ix < st.Length) ? st[ix] : '\0';
			sm = st.Substring(0, ix);
			if (!double.TryParse(sm, out lat)) return false;
			if (!double.TryParse(sa[4], out d)) return false;
			lat += (d / 100.0);
			if (c == 'S') lat = -lat;
			st = sa[5];
			ix = 0; while (char.IsDigit(st[ix])) ix += 1;
			c = (ix < st.Length) ? st[ix] : '\0';
			sm = st.Substring(0, ix);
			if (!double.TryParse(sm, out lon)) return false;
			if (!double.TryParse(sa[8], out d)) return false;
			lon += (d / 100.0);
			if (c == 'W') lon = -lon;
			wp = new CWayPoint(sname, lat, lon);
			return true;
		}
		//================================================================
		public string CSV
		{
			get
			{
				return '\"' + sname + "\"," +
					Lat.ToString("f5") + ',' + Lon.ToString("f5") + ',' + FloodHdgDeg.ToString("F3");
			}
		}
		//................................................................
		public static CWayPoint ReadCSV(string s)
		{
			//Reads the CSV generated line above: sa[0]:<name>,sa[1]:<lat>,sa[2]:<lon>,sa[3]:<hdg>
			CWayPoint wp = null; double lat, lon, hdg;

			s.Trim();
			string[] sa = CDASlib.CSVSplit(s, ',');
			if (sa.Length < 3) return null;
			if (!double.TryParse(sa[1], out lat)) return null;
			if (!double.TryParse(sa[2], out lon)) return null;
			if (sa.Length < 4) hdg = 0;
			else { if (!double.TryParse(sa[3], out hdg)) return null; }
			wp = new CWayPoint(sa[0], lat, lon, hdg);
			return wp;
		}
		//================================================================
		public override string ToString() { return sname + sLatLon; }
		public string sLatLon { get { return " [" + Lat.ToString("f5") + ',' + Lon.ToString("f5") + ']'; } }
	}
	//////////////////////////////////////////////////////////////////////
	[Serializable]
	public class CWPList
	{
		private string spath = "";
		private List<CWayPoint> list = null;
		private DateTime DOC;
		private DateTime DOLU;
		private bool bsaved = false;
		//................................................................
		private int ByName(CWayPoint s1, CWayPoint s2) { return s1.sWP.CompareTo(s2.sWP); }
		private int ByLat(CWayPoint s1, CWayPoint s2)
		{
			if (s1.Lat > s2.Lat) return 1;
			if (s1.Lat < s2.Lat) return -1;
			return 0;
		}
		private int ByLon(CWayPoint s1, CWayPoint s2)
		{
			if (s1.Lon > s2.Lon) return 1;
			if (s1.Lon < s2.Lon) return -1;
			return 0;
		}

		//================================================================
		private void Init(uint nwps)
		{
			list = new List<CWayPoint>((int)nwps);
			DOC = DateTime.Now;
			DOLU = DOC;
		}
		public CWPList() { Init(100); }
		public CWPList(uint nWPs) { Init(nWPs); }
		public void Add(CWayPoint wp) { list.Add(wp); DOLU = DateTime.Now; }
		public int nWP { get { return list.Count; } }
		public CWayPoint WP(int ix) { return (ix < list.Count) ? list[ix] : null; }
		public bool isSaved { get { return bsaved; } }
		public string sPath { get { return spath; } }
		public void SortByName() { list.Sort(ByName); }
		public void SortByLat() { list.Sort(ByLat); }
		public void SortByLon() { list.Sort(ByLon); }
		private string sFindWP = "";
		//================================================================
		public bool Delete(int ix)
		{
			if (ix >= nWP) return false;
			list.RemoveAt(ix);
			return true;
		}
		//================================================================
		public int FindIndex(string sWP)
		{
			sFindWP = sWP;
			return list.FindIndex(delegate (CWayPoint WP) { return (WP.sWP == sFindWP); });
		}
		//================================================================
		public CWayPoint Find(string sWP)
		{
			sFindWP = sWP;
			return list.Find(delegate (CWayPoint WP) { return (WP.sWP == sFindWP); });
		}
		//................................................................
		public static bool SaveWPList(string fns, CWPList list)
		{
			//Opens a file and serializes the object into it in binary format.
			string se; Stream stream = null;

			list.bsaved = false;
			try
			{
				if (fns == "") throw new Exception("Invalid file name");
				if (!Directory.Exists(Path.GetDirectoryName(fns))) throw new Exception("Invalid directory");

				stream = File.Open(fns, FileMode.Create, FileAccess.Write, FileShare.None);
				BinaryFormatter formatter = new BinaryFormatter();
				formatter.AssemblyFormat = System.Runtime.Serialization.Formatters.FormatterAssemblyStyle.Simple;
				formatter.Serialize(stream, list);      //save file 
			}
			catch (Exception exp)
			{
				se = "An error has occurred while trying to write to file:\n\"" + fns + '\"' +
					"\nError: " + exp.Message;
				MessageBox.Show(se);
				if (stream != null) stream.Close();
				return false;
			}
			if (stream != null) stream.Close();
			list.spath = fns;
			list.bsaved = true;
			return true;
		}
		//................................................................
		public static CWPList ReadWPList(string fns)
		{
			if (!File.Exists(fns)) return null;
			//Opens a file a deserializes it into the object in binary format
			CWPList list = null; string se;
			Stream stream = File.Open(fns, FileMode.Open, FileAccess.Read, FileShare.Read);
			BinaryFormatter formatter = new BinaryFormatter();
			formatter.AssemblyFormat = System.Runtime.Serialization.Formatters.FormatterAssemblyStyle.Simple;
			try
			{
				list = (CWPList)formatter.Deserialize(stream);
			}
			catch (Exception exp)
			{
				se = "Invalid CWPList File:\n" + fns +
					"\nIO Error: " + exp.Message;
				list = null;
			}
			finally
			{
				stream.Close();
			}
			if (list == null)
			{
				se = "Couldn't retrieve Waypoint list file:\n" + fns +
				"\nreturning empty.";
				MessageBox.Show(se, "***Error***", MessageBoxButtons.OK, MessageBoxIcon.Error);
				return null;
			}
			list.spath = fns;
			list.bsaved = true;
			return list;
		}
	};
	//////////////////////////////////////////////////////////////////////
	[Serializable]
	public class CMark
	{
		private string smark;                   //name of mark
		private string swp;                     //name of waypoint in waypoint list
		private bool bport = true;              //true if mark is to be taken to port
												//====================================================================
		public CMark(string sN, string wpname, bool bToPort) { smark = sN; swp = wpname; bport = bToPort; }
		public CMark(string sN, CWayPoint wp, bool bToPort) { smark = sN; swp = wp.sWP; bport = bToPort; }
		//====================================================================
		public string sMark { get { return smark; } }
		public string sWP { get { return swp; } }
		public bool isPort { get { return bport; } set { bport = value; } }
		public override string ToString() { return smark + '[' + swp + ']' + '(' + (bport ? 'P' : 'S') + ')'; }

	}
	//////////////////////////////////////////////////////////////////////
	[Serializable]
	public class CRoute
	{
		private const int nm = 10;
		private List<CMark> marks = null;       //list of marks
		private string sname;
		private DateTime DOC;
		private DateTime DOLU;
		//================================================================
		private void Init(int nmarks, string st)
		{
			if (nmarks < 10) nmarks = 10;
			marks = new List<CMark>(nmarks);
			sname = st;
			DOC = DateTime.Now;
			DOLU = DOC;
		}
		//================================================================
		public CRoute(string sN) { Init(nm, sN); }
		public int nMarks { get { return marks.Count; } }
		public DateTime Created { get { return DOC; } }
		public DateTime Updated { get { return DOLU; } }
		public string sRoute { get { return sname; } set { sname = value; } }
		public void Add(CMark mark) { marks.Add(mark); }
		public void Add(string smark, CWayPoint wp, bool bToPort)
		{
			string sm = "Mark " + marks.Count.ToString("D");
			if ((smark == null) || (smark.Length == 0)) smark = sm;
			Add(new CMark(smark, wp.sWP, bToPort));
		}
		public CMark Mark(int ix) { return (ix < marks.Count) ? marks[ix] : null; }
		public int FindIndex(string sMark) { return marks.FindIndex(delegate (CMark mark) { return (mark.sMark == sMark); }); }
		public bool DeleteMark(string sMark)
		{
			int ix = FindIndex(sMark);
			if ((ix < 0) || (ix >= marks.Count)) return false;
			marks.RemoveAt(ix);
			return true;
		}
		public bool DeleteMark(int ix)
		{
			if ((ix < 0) || (ix >= marks.Count)) return false;
			marks.RemoveAt(ix);
			return true;
		}
		//====================================================================
		public CMark FromWP(CWayPoint wp)
		{
			for (int ix = 0; ix < marks.Count; ix++) { if (marks[ix].sWP.Equals(wp.sWP)) return marks[ix]; }
			return null;
		}
	};
}   //end namespace