//////////////////////////////////////////////////////////////////////////
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;

namespace SNL {
	///////////////////////////////////////////////////////////////////////
	[Serializable]
	public class DrawingCanvas : Canvas /*from 2010 book*/{
		private static int count = 0;
		private static string sFont = "Times New Roman";
		private List<Visual> visuals = new List<Visual>();
		protected override int VisualChildrenCount => visuals.Count;
		protected override Visual GetVisualChild(int index) {
			return index < visuals.Count ? visuals[index] : null;
		}

		public int Count => count;
		//==================================================================
		public void AddVisual(Visual v) /*adds visual to list*/	{
			visuals.Add(v);
			base.AddVisualChild(v);
			base.AddLogicalChild(v);
			count += 1;
		}
		//==================================================================
		public void DeleteVisual(Visual v) /*removes visual from list*/	{
			visuals.Remove(v);
			base.RemoveLogicalChild(v);
			base.RemoveVisualChild(v);
			count -= 1;
		}
		//==================================================================
		public DrawingVisual GetVisual(Point pt)   /*retrieve the visial for near this point*/	{
			HitTestResult tr = VisualTreeHelper.HitTest(this, pt);
			return tr.VisualHit as DrawingVisual;
		}
		//==================================================================
		public void ClearAllVisuals()   /*clears out the visuals for this element*/		{
			while (VisualChildrenCount > 0) {
				Visual v = visuals[0];
				DeleteVisual(v);
			}
		}
		//..............................................................
		public static FormattedText DrawTextInRect(
			DrawingContext dc, Rect rText, string sText, double emPC,
			string fontfamily, FontStyle fontstyle, FontWeight FW,
			Brush br, TextAlignment ta, double PixelsPerDip) {

			if (rText.IsEmpty) return null;
			double emsize = rText.Height * 0.5;
			if (rText.Width * rText.Height == 0) return null;

			double pc, wwd, delta, twd = rText.Width, epsilon = 0.2, wd, low, high, oldem, em = rText.Height * emPC;
			if (em < 0.05) return null;
			if (em > rText.Height) em = rText.Height;

			FlowDirection flow = FlowDirection.LeftToRight;
			Typeface face = new Typeface(fontfamily);
			CultureInfo CI = new CultureInfo("en-us");
			FormattedText FT = new FormattedText(sText, CI, flow, face, em, br, PixelsPerDip);

			FT.SetFontWeight(FontWeights.Bold);
			FT.SetFontStyle(fontstyle);
			FT.TextAlignment = ta;
			FT.MaxTextWidth = rText.Width;
			FT.MaxTextHeight = rText.Height;
			FT.SetFontWeight(FW);
			FT.Trimming = TextTrimming.CharacterEllipsis;

			int ni = 0;
			wwd = rText.Width * emPC;
			low = 0; high = rText.Height;
			do {
				oldem = em;
				em = low + (high - low) / 2.0;
				ni += 1;
				FT.SetFontSize(em);
				wd = FT.WidthIncludingTrailingWhitespace;
				if (wd < 2) { FT.SetFontSize(oldem); break; }
				delta = wwd - wd;
				pc = delta / wwd;

				if (wd >= wwd)/*too big*/ { high = em; continue; }
				//is it within range, wd is < rText.Width ==> delta is always positive from here
				if (pc < epsilon) break;
				low = em;
			} while (true);

			Point ul = new Point(rText.Left, rText.Top);
			delta = rText.Height - FT.Height;
			Point pt = new Point(rText.Left, rText.Top + delta / 2);
			dc.DrawText(FT, pt);
			return FT;
		}

		//..................................................................
		public static void DrawTextAtPoint(
			DrawingContext dc, Point pt, string sText, double emSize,
			string fontfamily, FontStyle fontstyle, FontWeight FW,
			Brush br, TextAlignment ta, double PixelsPerDip) {

			FlowDirection flow = FlowDirection.LeftToRight;
			Typeface face = new Typeface(fontfamily);
			CultureInfo CI = new CultureInfo("en-us");
			FormattedText FT = new FormattedText(sText, CI, flow, face, emSize, br, PixelsPerDip);

			FT.TextAlignment = ta;
			FT.SetFontWeight(FW);
			FT.SetFontStyle(fontstyle);
			FT.Trimming = TextTrimming.CharacterEllipsis;
			dc.DrawText(FT, new System.Windows.Point(pt.X, pt.Y));
		}
		//...................................................................
		public static void DrawCross(DrawingContext dc, Rect R, Brush BR)    //draw cross lines
		{
			Pen pen = new Pen(BR, 1.0);
			dc.DrawLine(pen, new Point(0, 0), new Point(R.Right, R.Bottom));
			dc.DrawLine(pen, new Point(0, R.Bottom), new Point(R.Right, 0));
		}
		//................................................................
		public static void DrawCross(DrawingContext dc,Point pt,Brush br){
			const double delta = 4.0;
			Pen pen = new Pen(br,1.0);
			dc.DrawLine(pen,new Point(pt.X-delta,pt.Y),new Point(pt.X+delta,pt.Y));
			dc.DrawLine(pen,new Point(pt.X,pt.Y-delta),new Point(pt.X,pt.Y+delta));
		}
		//.....................................................................
		public static void DrawMark(DrawingContext dc, Point P, Pen pen)    //draws a mark
		{
			const double delta = 4.0;
			Point p1 = new Point(P.X - delta, P.Y), p2 = new Point(P.X + delta, P.Y);
			dc.DrawLine(pen, p1, p2);
			dc.DrawLine(pen, new Point(P.X, P.Y - delta), new Point(P.X, P.Y + delta));
		}
	}//end class
}

