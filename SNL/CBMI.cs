using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace SNL {
	public static class CBMI {
		//................................................................
		public static byte[] BMI2BA(BitmapSource bmi, int row, int col) {   //assume row 0 at bottom,col0 at left 
			Size isz = new System.Windows.Size(bmi.PixelWidth, bmi.PixelHeight);
			int stride = (int)bmi.PixelWidth * (bmi.Format.BitsPerPixel + 7 / 8);
			byte[] pixels = new byte[(int)bmi.PixelHeight * stride];
			bmi.CopyPixels(pixels, stride, 0);

			return pixels;
		}
		//.................................................................
		public static WriteableBitmap BA2BMI(byte[] pixels, int width, int height) {
			WriteableBitmap bmi = new WriteableBitmap(width, height, 360, 360, PixelFormats.Bgra32, null);
			int stride = width * (bmi.Format.BitsPerPixel + 7 / 8);
			bmi.WritePixels(new Int32Rect(0, 0, width, height), pixels, stride, 0);
			return bmi;
		}
		//................................................................
		public static System.Windows.Media.Imaging.BitmapImage ToWpfImage(this System.Drawing.Image image) {
			if (image == null) return null;
			MemoryStream ms = new MemoryStream();  // no using here! BitmapImage will dispose the stream after loading
			image.Save(ms, System.Drawing.Imaging.ImageFormat.Bmp);

			System.Windows.Media.Imaging.BitmapImage bmi = new System.Windows.Media.Imaging.BitmapImage();
			bmi.BeginInit();
			bmi.CacheOption = System.Windows.Media.Imaging.BitmapCacheOption.OnLoad;
			bmi.StreamSource = ms;
			bmi.EndInit();
			return bmi;
		}
		//.................................................................
		public static System.Windows.Media.Imaging.BitmapImage Bitmap2Image(System.Drawing.Bitmap bmp, System.Drawing.Imaging.ImageFormat IF) {
			if (bmp == null) return null;
			System.Windows.Media.Imaging.BitmapImage bmi = new System.Windows.Media.Imaging.BitmapImage();
			using (MemoryStream memory = new MemoryStream()) {
				bmp.Save(memory, IF);
				memory.Position = 0;
				bmi.BeginInit();
				bmi.StreamSource = memory;
				bmi.CacheOption = System.Windows.Media.Imaging.BitmapCacheOption.OnLoad;
				bmi.EndInit();
			}
			return bmi;
		}
		//..................................................................
		public static System.Windows.Media.Imaging.BitmapImage GetBMI(string path) {
			if (!File.Exists(path)) return null;
			string se; System.Windows.Media.Imaging.BitmapImage bi = null;

			try {
				// Create source.
				bi = new System.Windows.Media.Imaging.BitmapImage();
				// BitmapImage.UriSource must be in a BeginInit/EndInit block.
				bi.BeginInit();
				bi.UriSource = new Uri(path, UriKind.RelativeOrAbsolute);
				bi.EndInit();
			} catch (Exception exp) {
				se = "Error while reading image file:\n" + path
				+ "\nEror code: " + exp.Message
				+ "\nStack: " + exp.StackTrace;
				MessageBox.Show(se);
				bi = null;
			}
			return bi;
		}
	}
}
