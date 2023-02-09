//************************************************************************
//This package holds classes that can be used for Geo navigation
//This includes classes for sailboat racing
//Copyrgiht(C) Equin Technology, Inc. 2011
//************************************************************************

using System;
using System.Collections.Generic;
using System.Drawing.Printing;
using System.IO;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Security;
using System.Security.Permissions;
using System.Security.Principal;
using Microsoft.Win32;
using static System.Net.Mime.MediaTypeNames;

namespace SNL {
	//////////////////////////////////////////////////////////////////
	[Serializable]
	public class CSN_Asm {
		private readonly string sAC = "System.Reflection.AssemblyCopyrightAttribute";
		private readonly string sAV = "System.Reflection.AssemblyFileVersionAttribute";
		private readonly string sCA = "System.Reflection.AssemblyCompanyAttribute";
		private readonly string sPA = "System.Reflection.AssemblyProductAttribute";
		private readonly string sTT = "System.Reflection.AssemblyTitleAttribute";
		private readonly string sTM = "System.Reflection.AssemblyTrademarkAttribute";

		private string scopyright = null;
		private string sversion = null;
		private string sproduct = null;
		private string scompany = null;
		private string stitle = null;
		private string strademark = null;
		private bool b64 = false;
		private bool bOS64 = false;
		private readonly int tier = 0;
		private byte[] PKT = null;
		private IList<CustomAttributeData> CA = null;
		private string[] RNames;
		
		//============================================================
		public CSN_Asm() { Initialize(Assembly.GetExecutingAssembly()); }
		public CSN_Asm(Assembly Asmb) { Initialize(Asmb); }
		//==============================================================
		private void Initialize(Assembly Asmb) {
			PKT = Asmb.GetName().GetPublicKeyToken();
			CA = Asmb.GetCustomAttributesData();
			scopyright = finddatavalue(sAC);
			AssemblyName  aname = Asmb.GetName();
			Version ver = aname.Version;
			sversion =  ver.ToString();

			if (scopyright == null) scopyright = "Copyright (C) Steve Napear";
//			sversion = finddatavalue(sAV);
			scompany = finddatavalue(sCA);
			sproduct = finddatavalue(sPA);
			stitle = finddatavalue(sTT);
			strademark = finddatavalue(sTM);
			b64 = Environment.Is64BitProcess;
			bOS64 = Environment.Is64BitOperatingSystem;
			RNames = Asmb.GetManifestResourceNames();
		}
		//==============================================================
		private string finddatavalue(string sT) {        //find copyright notice from assembly
			int ix, i; string sCAD, sub;
			for (ix = 0; ix < CA.Count; ix++) {
				if (ix == 1) i = 0;
				CustomAttributeData CAD = CA[ix];
				sCAD = CAD.ToString();
				i = sCAD.IndexOf(sT);
				if (i < 0) continue;
				sub = sCAD.Substring(i + sT.Length);
				i = sub.IndexOf('\"');
				if (i < 0) continue;
				sub = sub.Substring(i + 1);
				i = sub.IndexOf('\"');
				if (i < 0) continue;
				return sub.Substring(0, i);
			}
			return null;
		}
		//==============================================================
		public string[] ResourceNames=>RNames;
		public string sCopyright { get { return scopyright; } }
		public byte[] GetPublicKey { get { return PKT; } }
		public string sVersion { get { return sversion; } }
		public string sTitle { get { return stitle; } }
		public string sTrademark { get { return strademark; } }
		public bool Is64 { get { return b64; } }
		public bool IsOS64 { get { return bOS64; } }
	}

	//////////////////////////////////////////////////////////////////
	[Serializable]
	public class CSN_Security {
		private readonly string sexedir;
		private readonly string suser;
		private readonly string sdomain;
		private readonly bool bsystem;
		private readonly string sauthtype;
		private readonly string sappname;
		private readonly bool badministrator;
		private readonly bool badmin;
		private readonly bool bsales;
		private readonly string st, suserroot = "HKEY_CURRENT_USER";

		//============================================================
		public bool isAdministrator { get { return badministrator; } }
		public bool isAdmin { get { return badmin; } }
		public bool isSales { get { return bsales; } }
		public bool isAuthenticated { get { return bsystem; } }
		public string sDomain { get { return sdomain; } }
		public string sUser { get { return suser; } }
		public string sExeDirectory { get { return sexedir; } }
		public string sDomainUser { get { return sdomain + '\\' + suser; } }
		public string sAppName { get { return sappname; } }
		public string skeyName { get { return suserroot + '\\' + sappname; } }
		public string sUserRoot { get { return suserroot; } }
		//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		public CSN_Security() {
			try {
				WindowsIdentity Wid = WindowsIdentity.GetCurrent();
				bsystem = Wid.IsAuthenticated;
				sauthtype = Wid.AuthenticationType;
				WindowsPrincipal pWP = new WindowsPrincipal(Wid);
				IIdentity id = pWP.Identity;
				sdomain = Environment.UserDomainName;
				badministrator = pWP.IsInRole(WindowsBuiltInRole.Administrator);
				badmin = pWP.IsInRole(sdomain + "\\Admin");
				bsales = pWP.IsInRole(sdomain + "\\Sales");
				st = pWP.Identity.Name;
				suser = Path.GetFileName(st);
				st = Environment.GetEnvironmentVariable("windir");
				string[] cla = Environment.GetCommandLineArgs();
				if (cla != null && cla.Length > 0) {
					st = Environment.GetCommandLineArgs()[0];
					sexedir = Path.GetDirectoryName(st);
					sappname = Path.GetFileNameWithoutExtension(st);
				} else {
					sexedir = "?";
					sappname = "?";
				}
			} catch {
				badministrator = false;
				bsales = false;
				badmin = false;
				suser = Environment.UserName;
				sexedir = "";
			} finally { }
		}
	}
	///////////////////////////////////////////////////////////////
	[Serializable]
	public class CSN_Printing {
		//======================================================================
		[DllImport("kernel32.dll")]
		static extern IntPtr GlobalLock(IntPtr hMem);

		[DllImport("kernel32.dll")]
		static extern bool GlobalUnlock(IntPtr hMem);

		[DllImport("kernel32.dll")]
		static extern bool GlobalFree(IntPtr hMem);

		[DllImport("winspool.Drv", EntryPoint = "DocumentPropertiesW",
					SetLastError = true, ExactSpelling = true,
					CallingConvention = CallingConvention.StdCall)]
		static extern int DocumentProperties(IntPtr hwnd, IntPtr hPrinter,
			[MarshalAs(UnmanagedType.LPWStr)] string pDeviceName,
			IntPtr pDevModeOutput, IntPtr pDevModeInput, int fMode);

		private const Int32 DM_OUT_BUFFER = 14;

		//========================================================================
		public static void OpenPrinterPropertiesDialog(PrinterSettings printerSettings, System.IntPtr pHandle) {
			if (!printerSettings.IsValid) return;
			IntPtr hDevMode = printerSettings.GetHdevmode();
			IntPtr pDevMode = GlobalLock(hDevMode);
			Int32 fMode = 0;
			int sizeNeeded = DocumentProperties(pHandle, IntPtr.Zero, printerSettings.PrinterName, pDevMode, pDevMode, fMode);
			IntPtr devModeData = Marshal.AllocHGlobal(sizeNeeded);

			fMode = DM_OUT_BUFFER;

			DocumentProperties(pHandle, IntPtr.Zero, printerSettings.PrinterName, devModeData, pDevMode, fMode);
			GlobalUnlock(hDevMode);
			printerSettings.SetHdevmode(devModeData);
			printerSettings.DefaultPageSettings.SetHdevmode(devModeData);
			GlobalFree(hDevMode);
			Marshal.FreeHGlobal(devModeData);
		}
	}
}
