package wormTracker;

import ij.gui.GenericDialog;
import ij.process.AutoThresholder;

public class Parameters {
	static	String 	buildNo = "Build 220107";
	static int 	fontSize = 16;
	static int	minSize = 100;
	static int	maxSize = 400;
	static float   	maxVelocity = 10;
	static float	maxAreaChange = 30;
	static int 	minTrackLength = 50;
	static float	minAngle = 2;
	static float	binSize = 0;
	static int	plotBendTrack=0;
	static boolean 	bSaveResultsFile = false;
	static boolean 	bShowLabels = true;
	static boolean 	bShowPositions = true;
	static boolean 	bShowPaths = false;
	static boolean 	bShowPathLengths = true;
	static boolean 	bShowSummary = true;
	static boolean 	bRoundCoord = false;
	static boolean  bSmoothing = true;
	static boolean  bPlotBendTrack = false;

	static int 	rawData = 0;
	static int	bendType = 2;
	static double	fps =0;
	static int	backSub = 0;
	static String	threshMode = "Otsu";
	static int 	maxColumns=75;

	//KP
	static 	int 	kpFrm = 0, kpMax = 0;
	static 	String row = new String();
	//EOKP
	public static void initializeParametersGUI() {
		GenericDialog gd = new GenericDialog("wrMTrck by Jesper S. Pedersen, "+ Parameters.buildNo);
		gd.addNumericField("minSize - Minimum Object Area (pixels^2): ", Parameters.minSize, 0);
		gd.addNumericField("maxSize - Maximum Object Area (pixels^2): ", Parameters.maxSize, 0);
		gd.addNumericField("maxVelocity - Maximum Velocity (pixels/frame):", Parameters.maxVelocity, 0);
		gd.addNumericField("maxAreaChange - Maximum area change (%):", Parameters.maxAreaChange, 0);
		gd.addNumericField("minTrackLength - Minimum track length (frames):", Parameters.minTrackLength, 0);
		gd.addNumericField("bendThreshold - Threshold for turn :", Parameters.minAngle,1);
		gd.addNumericField("binSize - Size of bin for speed histogram (pixels/frame) (0=disable):", Parameters.binSize,1);
		gd.addCheckbox("saveResultsFile - Save Results File:", Parameters.bSaveResultsFile);
		gd.addCheckbox("showPathLengths - Display Path Lengths:", Parameters.bShowPathLengths);
		gd.addCheckbox("showLabels - Show Labels:", Parameters.bShowLabels);
		gd.addCheckbox("showPositions - Show Positions:", Parameters.bShowPositions);
		gd.addCheckbox("showPaths - Show Paths:", Parameters.bShowPaths);
		gd.addCheckbox("showSummary - show a summary of tracking", Parameters.bShowSummary);
		gd.addCheckbox("roundCoord - round off coordinates", Parameters.bRoundCoord);
		gd.addCheckbox("smoothing - point smoothing", Parameters.bSmoothing);
		gd.addCheckbox("plotBendTrack - Quality control plots for thrashing analysis", Parameters.bPlotBendTrack);
		gd.addNumericField("rawData - (0=off,1=XYcord,2=Ellipse,3=AreaPerimDist,4=Ellipse+Circ,5=BendCalc):", Parameters.rawData, 0) ;
		gd.addNumericField("bendDetect - (0=Off,1=Angle,2=AspectRatio,3=AR+Histogram):", Parameters.bendType,0);
		gd.addNumericField("FPS - frames/s (0=try to load from file):", Parameters.fps,0 );
		gd.addNumericField("backSub - On-the-fly background subtraction (0=off,1=F1RB15):", Parameters.backSub,0);
		gd.addChoice("threshMode - Thresholding method (only if backSub>0)",AutoThresholder.getMethods(), "Otsu" );
		gd.addNumericField("fontSize - Size of labeling font:", Parameters.fontSize,0);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		Parameters.minSize = (int)gd.getNextNumber();
		Parameters.maxSize = (int)gd.getNextNumber();
		Parameters.maxVelocity = (float)gd.getNextNumber();
		Parameters.maxAreaChange = (float)gd.getNextNumber();
		Parameters.minTrackLength = (int)gd.getNextNumber();
		Parameters.minAngle = (float)gd.getNextNumber();
		Parameters.binSize = (float)gd.getNextNumber();
		Parameters.bSaveResultsFile = gd.getNextBoolean();
		Parameters.bShowPathLengths = gd.getNextBoolean();
		Parameters.bShowLabels = gd.getNextBoolean();
		Parameters.bShowPositions = gd.getNextBoolean();
		Parameters.bShowPaths = gd.getNextBoolean();
		Parameters.bShowSummary = gd.getNextBoolean();
		Parameters.bRoundCoord = gd.getNextBoolean();
		Parameters.bSmoothing= gd.getNextBoolean();
		Parameters.bPlotBendTrack = gd.getNextBoolean();
		Parameters.rawData = (int)gd.getNextNumber();
		Parameters.bendType = (int)gd.getNextNumber();
		Parameters.fps= (float)gd.getNextNumber();
		Parameters.backSub =(int)gd.getNextNumber();
		Parameters.threshMode = gd.getNextChoice();
		Parameters.fontSize = (int)gd.getNextNumber();
		if (Parameters.bShowPositions)
			Parameters.bShowLabels =true;
	}
}
