package wormTracker;

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
}
