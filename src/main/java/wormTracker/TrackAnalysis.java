package wormTracker;

import java.awt.Color;
import java.awt.Frame;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.DoubleSummaryStatistics;
import java.util.HashMap;
import java.util.List;
import java.util.ListIterator;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.measure.ResultsTable;
import ij.process.ImageProcessor;
import ij.text.TextWindow;

// analyze tracks
// Discard tracks that are smaller than minTrackLength
// Calculate length of track
// Calculate bends
// Calculate how many animals are tracked on each frame

public class TrackAnalysis {

	ArrayList<ArrayList<particle>> theTracks;
	int trackCount;
	String directory;
	String filename;
	int nFrames;
	ImagePlus imp;
	HashMap<Integer, ArrayList<particle>> theParticles;
	boolean suppressHeader;
	double pixelWidth;

	// KP
	int nMax = 0;
	int nMin = 99999;
	int nMaxFrm = 0;
	int nMinFrm = 0;
	// EOKP

	String rawFilename;
	private static String prevHdr;
	private String summaryHdr = "File\tnObj\tnObjFrm\tnFrames\tnTracks\ttotLength\tObjFrames\tObjSeconds\tavgSpeed\tavgArea\tavgPerim\tstdSpeed\tstdArea\tstdPerim"; // (KP)
	private static String prevhHdr;
	private String histogramHdr = "";
	private TextWindow tw;
	private TextWindow tw2;

	//computed data
	double[][] angles;
	double[][] dAngles;
	double[][] sumDAngles;
	double[][] bendCounter;
	double[][] dAAreas;
	double[] times;
	double[][] bendTimes;
	
	//path length data, optional
	double[] lengths;
	double[] distances;
	double[] maxspeeds;
	double[] areas;
	double[] areaStdev;
	double[] perims;
	double[] perimsStdev;
	double[] bends;
	double[] avgX;
	double[] avgY;
	int[] frames;
	int[] firstFrames;
	
	int[][] bins; // used to store speed histograms
	int[][] bBins; // used to store bend-histograms
	int[] FrameTrackCount; // store how many objects are tracked in each frame
	

	public double[][] getBendCounter() {
		return bendCounter;
	}

	/**
	 * 
	 * @param directory : where the outputs are saved. 
	 * @param filename : file name of the original image.
	 * @param imp : ImagePlus object of the original image. 
	 * @param theTracks : a list of tracks (each track is an ArrayList of particles object)
	 * @param theParticles: HashMap with key = frame number (starting from 1) and value = ArrayList of detected particles in that frame.  
	 * @param suppressHeader : an output option, whether to ignore printing column headers. 
	 */
	public TrackAnalysis(String directory, String filename, ImagePlus imp, ArrayList<ArrayList<particle>> theTracks,
			HashMap<Integer, ArrayList<particle>> theParticles, boolean suppressHeader) {
		super();
		this.theTracks = theTracks;
		// this.trackCount = trackCount;
		this.directory = directory;
		this.filename = filename;
		// this.nFrames = nFrames;
		this.imp = imp;
		this.theParticles = theParticles;
		this.suppressHeader = suppressHeader;
		rawFilename = imp.getTitle();
		trackCount = theTracks.size();
		nFrames = imp.getStackSize();
		pixelWidth = imp.getCalibration().pixelWidth;
		
		angles = new double[trackCount + 1][nFrames];
		dAngles = new double[trackCount + 1][nFrames];
		sumDAngles = new double[trackCount + 1][nFrames];
		bendCounter = new double[trackCount + 1][nFrames];
		dAAreas = new double[trackCount + 1][nFrames];
		times = new double[nFrames];
		bendTimes = new double[trackCount + 1][nFrames];
		
		//path length data, optional
		lengths = new double[trackCount];
		distances = new double[trackCount];
		maxspeeds = new double[trackCount];
		areas = new double[trackCount];
		areaStdev = new double[trackCount];
		perims = new double[trackCount];
		perimsStdev = new double[trackCount];
		bends = new double[trackCount];
		avgX = new double[trackCount];
		avgY = new double[trackCount];
		frames = new int[trackCount];
		firstFrames = new int[trackCount];
		
	}

	public void setKPVlues(int nMax, int nMaxFrm, int nMin, int nMinFrm) {
		this.nMax = nMax;
		this.nMin = nMin;
		this.nMaxFrm = nMaxFrm;
		this.nMinFrm = nMinFrm;
	}

	/**
	 * initialize output textfile if it is not there, or exit if not writable. 
	 * @param outputfile
	 * @return
	 */
	public boolean checkOutputTextFile(File outputfile) {
		boolean writefile = false;
		if (!outputfile.canWrite()) {
			try {
				outputfile.createNewFile();
			} catch (IOException e) {
				IJ.showMessage("Error", "Could not create " + directory + filename);
				return false;
			}
		}
		if (outputfile.canWrite())
			writefile = true;
		else {
			IJ.showMessage("Error", "Could not write to " + directory + filename);
			return false;
		}
		return writefile;
	}
	/**
	 * the main routine
	 */
	public void run() {
		boolean writefile = false;
		// Initiate the writing of an output textfile if a filename is given.
		if (filename != null) {
			File outputfile = new File(directory, filename);
			writefile = checkOutputTextFile(outputfile);
		} 

		// Calculate length of paths, area of objects and number of body bends
		if (Parameters.bShowPathLengths) {
			if (Parameters.bSmoothing)
				trackSmoothener(theTracks, nFrames);
			
			int displayTrackNr = computeParameters();

			//Generate a stack of plots, one for thrashing analysis of each track showing 
			//the raw angle or aspect ratio as well as the kpSum of the dA/dt
			if (Parameters.bPlotBendTrack && Parameters.bendType > 0) {
				IJ.showStatus("Plotting bend calculation plots");
				ImagePlus psimp = plotBendCalculation( displayTrackNr );
				psimp.show();
				imp.show();
				psimp.updateAndDraw();
			}

			if (writefile)
				writeTrackDatatoFile( displayTrackNr );
			else
				showTrackDataInGUI(  displayTrackNr );
			
			// summarize the speed histogram for the movie
			if (Parameters.binSize > 0)
				showSpeedHistogramSummaryinTextWindow(displayTrackNr);

			//store the frame number with the maximum number of trackCounts. 
			compute_kpMax();

			// Generate summarized output
			// filename,N,nTracks,TotLength,totFrames,TotTime,AvgSpeed,AvgArea,AvgPerim,StdSpeed,StdArea,StdPerim
			if (Parameters.bShowSummary)
				showTrackSummaryinTextWindow( displayTrackNr );

		}

		//the following loop is only for updating trackCount, which now should become lower by:
		// eliminating short tracks. 
		trackCount = 1;
		for (ArrayList<particle> bTrack : theTracks)
			if (bTrack.size() >= Parameters.minTrackLength)
				trackCount++;

		if (!writefile && (Parameters.rawData > 0)) 
			showParticleDetailsinReulsts( trackCount );

		// and now when we write to file
		if (writefile && (Parameters.rawData > 0)) 
			writeParticleData( trackCount);

		//this.bendCounter = bendCounter;
	}

	/**
	 * generate smoothed coordinates to eliminate some digital noise - without
	 * coordinate smoothing lenght of track for immobile particles can actually be
	 * quite long for now only 5-point smoothing is possible, but the number of
	 * points should be user definable in later releases
	 */
	public void trackSmoothener(ArrayList<ArrayList<particle>> theTracks, int nFrames) {
		IJ.showStatus("Smoothing tracks...");
		double temp;
		double temp2;
		double[] rawX = new double[nFrames + 2];
		double[] rawY = new double[nFrames + 2];
		double[] smoothX = new double[nFrames + 2];
		double[] smoothY = new double[nFrames + 2];
//		for (ListIterator iT = theTracks.listIterator(); iT.hasNext();) {
//			List bTrack = (ArrayList) iT.next();
		for (ArrayList<particle> bTrack : theTracks) {
			int displayTrackNr = 0;
			if (bTrack.size() >= Parameters.minTrackLength) {
				displayTrackNr++;
				ListIterator jT = bTrack.listIterator();
				particle oldParticle = (particle) jT.next();
				particle oldParticle2 = oldParticle;
				rawX[1] = oldParticle.x;
				rawY[1] = oldParticle.y;
				smoothX[1] = oldParticle.x;
				smoothY[1] = oldParticle.y;
				oldParticle.sx = oldParticle.x;
				oldParticle.sy = oldParticle.y;
				int i = 1;
				int t = 0;
				for (; jT.hasNext();) {
					i++;
					particle newParticle = (particle) jT.next();
					rawX[i] = newParticle.x;
					rawY[i] = newParticle.y;
					newParticle.sx = newParticle.x;
					newParticle.sy = newParticle.y;
					smoothX[i] = newParticle.x;
					smoothY[i] = newParticle.y;
					if (i >= 3) {
						temp = 0;
						temp2 = 0;
						for (t = 0; t < 3; t++) {
							temp += rawX[t + i - 2];
							temp2 += rawY[t + i - 2];
						}
						smoothX[i - 1] = (float) temp / 3;
						smoothY[i - 1] = (float) temp2 / 3;
						oldParticle.sx = (float) temp / 3;
						oldParticle.sy = (float) temp2 / 3;

					}
					if (i >= 5) {
						temp = 0;
						temp2 = 0;
						for (t = 0; t < 5; t++) {
							temp += rawX[t + i - 4];
							temp2 += rawY[t + i - 4];
						}
						smoothX[i - 2] = (float) temp / 5;
						smoothY[i - 2] = (float) temp2 / 5;
						oldParticle2.sx = (float) temp / 5;
						oldParticle2.sy = (float) temp2 / 5;

					} else {
					}
					oldParticle2 = oldParticle;
					oldParticle = newParticle;
				}
			}

		}
	}
	
	public int computeParameters() {
		
		double distance, deltaAngle, oldAngle, sumDAngle, oldDAngle1, oldDAngle2, temp, temp2, avgArea, sumAreaSq,
		avgPerim, sumPerimSq, sumX, sumY, oldX, oldY;
		
		int binNum = 0;
		if (Parameters.binSize > 0)
			binNum = (int) (Parameters.maxVelocity / Parameters.binSize + 1);
		else
			binNum = 100;
		bins = new int[trackCount][(int) binNum]; // used to store speed histograms
		bBins = new int[trackCount][101]; // used to store bend-histograms
		FrameTrackCount = new int[nFrames]; // store how many objects are tracked in each frame

		int trackNr = 0;
		int displayTrackNr = 0;
		int bendCount = 0;
		int oldBendFrame = 0;
		// for (ListIterator iT = theTracks.listIterator(); iT.hasNext();) {
		for (List bTrack : theTracks) {
			trackNr++;
			// List bTrack = (ArrayList) iT.next();
			if (bTrack.size() >= Parameters.minTrackLength) {
				displayTrackNr++;
				maxspeeds[displayTrackNr - 1] = 0;
				ListIterator jT = bTrack.listIterator();
				particle oldParticle = (particle) jT.next();
				particle firstParticle = new particle();
				FrameTrackCount[firstParticle.z]++;
				firstParticle.copy(oldParticle);
				avgArea = oldParticle.area;
				avgPerim = oldParticle.perimeter;
				sumPerimSq = 0;
				sumAreaSq = 0;
				sumX = oldParticle.x;
				sumY = oldParticle.y;
				firstFrames[displayTrackNr - 1] = oldParticle.z;
				oldAngle = sumDAngle = oldDAngle1 = oldDAngle2 = 0;
				if (Parameters.bendType == 1)
					oldAngle = oldParticle.angle;
				if (Parameters.bendType > 1)
					oldAngle = oldParticle.ar;
				angles[displayTrackNr][0] = oldAngle;
				bendCount = 0;
				int i = 1;
				int t = 0;
				sumDAngle = oldDAngle1 = oldDAngle2 = 0;
				frames[displayTrackNr - 1] = bTrack.size();
				for (; jT.hasNext();) {
					i++;
					particle newParticle = (particle) jT.next();
					FrameTrackCount[newParticle.z]++; // add one object to FrameCount at the given frame
					if (Parameters.bSmoothing) {
						distance = newParticle.sdistance(oldParticle);
					} else
						distance = newParticle.distance(oldParticle);

					// Calculate average and standard deviation of object area
					temp = (avgArea + (newParticle.area - avgArea) / i);
					sumAreaSq += (newParticle.area - avgArea) * (newParticle.area - temp);
					avgArea = temp;

					// Calculate average and standard deviation of object perimeter
					temp = (avgPerim + (newParticle.perimeter - avgPerim) / i);
					sumPerimSq += (newParticle.perimeter - avgPerim) * (newParticle.perimeter - temp);
					avgPerim = temp;

					// Calculate the average position of the animals - e.g. useful if a single movie
					// cointains for seprate wells
					sumX += newParticle.x;
					sumY += newParticle.y;

					if (Parameters.bendType > 0) {

						// detect wobbling of particle or bending of objects

						if (Parameters.bendType == 1) { // use the angle of the ellipse for the fitting

							deltaAngle = newParticle.angle - oldAngle;
							oldAngle = newParticle.angle;

							if (deltaAngle > 100) {
								deltaAngle = deltaAngle - 180;
							}
							; // motion probably transverse 0-180 degree boundery
							if (deltaAngle < -100) {
								deltaAngle = deltaAngle + 180;
							}
							; // ....

						} else if (Parameters.bendType > 1) { // use the aspect ratio of the ellipse for counting
																// body-bends
							deltaAngle = newParticle.ar - oldAngle;
							oldAngle = newParticle.ar;
						} else
							deltaAngle = 0;
						angles[displayTrackNr][newParticle.z] = oldAngle;

						if (oldDAngle1 * deltaAngle <= 0) { // A change in sign indicates object started turning the
															// other direction direction
							if (sumDAngle > Parameters.minAngle) { // if area of last turn angle was over threshold
																	// degrees then count as a bend.
								binNum = newParticle.z - oldBendFrame;
								// IJ.log("binNum:"+binNum);
								if (binNum < 1)
									binNum = newParticle.z;
								if (binNum > 100)
									binNum = 100;
								bBins[displayTrackNr - 1][binNum]++;
								dAAreas[displayTrackNr][bendCount] = sumDAngle;
								bendTimes[displayTrackNr][bendCount] = newParticle.z;
								bendCount++;
								oldBendFrame = newParticle.z;
							}
							;
							sumDAngle = deltaAngle; // start summing new motion
						} else {
							sumDAngle += deltaAngle;
						}
						bendCounter[displayTrackNr][newParticle.z] = bendCount;
						sumDAngles[displayTrackNr][newParticle.z] = sumDAngle;

						oldDAngle1 = deltaAngle;

						dAngles[displayTrackNr][newParticle.z] = deltaAngle;
						times[newParticle.z] = newParticle.z;

					}
					// Generate a speed histogram
					if (!(newParticle.flag || oldParticle.flag)) { // only accecpt maximum speeds for non-flagged
																	// particles
						if (Parameters.binSize > 0) {
							int binNr = (int) (distance / pixelWidth / Parameters.binSize);
							bins[displayTrackNr - 1][binNr]++;
						}
						if (distance > maxspeeds[displayTrackNr - 1]) {
							maxspeeds[displayTrackNr - 1] = distance;
						}
						;
					}
					lengths[displayTrackNr - 1] += distance;
					oldParticle = newParticle;
				}
				distances[displayTrackNr - 1] = oldParticle.distance(firstParticle);// Math.sqrt(sqr(oldParticle.x-firstParticle.x)+sqr(oldParticle.y-firstParticle.y));
				areas[displayTrackNr - 1] = avgArea;
				areaStdev[displayTrackNr - 1] = Math.sqrt(sumAreaSq / (bTrack.size() - 1)); // Calculate Stdev of
																							// the areas
				perims[displayTrackNr - 1] = avgPerim;
				perimsStdev[displayTrackNr - 1] = Math.sqrt(sumPerimSq / (bTrack.size() - 1)); // Calculate the
																								// Stdev of the
																								// perimeters
				avgX[displayTrackNr - 1] = sumX / bTrack.size();
				avgY[displayTrackNr - 1] = sumY / bTrack.size();

				if (Parameters.bendType == 1)
					bends[displayTrackNr - 1] = (float) bendCount;
				else if (Parameters.bendType > 1)
					bends[displayTrackNr - 1] = (float) bendCount / 2;
			}
		}
		return displayTrackNr;
	}
	/**
	 * @TODO replace the Plot constructor https://javadoc.io/doc/net.imagej/ij/latest/ij/ij/gui/Plot.html
	 * @param displayTrackNr
	 * @return
	 */
	public ImagePlus plotBendCalculation(int displayTrackNr) {
		// plot the first bend calculation for the first track
		Parameters.plotBendTrack = 1;
		String YaxisLabel = "Aspect ratio";
		if (Parameters.bendType == 1)
			YaxisLabel = "Degrees";
		Plot plot = new Plot("Bend detection plot for track " + (int) Parameters.plotBendTrack, "Frame#",
				YaxisLabel, times, sumDAngles[Parameters.plotBendTrack]); // angles
		plot.setSize(nFrames + 150, 300);
		if (Parameters.bendType > 1)
			plot.setLimits(firstFrames[Parameters.plotBendTrack - 1],
					firstFrames[Parameters.plotBendTrack - 1] + frames[Parameters.plotBendTrack - 1], 0,
					Math.round(max(angles[Parameters.plotBendTrack]) + 2));
		if (Parameters.bendType == 1)
			plot.setLimits(firstFrames[Parameters.plotBendTrack - 1],
					firstFrames[Parameters.plotBendTrack - 1] + frames[Parameters.plotBendTrack - 1], -200,
					200);
		plot.setColor(Color.red);
		plot.addPoints(bendTimes[Parameters.plotBendTrack], dAAreas[Parameters.plotBendTrack], Plot.X);
		plot.setColor(Color.green);
		plot.addPoints(times, angles[Parameters.plotBendTrack], PlotWindow.LINE); // dAngles sumDAngles
		float x1[] = { 0, nFrames };
		float y1[] = { Parameters.minAngle, Parameters.minAngle };
		plot.setColor(Color.blue);
		plot.addPoints(x1, y1, PlotWindow.LINE);
		plot.setColor(Color.black);
		plot.draw();
		ImagePlus pimp = plot.getImagePlus();
		ImageStack plotstack = pimp.createEmptyStack();
		ImageProcessor nip = plot.getProcessor();
		plotstack.addSlice("Bend detection plot for track " + (int) Parameters.plotBendTrack, nip.crop());
		ImageProcessor pip = plotstack.getProcessor(1);

		// Repeat for the rest of the plots
		for (int i = Parameters.plotBendTrack; i < displayTrackNr; i++) {
			Parameters.plotBendTrack++;
			IJ.showProgress((double) i / displayTrackNr);
			if (IJ.escapePressed()) {
				IJ.beep();
				// done = true;
				return null;
			}
			Plot plot2 = new Plot("Bend detection plot for track " + (int) Parameters.plotBendTrack, "Frame#",
					YaxisLabel, times, sumDAngles[Parameters.plotBendTrack]); // angles
			plot2.setSize(nFrames + 150, 300);
			if (Parameters.bendType > 1)
				plot2.setLimits(firstFrames[Parameters.plotBendTrack - 1],
						firstFrames[Parameters.plotBendTrack - 1] + frames[Parameters.plotBendTrack - 1], 0,
						Math.round(max(angles[Parameters.plotBendTrack]) + 2));
			if (Parameters.bendType == 1)
				plot2.setLimits(firstFrames[Parameters.plotBendTrack - 1],
						firstFrames[Parameters.plotBendTrack - 1] + frames[Parameters.plotBendTrack - 1], -200,
						200);
			plot2.setColor(Color.red);
			plot2.addPoints(bendTimes[Parameters.plotBendTrack], dAAreas[Parameters.plotBendTrack], Plot.X);
			plot2.setColor(Color.green);
			plot2.addPoints(times, angles[Parameters.plotBendTrack], PlotWindow.LINE); // dAngles sumDAngles
			plot2.setColor(Color.blue);
			plot2.addPoints(x1, y1, PlotWindow.LINE);
			plot2.setColor(Color.black);
			plot2.draw();
			nip = plot2.getProcessor();
			plotstack.addSlice("Bend detection plot for track " + (int) Parameters.plotBendTrack, nip.crop());
		}

		ImagePlus psimp = new ImagePlus(imp.getTitle() + " plots", plotstack);		
		return psimp;
	}
	
	public String generateSummarizedOutputString(int displayTrackNr) {
		float sumLengths = 0;
		float sumFrames = 0;
		//float sumTimes = 0;
		double avgArea = 0;
		double sumAreaSq = 0;
		double avgPerim = 0;
		double sumPerimSq = 0;
		double avgSpeed = 0;
		double sumSpeedSq = 0;
		double speed = 0;
		double sumBends = 0;
		double avgBBPS = 0;
		double sumBBPSSq = 0;
		double BBPS = 0;
		double temp = 0;

		for (int i = 0; i < displayTrackNr; i++) {
			sumLengths += lengths[i];
			sumFrames += frames[i];
			sumBends += bends[i];

			// Calculate average and standard deviation of object area
			temp = (avgArea + (areas[i] - avgArea) / (i + 1));
			sumAreaSq += (areas[i] - avgArea) * (areas[i] - temp);
			avgArea = temp;

			// Calculate average and standard deviation of object perimeter
			temp = (avgPerim + (perims[i] - avgPerim) / (i + 1));
			sumPerimSq += (perims[i] - avgPerim) * (perims[i] - temp);
			avgPerim = temp;

			// Calculate average and standard deviation of object Speed
			speed = lengths[i] / (frames[i] / Parameters.fps);
			temp = (avgSpeed + (speed - avgSpeed) / (i + 1));
			sumSpeedSq += (speed - avgSpeed) * (speed - temp);
			avgSpeed = temp;

			// Calculate average and standard deviation of BodyBends per seconds
			BBPS = bends[i] / (frames[i] / Parameters.fps);
			temp = (avgBBPS + (BBPS - avgBBPS) / (i + 1));
			sumBBPSSq += (BBPS - avgBBPS) * (BBPS - temp);
			avgBBPS = temp;

		}
		String aLine = null;
		summaryHdr = "File\tnObjMax\tnObjMaxFrm\tnObjMin\tnObjMinFrm\tkpMax\tkpMaxFrm\tnFrames\tnTracks\ttotLength\tObjFrames\tObjSeconds\tavgSpeed\tavgArea\tavgPerim\tstdSpeed\tstdArea\tstdPerim";
		if (Parameters.bendType > 0)
			summaryHdr += "\tBends\tavgBBPS\tstdBBPS";

		aLine = rawFilename + "\t" + (int) nMax // number of objects (N-value)
				+ "\t" + (int) nMaxFrm // frame at nMax (KP)
				+ "\t" + (int) nMin // nMin (KP)
				+ "\t" + (int) nMinFrm // frame at nMin (KP)
				+ "\t" + (int) Parameters.kpMax // kpMax objects (KP)
				+ "\t" + (int) Parameters.kpFrm // frame at kpMax objects (KP)
				+ "\t" + (int) nFrames // number of frames
				+ "\t" + (int) displayTrackNr // number of tracks
				+ "\t" + (float) sumLengths // total distance covered by all objects
				+ "\t" + (float) sumFrames // total number of object*frames
				+ "\t" + (float) sumFrames / Parameters.fps // total obj*seconds
				+ "\t" + (float) avgSpeed // average speed (pixels/seconds)
				+ "\t" + (float) avgArea // average worm area
				+ "\t" + (float) avgPerim // average worm perimeter
				+ "\t" + (float) Math.sqrt(sumSpeedSq / (displayTrackNr - 1))// average speed (pixels/seconds)
				+ "\t" + (float) Math.sqrt(sumAreaSq / (displayTrackNr - 1))// average worm area
				+ "\t" + (float) Math.sqrt(sumPerimSq / (displayTrackNr - 1)) // average worm perimeter
		;
		if (Parameters.bendType > 0) {
			aLine += "\t" + (float) sumBends // total number of body bends in the movie
					+ "\t" + (float) avgBBPS // average BBPS per track
					+ "\t" + (float) Math.sqrt(sumBBPSSq / (displayTrackNr - 1)) // standard deviation of BBPS
																					// per track
			;

		}	
		return aLine;
	}
	
	public void writeTrackDatatoFile(int displayTrackNr) {
		try {
			File outputfile = new File(directory, filename);
			BufferedWriter dos = new BufferedWriter(new FileWriter(outputfile)); // append
			if (Parameters.bendType > 0) {
				dos.write(
						"Track \tLength\tDistance\t#Frames\t1stFrame\ttime(s)\tMaxSpeed\tAvgArea\tStdArea\tAvgPerim\tStdPerim\tAvgSpeed\tBLPS\tavgX\tavgX\tBends\tBBPS");
			} else {
				dos.write(
						"Track \tLength\tDistance\t#Frames\t1stFrame\ttime(s)\tMaxSpeed\tAvgArea\tStdArea\tAvgPerim\tStdPerim\tAvgSpeed\tBLPS\tavgX\tavgX");
			}

			dos.newLine();
			for (int i = 0; i < displayTrackNr; i++) {
				String str = "" + (i + 1) + "\t" + (float) lengths[i] + "\t" + (float) distances[i] + "\t"
						+ (int) frames[i] + "\t" + (int) firstFrames[i] + "\t"
						+ (float) (frames[i] / Parameters.fps) + "\t" + (float) Parameters.fps * maxspeeds[i]
						+ "\t" + (float) areas[i] + "\t" + (float) areaStdev[i] + "\t" + (float) perims[i]
						+ "\t" + (float) perimsStdev[i] + "\t"
						+ (float) (Parameters.fps * lengths[i] / frames[i]) + "\t"
						+ (float) (Parameters.fps * 2 * lengths[i] / perims[i] / frames[i]) + "\t"
						+ (float) avgX[i] + "\t" + (float) avgY[i];

				if (Parameters.bendType > 0) {
					str += "\t" + (float) bends[i] + "\t" + (float) Parameters.fps * bends[i] / frames[i];
				}

				dos.write(str);
				dos.newLine();
			}
			if (Parameters.binSize > 0) {
				dos.newLine();
				dos.write("Bin\t");
				for (int t = 0; t < displayTrackNr; t++) {
					dos.write("T#" + (int) (t + 1) + "\t");
				}
				;
				dos.newLine();

				for (int i = 0; i < (int) (Parameters.maxVelocity / Parameters.binSize); i++) {
					String str = "" + (float) i * Parameters.binSize + "\t";
					for (int t = 0; t < displayTrackNr; t++) {
						str = str + bins[t][i] + "\t";
					}
					;
					dos.write(str);
					dos.newLine();
				}
			}
			if (Parameters.bendType > 2) {
				dos.newLine();
				String str = "Bin";
				for (int t = 0; t < displayTrackNr; t++) {
					str = str + ("\tTrack" + (int) (t + 1));
				}
				;
				dos.write(str);
				dos.newLine();
				for (int i = 0; i < 101; i++) {
					int mySum = 0;
					str = "" + (float) i;
					for (int t = 0; t < displayTrackNr; t++) {
						str += "\t" + bBins[t][i];
						// mySum+=bBins[t][i];
					}
					;
					dos.write(str);
					dos.newLine();
				}

			}
			dos.close();
		} catch (IOException e) {
			if (filename != null)
				IJ.error("An error occurred writing the file. \n \n " + e);
		}	
	}
	
	public void showTrackDataInGUI( int displayTrackNr) {
		ResultsTable mrt = ResultsTable.getResultsTable();
		mrt.reset();

		for (int i = 0; i < displayTrackNr; i++) {
			mrt.incrementCounter();
			mrt.addValue("Track", (i + 1));
			mrt.addValue("Length", (float) lengths[i]);
			mrt.addValue("Distance", (float) distances[i]);
			mrt.addValue("#Frames", (float) frames[i]);
			mrt.addValue("1stFrame", (float) firstFrames[i]);
			mrt.addValue("Time(s)", (float) frames[i] / Parameters.fps);
			mrt.addValue("MaxSpeed", (float) Parameters.fps * maxspeeds[i]);
			mrt.addValue("Area", (float) areas[i]);
			mrt.addValue("sdArea", (float) areaStdev[i]);
			mrt.addValue("Perim", (float) perims[i]);
			mrt.addValue("sdPerim", (float) perimsStdev[i]);
			mrt.addValue("avgSpeed", (float) (Parameters.fps * lengths[i] / frames[i]));
			mrt.addValue("BLPS", (float) (Parameters.fps * 2 * lengths[i] / perims[i] / frames[i]));
			mrt.addValue("avgX", (float) avgX[i]);
			mrt.addValue("avgY", (float) avgY[i]);
			if (Parameters.bendType > 0) {
				mrt.addValue("Bends", (float) bends[i]);
				mrt.addValue("BBPS", (float) Parameters.fps * bends[i] / frames[i]);
			}
		}
		mrt.show("Results");

		if (Parameters.binSize > 0) {

			// Write histogram for individual track in the movie
			IJ.write("");
			String str = "Bin\t";
			for (int t = 0; t < displayTrackNr; t++) {
				str = str + ("Track" + (int) (t + 1) + "\t");
			}
			;
			str = str + "\n";
			IJ.write(str);

			for (int i = 0; i < (int) (Parameters.maxVelocity / Parameters.binSize); i++) {
				str = "" + (float) i * Parameters.binSize + "\t";
				for (int t = 0; t < displayTrackNr; t++) {
					str = str + bins[t][i] + "\t";
				}
				;
				IJ.write(str + "\n");

			}
		}
		// thrashing histogram
		if (Parameters.bendType > 2) {
			IJ.write("");
			String str = "#Frames";
			for (int t = 0; t < displayTrackNr; t++) {
				str = str + ("\tTrack" + (int) (t + 1));
			}
			;
			IJ.write(str + "\n");
			for (int i = 0; i < 101; i++) {
				int mySum = 0;
				str = "" + (float) i;
				for (int t = 0; t < displayTrackNr; t++) {
					str += "\t" + bBins[t][i];
					// mySum+=bBins[t][i];
				}
				;
				IJ.write(str /* +(int)mySum */ + "\n");
			}

		}		
	}
	
	public void showSpeedHistogramSummaryinTextWindow(int displayTrackNr) {
		histogramHdr = "Bin\t";
		String str2 = rawFilename + "\t";
		for (int i = 0; i < (int) (Parameters.maxVelocity / Parameters.binSize); i++) {
			int mySum = 0;
			histogramHdr += (float) Parameters.fps * pixelWidth * i * Parameters.binSize + "\t"; // pixels/frame
																									// *
																									// ｵm/pixels
																									// *frame/s
																									// pixelWidth*
			for (int t = 0; t < displayTrackNr; t++) {
				mySum += bins[t][i];
			}
			; // str = str +bins[t][i]+"\t";};
			str2 += (float) mySum + "\t";
		}

		Frame frame = WindowManager.getFrame("SpeedHistogram");
		if (frame != null && (frame instanceof TextWindow) && histogramHdr.equals(prevhHdr))
			tw2 = (TextWindow) frame;

		if (tw2 == null) {
			String title = "SpeedHistogram";
			tw2 = new TextWindow(title, histogramHdr, str2, 450, 300);
			prevhHdr = histogramHdr;
		} else
			tw2.append(str2);		
	}

	public void compute_kpMax() {
		/*
		 * //KP MAX CONCURRENT VALUES // Matrix that is nFrames for columns, and
		 * displayTrackNr for Rows int[][] kpMatrix = new int[nFrames][displayTrackNr];
		 * 
		 * // Calculates the the lastFrame by adding each frame value to each firstFrame
		 * value int[] lastFrames = new int[displayTrackNr]; for (int i=0; i <
		 * displayTrackNr; i++) { lastFrames[i] = firstFrames[i] + frames [i]; }
		 * //Populates the kpMatrix with 1s from firstFrames[n] to lastFrames[n] one row
		 * at a time for(int n = 0; n < displayTrackNr; n++ ) { for(int i =
		 * firstFrames[n]; i < lastFrames[n]; i++) { kpMatrix[i][n] = 1; } }
		 */
		// Searches through kpMatrix one column at a time to find the maximum value,
		// stores it in kpMax and kpFrm
		int kpSum;
		Parameters.kpMax = 0;
		for (int i = 1; i < FrameTrackCount.length; i++) {
			kpSum = 0;
//		for(int n = 0; n < displayTrackNr; n++) {
//			kpSum += kpMatrix[i][n];
			kpSum = FrameTrackCount[i];
//			if (n == displayTrackNr-1) {
			if (kpSum > Parameters.kpMax) {
				Parameters.kpMax = kpSum;
				Parameters.kpFrm = i;
//			}
//			}
			}
			// IJ.log("frame,"+(int) i+ "," + (int)kpSum);
		}

		// For debugging. Outputs raw 1 and 0s to console
		// for(int i = 0; i < displayTrackNr; i++) {
		// row = "";
		// for(int n = 0; n < nFrames; n++) {
		// row = row + kpMatrix[n][i];
		// if (n == nFrames-1) {
		// System.out.println(row);
		// }
		// }
		// }
		// EOKP		
	}
	
	public void showTrackSummaryinTextWindow( int displayTrackNr ) {
		String aLine = generateSummarizedOutputString( displayTrackNr);
		// IJ.log(aLine+"\n");
		Frame frame = WindowManager.getFrame("Summary");
		if (frame != null && (frame instanceof TextWindow) && summaryHdr.equals(prevHdr))
			tw = (TextWindow) frame;

		if (tw == null) {
			String title = "Summary";
			tw = new TextWindow(title, summaryHdr, aLine, 450, 300);
			prevHdr = summaryHdr;
		} else
			tw.append(aLine);		
	}
	//@TODO
	public String generateParticleTableHeadings() {
		// Output raw values based on user selection.
		// Create the column headings based on the number of tracks
		// with length greater than minTrackLength
		// since the number of tracks can be larger than can be accomodated by Excell,
		// we deliver the tracks in chunks of maxColumns
		// As a side-effect, this makes the code quite complicated

		// display the table with particle positions
		// first when we only write to the screen

		String strHeadings = "Frame";
		// @ToDo kota: this trackCount is better be a different variable as it counts tracks longer
		// than minTracklength
		// in addition, the output file can be a CSV file - and do not separate data for each track. 
		// rather add a track ID number as another column. long and narrow table is preferred rather than short and wide for data wrangling. 
		int trackCount = 1;
		//for (ListIterator iT = theTracks.listIterator(); iT.hasNext();) {
		for (List bTrack : theTracks) {
			//List bTrack = (ArrayList) iT.next();
			if (bTrack.size() >= Parameters.minTrackLength) {
				if (trackCount <= Parameters.maxColumns) {
					if (Parameters.rawData == 1)
						strHeadings += "\tX" + trackCount + "\tY" + trackCount + "\tFlag" + trackCount;
					if (Parameters.rawData == 2)
						strHeadings += "\tMajor" + trackCount + "\tMinor" + trackCount + "\tAngle" + trackCount;
					if (Parameters.rawData == 3)
						strHeadings += "\tArea" + trackCount + "\tPerimeter" + trackCount + "\tDistance" + trackCount;
					if (Parameters.rawData == 4)
						strHeadings += "\tMajor" + trackCount + "\tMinor" + trackCount + "\tCircularity" + trackCount;
					if (Parameters.rawData == 5)
						strHeadings += "\tShape" + trackCount + "\tRateOfChange" + trackCount + "\tSumOfChange"
								+ trackCount;
					if (Parameters.rawData == 6)
						strHeadings += "\tX" + trackCount + "\tY" + trackCount + "\tArea" + trackCount + "\tPerim"
								+ trackCount + "\tAngle" + trackCount + "\tAR" + trackCount + "\tFlag" + trackCount;
					if (Parameters.rawData == 7)
						strHeadings += "\tX" + trackCount + "\tY" + trackCount + "\tMajor" + trackCount + "\tMinor"
								+ trackCount;
				}
				trackCount++;
			}
		}
		return strHeadings;
	}
	
	public void writeParticleData( int trackCount ) {
		String flags;
		String strHeadings = generateParticleTableHeadings();
		//int repeat = (int) (Math.ceil(trackCount / Parameters.maxColumns));
		int repeat=(int) ( (trackCount/Parameters.maxColumns) );
		float reTest = (float) trackCount/ (float) Parameters.maxColumns;
		if (reTest > repeat)
			repeat++;
		try {

			File outputfile = new File(directory, filename.substring(0, filename.length() - 4) + "_raw.txt");

//		File outputfile=new File (directory,filename);

			BufferedWriter dos = new BufferedWriter(new FileWriter(outputfile, true));
			if (!suppressHeader) {
				dos.write(strHeadings);
				dos.newLine();
			}
			for (int j = 1; j <= repeat; j++) {
				int to = j * Parameters.maxColumns;
				if (to > trackCount - 1)
					to = trackCount - 1;
				String stLine = "Tracks " + ((j - 1) * Parameters.maxColumns + 1) + " to " + to;
				if (!suppressHeader) {
					dos.write(stLine);
					dos.newLine();
				}
				for (int i = 0; i <= (nFrames - 1); i++) {
					String strLine = "" + (i + 1);
					int trackNr = 0;
					int listTrackNr = 0;
					for (ListIterator iT = theTracks.listIterator(); iT.hasNext();) {
						trackNr++;
						List bTrack = (ArrayList) iT.next();
						boolean particleFound = false;
						if (bTrack.size() >= Parameters.minTrackLength) {
							listTrackNr++;
							if ((listTrackNr > ((j - 1) * Parameters.maxColumns))
									&& (listTrackNr <= (j * Parameters.maxColumns))) {
								for (ListIterator k = theParticles.get(i + 1).listIterator(); k.hasNext()
										&& !particleFound;) {
									particle aParticle = (particle) k.next();
									if (aParticle.trackNr == trackNr) {
										particleFound = true;
										if (aParticle.flag)
											flags = "*";
										else
											flags = " ";
										if (Parameters.rawData == 3) {
											strLine += "\t" + aParticle.area + "\t" + aParticle.perimeter + "\t"
													+ aParticle.dist;
										} else if (Parameters.rawData == 5) {
											strLine += "\t" + angles[listTrackNr][i] + "\t"
													+ dAngles[listTrackNr][i] + "\t" + sumDAngles[listTrackNr][i];
										} else if (Parameters.rawData == 6) {
											strLine += "\t" + aParticle.x + "\t" + aParticle.y + "\t"
													+ aParticle.area + "\t" + aParticle.perimeter + "\t"
													+ aParticle.angle + "\t" + aParticle.ar + "\t" + flags;
										} else if (Parameters.rawData == 7) {
											strLine += "\t" + aParticle.x + "\t" + aParticle.y + "\t"
													+ aParticle.major + "\t" + aParticle.minor;
										} else if (Parameters.rawData == 4) {
											strLine += "\t" + aParticle.major + "\t" + aParticle.minor + "\t"
													+ aParticle.circularity;
										} else if (Parameters.rawData == 2) {
											strLine += "\t" + aParticle.major + "\t" + aParticle.minor + "\t"
													+ aParticle.angle;
										} else {
											strLine += "\t" + aParticle.x + "\t" + aParticle.y + "\t" + flags;

										}
									}
								}
								if (!particleFound)
									strLine += "\t\t\t";
							}
						}
					}
					dos.write(strLine);
					dos.newLine();
				}
			}

			dos.newLine();
			dos.close();
		} catch (IOException e) {
			if (filename != null)
				IJ.error("An error occurred writing the file. \n \n " + e);
		}
	}
	
	public void showParticleDetailsinReulsts(int trackCount) {
		float flag;
//		int repeat = (int) (Math.ceil(trackCount / Parameters.maxColumns));
		int repeat=(int) ( (trackCount/Parameters.maxColumns) );
		float reTest = (float) trackCount/ (float) Parameters.maxColumns;
		if (reTest > repeat)
			repeat++;
		ResultsTable rrt = new ResultsTable();
		for (int j = 1; j <= repeat; j++) {
			int to = j * Parameters.maxColumns;
			if (to > trackCount - 1)
				to = trackCount - 1;
			rrt.reset();
			String stLine = "Raw Tracks " + ((j - 1) * Parameters.maxColumns + 1) + " to " + to;
			// IJ.write(stLine);
			for (int i = 0; i <= (nFrames - 1); i++) {
				// String strLine = "" + (i+1);
				int trackNr = 0;
				int listTrackNr = 0;
				rrt.incrementCounter();
				for (ListIterator iT = theTracks.listIterator(); iT.hasNext();) {
					trackNr++;
					List bTrack = (ArrayList) iT.next();
					boolean particleFound = false;
					if (bTrack.size() >= Parameters.minTrackLength) {
						listTrackNr++;
						if ((listTrackNr > ((j - 1) * Parameters.maxColumns))
								&& (listTrackNr <= (j * Parameters.maxColumns))) {
							if (theParticles.size() > 0)
								for (ListIterator k = theParticles.get(i + 1).listIterator(); k.hasNext()
										&& !particleFound;) {
									particle aParticle = (particle) k.next();
									if (aParticle.trackNr == trackNr) {
										particleFound = true;
										if (aParticle.flag)
											flag = 1;
										else
											flag = 0;
										if (Parameters.rawData == 3) {
											rrt.addValue("Area" + listTrackNr, (float) aParticle.area);
											rrt.addValue("Perimeter" + listTrackNr, (float) aParticle.perimeter);
											rrt.addValue("Distance" + listTrackNr, (float) aParticle.dist);
										} else if (Parameters.rawData == 4) {
											rrt.addValue("Major" + listTrackNr, (float) aParticle.major);
											rrt.addValue("Minor" + listTrackNr, (float) aParticle.minor);
											rrt.addValue("Circularity" + listTrackNr,
													(float) aParticle.circularity);
										} else if (Parameters.rawData == 5) {
											rrt.addValue("Shape" + listTrackNr, (float) angles[listTrackNr][i]);
											rrt.addValue("RateOfChange" + listTrackNr,
													(float) dAngles[listTrackNr][i]);
											rrt.addValue("SumOfChange" + listTrackNr,
													(float) sumDAngles[listTrackNr][i]);
										} else if (Parameters.rawData == 6) {
											rrt.addValue("X" + listTrackNr, (float) aParticle.x);
											rrt.addValue("Y" + listTrackNr, (float) aParticle.y);
											rrt.addValue("Area" + listTrackNr, (float) aParticle.area);
											rrt.addValue("Perim" + listTrackNr, (float) aParticle.perimeter);
											rrt.addValue("Angle" + listTrackNr, (float) aParticle.angle);
											rrt.addValue("AR" + listTrackNr, (float) aParticle.ar);
											rrt.addValue("Flag" + listTrackNr, flag);// aParticle.flag);
										} else if (Parameters.rawData == 7) {
											rrt.addValue("X" + listTrackNr, (float) aParticle.x);
											rrt.addValue("Y" + listTrackNr, (float) aParticle.y);
											rrt.addValue("Major" + listTrackNr, (float) aParticle.major);
											rrt.addValue("Minor" + listTrackNr, (float) aParticle.minor);
										} else if (Parameters.rawData == 2) {
											rrt.addValue("Major" + listTrackNr, (float) aParticle.major);
											rrt.addValue("Minor" + listTrackNr, (float) aParticle.minor);
											rrt.addValue("Angle" + listTrackNr, (float) aParticle.angle);
										} else {
											rrt.addValue("X" + listTrackNr, (float) aParticle.x);
											rrt.addValue("Y" + listTrackNr, (float) aParticle.y);
											rrt.addValue("Flag" + listTrackNr, (int) flag);// aParticle.flag);
										}
									}
								}
						}
					}
				}
			}
			rrt.show(stLine);
		}		
	}

	// ===================================================== max
	public static double max(double[] t) {
		DoubleSummaryStatistics stat = Arrays.stream(t).summaryStatistics();
		double maximum = stat.getMax();
		return maximum;
	}// end method max

}
