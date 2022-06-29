package wormTracker;

import java.awt.Frame;
import java.util.ArrayList;
import java.util.HashMap;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.io.SaveDialog;
import ij.measure.Calibration;
import ij.measure.Measurements;
import ij.measure.ResultsTable;
import ij.plugin.filter.BackgroundSubtracter;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.filter.PlugInFilter;
import ij.process.AutoThresholder;
import ij.process.Blitter;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import ij.text.TextWindow;

/**
	wrMTrck by Jesper Søndergaard Pedersen (JSP)
	
	Uses ImageJ's particle analyzer to track the movement of
	multiple objects through a stack. Uses the changes in aspect
	ratio of elipse fitting to quantify "thrashing" or
	swimming motion of C. elegans worms
	
	Based on MTrack2 by 2003/06/28 Nico Stuurman

	History

	Build 090408
		Implemented dos.newLine() instead of hard-coded "\n"
		Changed maxAreaChange to a percentage instead of pixel value
	Build 090707
		Added summarize feature
		Added Body-bend display on showPositions
		Added posibility to make an FPS override
		Now Works with scale-calibrated images
	Build 090804
		Fixed a bug in output of raw data.
	Build 091115
		Added x-y point smoothing function to decrease the "noise" movement of immobile objects
	Build 100201
		Corrected spelling mistake in input parameter "minTrackLength" (was "minTrackLenght") - this may affect scripts.
	Build 100411
		Corrected names in some fields on the dialog
	Build 100706
		Option to output plots of thrashing analysis for quality control purposes to rapidly asses if threshold value is good
	Build 101004
		Now works on both binary and thresholded movies.
		Posibility to enable on the fly deflickering and simple background subtraction based on RB15 on first frame.
		This means that wrMTrck can now track animals in raw movies and virtual stacks.
		Hitting the ESC key aborts analysis
	Build 110617
		Added an extra input field to allow adjustment of label font size
	Build 110622
		Added output of bend-time histogram (in frames) for each animal in text file output. Switch on/off histogram by adjusting bendDetect parameter
	    Build 111031
		Corrected undocumented raw data output with rawData =6
		With rawData <0 headers are suppressed.
		(Konstantine Palanski) Summary repors frames with maximum and minimum number of objects and the frame with maximum number of tracked objects
	Build 170303
		Added undocumented raw data output with rawData = 7
	Build 220107 (Kota)
			Mavenized
	TODO
	(Allow user to provide a thresholded image-stack instead of requiring binary-stack) - done
	Better support for multicore CPU to speed up analysis (Not sure how much can be gained)
	Make output of histograms in more user-friendly format
	Flexible user defineable point smoothing setting instead of hard coded 5-point smoothing
*/
public class wrMTrck_ implements PlugInFilter, Measurements {

	private boolean verbose = IJ.debugMode;
	ImagePlus imp;
	int nParticles;

	// KP
	int nMax = 0;
	int nMin = 99999;
	int nMaxFrm = 0;
	int nMinFrm = 0;
	// EOKP

	private TextWindow tw;
	private TextWindow tw2;
	private String directory, filename, rawFilename;

	private static String prevHdr;
	private String summaryHdr = "File\tnObj\tnObjFrm\tnFrames\tnTracks\ttotLength\tObjFrames\tObjSeconds\tavgSpeed\tavgArea\tavgPerim\tstdSpeed\tstdArea\tstdPerim"; // (KP)
	private static String prevhHdr;
	private String histogramHdr = "";
	private double pixelWidth = 1.0, pixelHeight = 1.0;
	private boolean suppressHeader;
	boolean done;

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		if (IJ.versionLessThan("1.17y"))
			return DONE;
		else
			return DOES_8G + NO_CHANGES;
	}

	public void run(ImageProcessor ip) {
		Parameters.initializeParametersGUI();

		/* Extract file name of movie file opened in ImageJ */
		if (Parameters.rawData < 0) {
			suppressHeader = true;
			Parameters.rawData = -Parameters.rawData;
		} else {
			suppressHeader = false;
		}
		;
		rawFilename = prepareRawFilename( imp );
		if (Parameters.bSaveResultsFile) {
			SaveDialog sd = new SaveDialog("Save Track Results", rawFilename, ".txt");
			directory = sd.getDirectory();
			filename = sd.getFileName();
		}
		Frame frame = WindowManager.getFrame("Summary"); //?
		if (done)  //?
			return;

		track(imp, directory, filename);
	}
	
	public String prepareRawFilename(ImagePlus imp) {
		String rawFilename = imp.getTitle();
		int dotPos = rawFilename.lastIndexOf('.');
		if (0 < dotPos && dotPos < rawFilename.length() - 2) {
			rawFilename = rawFilename.substring(0, dotPos);
		}
		return rawFilename;
	}

	public void track(ImagePlus imp, String directory, String filename) {
		int minSize = Parameters.minSize;
		int maxSize = Parameters.maxSize;
		float maxVelocity = Parameters.maxVelocity;

		int nFrames = imp.getStackSize();
		if (nFrames < 2) {
			IJ.showMessage("Tracker", "Stack required");
			return;
		}
		Calibration cal = imp.getCalibration();
		pixelWidth = cal.pixelWidth;
		pixelHeight = cal.pixelHeight;

		// Attempt to extract the frame rate for the movie if user did not provide a fps
		// value

		if (Parameters.fps <= 0) {
			Parameters.fps = cal.fps;
		}

		// Get information on the thresholding used

		ImageStack stack = imp.getStack();
		ImageProcessor ip = imp.getProcessor();
		double minThresh = ip.getMinThreshold();
		double maxThresh = ip.getMaxThreshold();

		if (verbose)
			IJ.log("min=" + minThresh + ",max=" + maxThresh);

		// Set options for particle analyser and measurements

		int options = (ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES
				// |ParticleAnalyzer.INCLUDE_HOLES
				| ParticleAnalyzer.SHOW_NONE); // set all PA options false
		int measurements = CENTROID + AREA + ELLIPSE + PERIMETER + CIRCULARITY;

		// Initialize results table
		ResultsTable rt = new ResultsTable();
		rt.reset();

		// Code required for deflickering and background subtraction on the fly
		double roiAvg[] = new double[nFrames + 1];
		BackgroundSubtracter bs = new BackgroundSubtracter();

		ip = stack.getProcessor(1);
		ip.resetRoi();
		if (verbose) {
			if (ip.isBinary()) {
				IJ.log("The image is binary");
			} else {
				IJ.log("The image is NOT binary");
			}
		}
		ImageProcessor bip = ip.duplicate();// .crop();
		if ((Parameters.backSub > 0) & !ip.isBinary()) { // subtract background if requested, but not if the image is
															// already binary.
			bs.rollingBallBackground(bip, 15, true, true, false, false, true);
		}

		
////////============ particle segmentation by thresholding 		
		IJ.showStatus("Finding all particles...");

		// record particle positions for each frame in an ArrayList	
		// create storage for particle positions
		// List[] theParticles = new ArrayList[nFrames];
		// supposed to be <framenumber starting from 1, an arraylist of particle
		// objects>
		HashMap<Integer, ArrayList<particle>> theParticles = new HashMap<Integer, ArrayList<particle>>();
		
		for (int iFrame = 1; iFrame <= nFrames; iFrame++) {			
			// Do frame normalization and background subtraction on the fly
			ip = stack.getProcessor(iFrame);
			// key steps for backgound subtraction
			if (Parameters.backSub > 0) {
				ImageStatistics is = ImageStatistics.getStatistics(ip, Measurements.MEAN, imp.getCalibration());
				roiAvg[iFrame] = is.mean;
				ip.resetRoi();
				ip.multiply(roiAvg[1] / roiAvg[iFrame]);
				ip.copyBits(bip, 0, 0, Blitter.DIFFERENCE);
				// bs.rollingBallBackground(ip, 15, false, true, false, false , true);
			}
			if (!(minThresh == -808080.0)) { // NO_THRESHOLD )) {
				ip.setThreshold(minThresh, maxThresh, 0);
				// IJ.log("preset threshold");
			} else if ((Parameters.backSub > 0) & !ip.isBinary()) {
				AutoThresholder at = new AutoThresholder();
				ip.setThreshold(at.getThreshold(Parameters.threshMode, ip.getHistogram()), 255, 0);

				minThresh = ip.getMinThreshold();
				maxThresh = ip.getMaxThreshold();
				if (verbose)
					IJ.log("AutoThresholding with " + Parameters.threshMode + ": min=" + minThresh + ",max="
							+ maxThresh);
			}

			// theParticles[iFrame - 1] = new ArrayList();

	////////============ particle localization 
			
			rt.reset();
			ParticleAnalyzer pa = new ParticleAnalyzer(options, measurements, rt, minSize, maxSize);
			if (!pa.analyze(imp, ip))
				return;
			float[] sxRes = rt.getColumn(ResultsTable.X_CENTROID);
			float[] syRes = rt.getColumn(ResultsTable.Y_CENTROID);
			float[] areaRes = rt.getColumn(ResultsTable.AREA);
			float[] angleRes = rt.getColumn(ResultsTable.ANGLE);
			float[] majorRes = rt.getColumn(ResultsTable.MAJOR);
			float[] minorRes = rt.getColumn(ResultsTable.MINOR);
			float[] perimRes = rt.getColumn(ResultsTable.PERIMETER);
			float[] circRes = rt.getColumn(ResultsTable.CIRCULARITY);

			ArrayList<particle> particleList = new ArrayList<particle>();
//			if (sxRes == null)
//				;
			// IJ.log("Frame"+iFrame+" contains no objects!");//return;// do not exit! there
			// may be objects on the next slice image!!!
			if (sxRes != null) {
				// IJ.log("Frame"+iFrame+" contains "+sxRes.length+" objects");
				for (int iPart = 0; iPart < sxRes.length; iPart++) {
					particle aParticle = new particle(sxRes[iPart], syRes[iPart], iFrame - 1, areaRes[iPart],
							angleRes[iPart], majorRes[iPart], minorRes[iPart], perimRes[iPart]);
					particleList.add(aParticle);
				}
				if (sxRes.length > nMax) {
					nMax = sxRes.length; // record the maximum number of particles found in one frame.
					nMaxFrm = iFrame; // record the frame at which the maximum number of particles is found. (KP)
				}
				if (sxRes.length < nMin) {
					nMin = sxRes.length; // record the minimum number of particles found in one frame. (KP)
					nMinFrm = iFrame; // record the frame at which the minimum number of particles is found. (KP)
				}
			}
			//frame number starts from 1
			theParticles.put(iFrame, particleList);
			
			IJ.showProgress((double) iFrame / nFrames);
			if (IJ.escapePressed()) {
				IJ.beep();
				done = true;
				return;
			}
			// IJ.log((int)iFrame+",");
		}

////////============ particle linking 
		
		// now assemble tracks out of the particle lists
		// Also record to which track a particle belongs in ArrayLists
		IJ.showStatus("Assembling tracks...");
		// IJ.log("Assembling tracks...");
		ArrayList<ArrayList<particle>> theTracks = particleLinking(nFrames, theParticles);

////////============ Analysis and Outputs
		
		// Initiate the writing of an output textfile if a filename is given.
		TrackAnalysis trackanal = new TrackAnalysis( directory, filename, imp, 
				theTracks, theParticles, suppressHeader);
		trackanal.setKPVlues( nMax, nMaxFrm, nMin, nMinFrm);
		trackanal.run();
		double[][] bendCounter = trackanal.getBendCounter();

////////============ Optional Outputs

		if (Parameters.bShowLabels) {
			IJ.showStatus("Generating movie with labels...");
			TrackNumberLabeling tnl = new TrackNumberLabeling(imp, theTracks, theParticles, bendCounter);
			ImagePlus nimp = tnl.run();
			if (nimp != null) {
			nimp.show();
			imp.show();
			nimp.updateAndDraw();
			} else {
				done = true;
				IJ.log("processing teminated during preparation of movie with labels");
			}
		}
		// 'map' of tracks
		if (Parameters.bShowPaths) {
			IJ.showStatus("Generating map of tracks...");
			TrackPlotter tpl = new TrackPlotter(imp, theTracks);
			tpl.run();
		}
	}

	public ArrayList<ArrayList<particle>> particleLinking(int nFrames, HashMap<Integer, ArrayList<particle>> theParticles) {
		ArrayList<ArrayList<particle>> theTracks = new ArrayList<ArrayList<particle>>();
		float maxVelocity = Parameters.maxVelocity;
		int trackCount = 0;
		for (int i = 0; i <= (nFrames - 1); i++) {
			IJ.showProgress((double) i / nFrames);
			if (IJ.escapePressed()) {
				IJ.beep();
				done = true;
				return null;
			}
			if (theParticles.size() > 0)
				for (particle aParticle : theParticles.get(i + 1)) {
					//particle aParticle = (particle) j.next();
					if (!aParticle.inTrack) {
						// This must be the beginning of a new track
						ArrayList<particle> aTrack = new ArrayList<particle>();
						trackCount++;
						aParticle.inTrack = true;
						aParticle.trackNr = trackCount;
						aTrack.add(aParticle);
						// search in next frames for more particles to be added to track
						boolean searchOn = true;
						particle oldParticle = new particle();
						particle tmpParticle = new particle();
						oldParticle.copy(aParticle);
						for (int iF = i + 1; iF <= (nFrames - 1); iF++) {
							boolean foundOne = false;
							particle newParticle = new particle();
							for (particle testParticle : theParticles.get(iF + 1)) {
								if (!searchOn)
									break;
								//particle testParticle = (particle) jF.next();
								float distance = testParticle.distance(oldParticle);
								float areaChange = 100 * Math.abs(1 - testParticle.area / oldParticle.area); // Calculate
																												// the
																												// area
																												// change
																												// in
																												// percent...

								// record a particle when it is within the search radius, and when it had not
								// yet been claimed by another track
								if ((distance / pixelWidth < maxVelocity) 
										&& (areaChange < Parameters.maxAreaChange)
										&& !testParticle.inTrack) {
									// if we had not found a particle before, it is easy
									if (!foundOne) {
										tmpParticle = testParticle;
										testParticle.inTrack = true;
										testParticle.trackNr = trackCount;
										testParticle.dist = distance;
										newParticle.copy(testParticle);
										foundOne = true;
									} else {
										// if we had one before, we'll take this one if it is closer. In any case, flag
										// these particles
										testParticle.flag = true;
										if (distance < newParticle.distance(oldParticle)) {
											testParticle.inTrack = true;
											testParticle.dist = distance;
											testParticle.trackNr = trackCount;
											newParticle.copy(testParticle);
											tmpParticle.inTrack = false;
											tmpParticle.trackNr = 0;
											tmpParticle.dist = 0;
											tmpParticle = testParticle;
										} else {
											newParticle.flag = true;
										}
									}
								} else if (distance / pixelWidth < maxVelocity) {
									// this particle is already in another track but could have been part of this
									// one
									// We have a number of choices here:
									// 1. Sort out to which track this particle really belongs (but how?)
									// 2. Stop this track
									// 3. Stop this track, and also delete the remainder of the other one
									// 4. Stop this track and flag this particle:
									testParticle.flag = true;
								}
							}
							if (foundOne)
								aTrack.add(newParticle);
							else
								searchOn = false;
							oldParticle.copy(newParticle);
						}
						theTracks.add(aTrack);
					}
				}
		}
		return theTracks;
	}
	
	// Utility functions

	int doOffset(int center, int maxSize, int displacement) {
		if ((center - displacement) < 2 * displacement) {
			return (center + 4 * displacement);
		} else {
			return (center - displacement);
		}
	}

	double s2d(String s) {
		Double d;
		try {
			d = new Double(s);
		} catch (NumberFormatException e) {
			d = null;
		}
		if (d != null)
			return (d.doubleValue());
		else
			return (0.0);
	}

	public static void main(String[] args) {
		String orgimagepath = "/Users/miura/Desktop/Celegans/sample01c-360-459crop.tif";
		ImagePlus imp = new ImagePlus(orgimagepath);
		wrMTrck_ wmt = new wrMTrck_();
		String  rawFilename = wmt.prepareRawFilename(imp);
		Parameters.minTrackLength = 10;
		Parameters.bShowPathLengths = true;
		Parameters.rawData = 2;
		Parameters.bSaveResultsFile = false;
		Parameters.bPlotBendTrack = true;
		Parameters.bShowPositions = true;
		wmt.suppressHeader = false;

		//SaveDialog sd = new SaveDialog("Save Track Results", rawFilename, ".txt");

		String directory = "/Users/miura/Desktop/testoutCUI/";
					// sd.getDirectory();
		String filename = rawFilename + ".txt"; //sd.getFileName();
		//IJ.log(directory);
		//IJ.log(filename);
		imp.getProcessor().setThreshold(0, 133, ImageProcessor.NO_LUT_UPDATE); //this should become adjustable. 
		wmt.track(imp, directory, filename);

	}
}
