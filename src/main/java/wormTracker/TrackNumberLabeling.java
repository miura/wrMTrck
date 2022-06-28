package wormTracker;

import java.awt.Color;
import java.awt.Font;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.ListIterator;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.process.ImageProcessor;

// makes a new stack with objects labeled with track nr
// optionally also displays centroid position
public class TrackNumberLabeling {

	ImagePlus imp;
	ArrayList<ArrayList<particle>> theTracks;
	HashMap<Integer, ArrayList<particle>> theParticles;
	double[][] bendCounter;

	public TrackNumberLabeling(ImagePlus imp, ArrayList<ArrayList<particle>> theTracks,
			HashMap<Integer, ArrayList<particle>> theParticles, double[][] bendCounter) {
		super();
		this.imp = imp;
		this.theTracks = theTracks;
		this.theParticles = theParticles;
		this.bendCounter = bendCounter;
		
	}
	
	public ImagePlus run() {
		IJ.showStatus("Generating movie with labels...");
			String strPart;
			Calibration cal = imp.getCalibration();
			double pixelWidth = cal.pixelWidth;
			double pixelHeight = cal.pixelHeight;
			ImageStack newstack = imp.createEmptyStack();
			int nFrames = imp.getStackSize();
			ImageStack stack = imp.getImageStack();
			int xHeight = newstack.getHeight();
			int yWidth = newstack.getWidth();

			for (int i = 0; i <= (nFrames - 1); i++) {
				if (IJ.escapePressed()) {
					IJ.beep();
					//done = true;
					return null;
				}
				int iFrame = i + 1;
				String strLine = "" + i;
				ImageProcessor ip = stack.getProcessor(iFrame);
				//newstack.addSlice(stack.getSliceLabel(iFrame), ip.crop());
				newstack.addSlice(stack.getSliceLabel(iFrame), ip.duplicate());
				ImageProcessor nip = newstack.getProcessor(iFrame);
				
				nip.setColor(Color.black);
//				Font f1 = Font.getFont(Font.FACE_SYSTEM, Font.STYLE_PLAIN, Font.SIZE_LARGE);
//				fontSize
				nip.setFont(new Font("SansSerif", Font.PLAIN, Parameters.fontSize));
				// nip.boldFont = true;
				// hack to only show tracks longerthan minTrackLength
				int trackNr = 0;
				int displayTrackNr = 0;
				for (ListIterator iT = theTracks.listIterator(); iT.hasNext();) {
					trackNr++;
					List bTrack = (ArrayList) iT.next();
					if (bTrack.size() >= Parameters.minTrackLength) {
						displayTrackNr++;
						for (ListIterator k = theParticles.get( i+ 1 ).listIterator(); k.hasNext();) {
							particle aParticle = (particle) k.next();
							if (aParticle.trackNr == trackNr) {
								nip.setColor(127);// Color.black);
								strPart = "" + displayTrackNr;
								if (Parameters.bShowPositions) {
									// if (bendType==0) strPart += "= "+(int)aParticle.x+","+(int)aParticle.y ; //
									// bend detection not enabled plot coordinates instead
									if (Parameters.bendType == 1)
										strPart += "= b" + bendCounter[displayTrackNr][i];
									if (Parameters.bendType > 1)
										strPart += "= b" + bendCounter[displayTrackNr][i] / 2;

									// we could plot a number of different infos for the particles tracked
									// {strPart+="="+/*(int)aParticle.angle+","*/+(int)aParticle.area+","+(int)aParticle.x+","+(int)aParticle.y;}
									// {strPart+="="+/*(int)aParticle.angle+","+*/(int)aParticle.area+","+bendCounter[displayTrackNr][i]/2;
									// };
								}
								// we could do someboundary testing here to place the labels better when we are
								// close to the edge
								nip.moveTo((int) (aParticle.x / pixelWidth + 0),
										doOffset((int) (aParticle.y / pixelHeight), yWidth, 5));
								// nip.moveTo(doOffset((int)aParticle.x,xHeight,5),doOffset((int)aParticle.y,yWidth,5)
								// );
								nip.drawString(strPart);
								if (Parameters.bendType == 1) {
									nip.setColor(Color.gray);
									double rad = Math.toRadians(aParticle.angle + 90);
									nip.moveTo((int) (aParticle.x / pixelWidth + 20 * Math.sin(rad)),
											(int) (aParticle.y / pixelHeight + 20 * Math.cos(rad)));
									nip.lineTo((int) (aParticle.x / pixelWidth - 20 * Math.sin(rad)),
											(int) (aParticle.y / pixelHeight - 20 * Math.cos(rad)));
								}
							}
						}
					}
				}
				IJ.showProgress((double) iFrame / nFrames);
			}
			ImagePlus nimp = new ImagePlus(imp.getTitle() + " labels", newstack);
			nimp.setCalibration(imp.getCalibration());
			
			//nimp.show();
			//imp.show();
			//nimp.updateAndDraw();
			return nimp;
		}		

	// Utility functions

	int doOffset(int center, int maxSize, int displacement) {
		if ((center - displacement) < 2 * displacement) {
			return (center + 4 * displacement);
		} else {
			return (center - displacement);
		}
	}
}
