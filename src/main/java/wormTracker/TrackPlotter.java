package wormTracker;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;

import ij.ImagePlus;
import ij.measure.Calibration;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;

public class TrackPlotter {
	ImagePlus imp;
	List theTracks;

	
	public TrackPlotter(ImagePlus imp, List theTracks) {
		super();
		this.imp = imp;
		this.theTracks = theTracks;
	}

	ImagePlus run(){
		Calibration cal = imp.getCalibration();
		double pixelWidth = cal.pixelWidth;
		double pixelHeight = cal.pixelHeight;
		ImageProcessor ip = new ByteProcessor(imp.getWidth(), imp.getHeight());
		ip.setColor(Color.white);
		ip.fill();
		int trackCount = 0;
		int color;
		for (ListIterator iT = theTracks.listIterator(); iT.hasNext();) {
			List bTrack = (ArrayList) iT.next();
			if (bTrack.size() >= Parameters.minTrackLength) {
				trackCount++;
				ListIterator jT = bTrack.listIterator();
				particle oldParticle = (particle) jT.next();
				for (; jT.hasNext();) {
					particle newParticle = (particle) jT.next();
					color = Math.min(trackCount + 1, 254);
					ip.setValue(color);
					ip.moveTo((int) (oldParticle.x / pixelWidth), (int) (oldParticle.y / pixelHeight));
					ip.lineTo((int) (newParticle.x / pixelWidth), (int) (newParticle.y / pixelHeight));
					oldParticle = newParticle;
				}
			}
		}
		ImagePlus newimp = new ImagePlus("Paths", ip);
		newimp.show();
		return newimp;
	}
}
