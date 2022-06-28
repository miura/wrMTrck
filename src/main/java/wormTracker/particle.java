package wormTracker;


public class particle {

	float	x;
	float	y;
	float	sx;		// smoothed x-coordinate
	float	sy;		// smoothed y-coordinate
	int	z;
	float	area;	// inserted by JSP
	float	major;
	float	minor;
	float	ar;
	float	angle;
	float	circularity;
	float	dist;  // contains the distance traveled since last frame.
	float   perimeter; // length of the perimeter of a object
	int	trackNr;
	boolean inTrack=false;
	boolean flag=false;
	
	public particle() {
		
	}
	public particle(float x, float y, int z, 
			float area, float angle, float major, 
			float minor, float perimeter) {
		super();
		this.x = x;
		this.y = y;
		if (Parameters.bRoundCoord) {
			this.x = Math.round( x );
			this.y = Math.round( y );			
		}
		this.z = z;
		this.area = area;
		this.major = major;
		this.minor = minor;
		this.angle = angle;
		//this.circularity = circularity;
		this.perimeter = perimeter;
		if (minor > 0)
			this.ar = major / minor;
	}
	
	public void copy(particle source) {
		this.x=source.x;
		this.y=source.y;
		this.z=source.z;  // frame
		this.major=source.major;
		this.minor=source.minor;
		this.ar=source.ar;
		this.area= source.area;
		this.angle = source.angle;
		this.circularity = source.circularity;
		this.dist = source.dist;
		this.perimeter = source.perimeter;
		this.inTrack=source.inTrack;
		this.flag=source.flag;
	}

	public float distance (particle p) {
		return (float) Math.sqrt(sqr(this.x-p.x) + sqr(this.y-p.y));
	}

	public float sdistance (particle p) {
		return (float) Math.sqrt(sqr(this.sx-p.sx) + sqr(this.sy-p.sy));
	}
	
	double sqr(double n) {
		return n*n;
	}
}
