#include "Geom.h"
    
namespace MaxLib {
namespace Geom {
    
Polar::Polar(double R, double Th)  
{ 
    r = R; th = Th; 
}

Polar::Polar(const Vec2& p) 
{ 
    r = hypot(p.x, p.y); 
    th = CleanAngle(atan2(p.x, p.y)); 
}

Vec2 Polar::Cartesian() 
{ 
    return Vec2 (r*cos(th - M_PI_2),  r*-sin(th - M_PI_2)); 
}

double RoundTo(double input, double roundto) {
    double x = input / roundto;
    return (round(x) * roundto);
}

// Returns the minimum distance between Line l and Point p
double DistanceBetween(const Vec2& l0, const Vec2& l1, const Vec2& p) 
{    
    const double l2 = DistanceBetween(l0, l1) * DistanceBetween(l0, l1); // i.e. |l1-l0|^2 -  avoid a sqrt
    if (l2 == 0.0) return DistanceBetween(p, l0);   // l0 == l1 case
    // Consider the line extending the segment, parameterized as l0 + t (l1 - l0).
    // We find projection of point p onto the line. 
    // It falls where t = [(p-l0) . (l1-l0)] / |l1-l0|^2
    // We clamp t from [0,1] to handle points outside the segment vw.
    const double t = std::max(0.0, std::min(1.0, DotProduct(p - l0, l1 - l0) / l2));
    const Vec2 projection = l0 + (l1 - l0) * t;  // Projection falls on the segment
    return DistanceBetween(p, projection);
} 

// Returns Angle between 0 - 2PI
double CleanAngle(double angle)
{
    double th = angle;
    while ( (th >= 2*M_PI)  ||  (th < 0) ) 
    {
        if (th < 0) 
            th += 2*M_PI;
        else if (th >= 2*M_PI) 
            th -= 2*M_PI;
    }
    return th;
}

void CleanAngles(double& startAngle, double& endAngle, Direction direction)
{
    // swap start and end for anticlockwise
    double start = (direction == Direction::CW) ? startAngle : endAngle;
    double end = (direction == Direction::CW) ? endAngle : startAngle;
    
    // if values are equal, make out of phase by 2M_PI
    if (start == end) {  // special condition required due to rounding error -    after adding 2M_PI to end, (end - start) !>= 2*M_PI)
        end += 2*M_PI;
    }
    else { // make start and end within 360degs
        while ( ((end - start) > 2*M_PI)  ||  ((end - start) < 0) )
        {
            if (((end - start) < 0))
                end += 2*M_PI;
            else if ((end - start) > 2*M_PI) 
                end -= 2*M_PI;
        }
    }
    // phase shift both by -2*M_PI until one is less than zero (for vals greater than 4M_PI)
    while((start>0) && (end>0)) {
        start -= 2*M_PI;
        end -= 2*M_PI;
    }
    // phase shift both by +2*M_PI so they are both positive (for values < 0)
    while((start<0) || (end<0)) {
        start += 2*M_PI;
        end += 2*M_PI;
    }
    // swap and return start and end
    startAngle = (direction == Direction::CW) ? start : end;
    endAngle = (direction == Direction::CW) ? end : start;
    
}
 

// returns angle from positive x axis in CW direction based on centre point and end point
std::optional<double> AngleBetween(const Vec2& p0, const Vec2& p1)
{
    // error if points are the same
    if( p0 == p1 ) { return {}; }
    return Polar(p1-p0).th;
}

// calculates angle between 3 points        p1 is start, p2 is centre, p3 is end
std::optional<double> AngleBetween(const Vec2& p1, const Vec2& p2, const Vec2& p3, Direction direction)
{
    // error if any points are the same
    if( p1 == p2 || p2 == p3 || p3 == p1 ) { return {}; }
    // get angles of lines (p1 -> p2) & (p2 -> p3)
    double startAngle = Polar(p1-p2).th;
    double endAngle = Polar(p3-p2).th;
    CleanAngles(startAngle, endAngle, direction);
    // clockwise / anticlockwise
    double output = (direction == Direction::CW) ? (endAngle - startAngle) : (startAngle - endAngle); 
    return CleanAngle(output);
}

Vec2 ArcCentreFromRadius(const Vec2& p0, const Vec2& p1, double r, Direction direction)
{            
    Vec2 dif = p1 - p0;
    // midpoint of start to end    
    Vec2 pMid = (p0 + p1) / 2.0;
    
    // length between start and end points
    double     L = sqrt(dif.x * dif.x + dif.y * dif.y);
    
    //    angle between x axis & line from start to end
    double theta_G = fabs(atan(dif.y / dif.x));
    double h = sqrt(std::max(0.0, r*r - (L/2.0)*(L/2.0))); 
    
    h = direction * h;
    // 2nd version of the curve (when the centrepoint is past the midway line between start and end) 
    if(r < 0.0) { h = -h; }
    // prevent flipping of centrepoint when p0 & p1 are horizontal & vertical
    if(p0.y == p1.y && p1.x > p0.x) { h = -h; }
    if(p0.x == p1.x && p1.y < p0.y) { h = -h; }
    
    Vec2 pCentre;
    Vec2 invert = { ((p0.y > p1.y) ? -1.0 : 1.0), ((p1.x > p0.x) ? -1.0 : 1.0) };
        
    // if start to end is vertical
    if(dif.x == 0.0) {
        pCentre.x = p0.x + invert.y * h;
        pCentre.y = pMid.y;
    }    
    // if start to end is horizontal
    else if(dif.y == 0.0) { 
        pCentre.x = pMid.x;
        pCentre.y = p0.y + invert.x * h;
    }
    else  {        
        Vec2 hyp = { h*sin(theta_G), h*cos(theta_G) };
        pCentre = pMid + invert * hyp;
    }


    return pCentre;
}

double DotProduct(const Vec2& v1, const Vec2& v2)
{
    return v1.x * v2.x + v1.y * v2.y;
}

/** Calculate determinant of matrix:
    [a b]
    [c d] */
double Determinant(double a, double b, double c, double d)
{
    return a*d - b*c;
} 


// returns true if point ios
bool LeftOfLine(const Vec2& p1, const Vec2& p2, const Vec2& pt)
{
    return ((p2.x-p1.x)*(pt.y-p1.y) - (p2.y-p1.y)*(pt.x-p1.x)) > 0.0;
}

// returns tangent point of circle to point p0, (circle has centre pC, and radius r)     side: 1 is left, -1 is right      TODO: ***** CHECK DIRECITON IS CORRECT******
Vec2 TangentOfCircle(const Vec2& p0, const Vec2& pC, double r, int side)
{
    Vec2 line = p0 - pC;
    Polar pol = Polar(line);
    // length of line
    double L = pol.r;
    // angle of line
    double th = asin(r / L);
    return p0 + Polar(sqrt(L*L-r*r), (pol.th - (double)side*th)).Cartesian();        
}


// calculates the intersection points between 2 circles and passes back in the 2 return pointers
// return 0 if no intersect,   1 if 1 intersect point(tangent)   and 2 if 2 intersect points
IntersectPair IntersectTwoCircles(const Vec2& c1, double r1, const Vec2& c2, double r2)
{
    auto err = std::make_pair(std::nullopt, std::nullopt);
    // error if centre is s
    if(c1.x == c2.x && c1.y == c2.y) { return err; }
    double xdif = c2.x - c1.x;
    double ydif = c2.y - c1.y;
    double L = hypot(xdif, ydif);

    if ((L <= r1 + r2) && L >= fabs(r2 - r1)) {
        double ex = (c2.x - c1.x) / L;
        double ey = (c2.y - c1.y) / L;
        double x = (r1*r1 - r2*r2 + L*L) / (2*L);
        double y = sqrt(r1*r1 - x*x);
        Vec2 p1, p2;
        p1.x = c1.x + x*ex - y*ey;
        p1.y = c1.y + x*ey + y*ex;
        p2.x = c1.x + x*ex + y*ey;
        p2.y = c1.y + x*ey - y*ex;
        if(p1 == p2) { return std::make_pair(p1, std::nullopt); }
        else         { return std::make_pair(p1, p2); }
    }
    // No Intersection, far outside or one circle within the other
    else { return err; }
}

// calculates the intersection points between a line and a circle and passes back in the 2 return pointers
// return 0 if no intersect,   1 if 1 intersect point(tangent)   and 2 if 2 intersect points
IntersectPair IntersectLineCircle(const Vec2& p1, const Vec2& p2, const Vec2& c, double r)
{
    double dx = p2.x - p1.x;
    double dy = p2.y - p1.y;
    double L = hypot(dx, dy);
    double D = Determinant(p1.x-c.x, p1.y-c.y, p2.x-c.x, p2.y-c.y);
    double q = sqrt(r*r*L*L - D*D);
    
    double disc = r*r*L*L - D*D;
    // no intersection
    if(disc < 0.0) { return std::make_pair(std::nullopt, std::nullopt); }
    // calculate intersection point
    Vec2 pIntersect1 = Vec2(c.x + ((D*dy  + Sign(dy)*dx*q) / (L*L)),   c.y + ((-D*dx + fabs(dy)*q) / (L*L)));
    Vec2 pIntersect2 = Vec2(c.x + ((D*dy  - Sign(dy)*dx*q) / (L*L)),   c.y + ((-D*dx - fabs(dy)*q) / (L*L)));
    // 1 intersection (tangent)
    if (disc == 0.0) { 
        return std::make_pair(pIntersect1, std::nullopt); 
    }
    // 2 intersections
    // else if( disc > 0.0 )
    return std::make_pair(pIntersect1, pIntersect2); 
}

// returns whether there is an intersect or not (faster alternative to IntersectLines however does not return location
bool IntersectLinesFast(const Vec2& p1, const Vec2& p2, const Vec2& p3, const Vec2& p4)
{
    double test1 = (p2.x - p1.x)*(p3.y - p2.y) - (p2.y - p1.y)*(p3.x - p2.x);
    double test2 = (p2.x - p1.x)*(p4.y - p2.y) - (p2.y - p1.y)*(p4.x - p2.x);
    return Sign(test1) != Sign(test2);
}
// returns 0 on success
///Calculate intersection of two lines.
Intersect IntersectLines(const Vec2& p1, const Vec2& p2, const Vec2& p3, const Vec2& p4)
{
    // http://mathworld.wolfram.com/Line-LineIntersection.html
    double detL1 = Determinant(p1.x, p1.y, p2.x, p2.y);
    double detL2 = Determinant(p3.x, p3.y, p4.x, p4.y);
    double x1mx2 = p1.x - p2.x;
    double x3mx4 = p3.x - p4.x;
    double y1my2 = p1.y - p2.y;
    double y3my4 = p3.y - p4.y;

    double xnom = Determinant(detL1, x1mx2, detL2, x3mx4);
    double ynom = Determinant(detL1, y1my2, detL2, y3my4);
    double denom = Determinant(x1mx2, y1my2, x3mx4, y3my4);
    
    // Intersection lines do not cross
    if(denom == 0.0) { return {}; }

    double ixOut = xnom / denom;    
    double iyOut = ynom / denom;
    // intersection is not finite, probably a numerical issue
    if(!std::isfinite(ixOut) || !std::isfinite(iyOut)) { return {}; }
        
    return Vec2(ixOut, iyOut);
}

} // end namespace Geom
} // end namespace Maxomous
