#include "Geom.h"
    
#include <sstream>
    
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



// ELEMENT RENDERING

LineString RenderLine(const Vec2& p0, const Vec2& p1)
{
    return { p0, p1 };
}

LineString RenderArc(const Vec2& pC, double radius, Direction direction, double th_Start, double th_End, int arcSegments)
{
    LineString linestring;
    // Clean up angles
    CleanAngles(th_Start, th_End, direction);
    // Calculate increment from n segments in 90 degrees
    double th_Incr = direction * (M_PI / 2.0) / arcSegments;
    // Calculate number of incrments for loop
    int nIncrements = floorf(fabsf((th_End - th_Start) / th_Incr));
    
    // from 'n == 1' because we have already calculated the first angle
    // to 'n == nIncremenets' to ensure last point is added
    for (int n = 0; n <= nIncrements; n++) {
        
        double th = (n == nIncrements) ? th_End : th_Start + n * th_Incr;
        // Calculate position from radius and angle
        Vec2 p = pC + Vec2(fabsf(radius) * sin(th), fabsf(radius) * cos(th));       
        
        // This prevents double inclution of point 
        if(!linestring.empty()) { 
            if(p == linestring.back()) { continue; }
        }
        
        //Add Line to output
        linestring.emplace_back(std::move(p));
    }
    return std::move(linestring);
}

LineString RenderArc(const Vec2& p0, const Vec2& p1, const Vec2& pC, MaxLib::Geom::Direction direction, int arcSegments)
{
    // get start and end points relative to the centre point
    Vec2 v_Start    = p0 - pC;
    Vec2 v_End      = p1 - pC;
    // get start and end angles
    double th_Start = atan2(v_Start.x, v_Start.y);
    double th_End   = atan2(v_End.x, v_End.y);
    double radius   = hypot(v_End.x, v_End.y);
    // draw arc between angles
    LineString arc = RenderArc(pC, radius, direction, th_Start, th_End, arcSegments);
    // adjust the front and back points to remove rounding errors
    if(arc.size() >= 2) {
        arc.front() = p0;
        arc.back() = p1;
    }
    return arc;
}

LineString RenderCircle(const Vec2& pC, double radius, int arcSegments) {
    // draw arc between angles
    LineString circle = RenderArc(pC, radius, Direction::CW, 0.0, 2.0 * M_PI, arcSegments);
    // make sure the first and last points match as 2PI produces a rounding error
    if(!circle.empty()) { circle.back() = circle.front(); }
    return std::move(circle);
}



 // This takes a std::string of 3 values seperated by commas (,) and will return a 3DPoint
// 4.000,0.000,0.000
Vec2 StringToVec2(const std::string& msg) 
{
    std::stringstream stream(msg);
    std::string segment;
    Vec2 p;
    
    getline(stream, segment, ',');
    p.x = stof(segment);
    getline(stream, segment);
    p.y = stof(segment);
    return p;
}
 
// This takes a std::string of 3 values seperated by commas (,) and will return a 3DPoint
// 4.000,0.000,0.000
Vec3 StringToVec3(const std::string& msg) 
{
    std::stringstream stream(msg);
    std::string segment;
    Vec3 p;
    
    std::getline(stream, segment, ',');
    p.x = stof(segment);
    std::getline(stream, segment, ',');
    p.y = stof(segment);
    std::getline(stream, segment);
    p.z = stof(segment);
    return p;
}

double RoundTo(double input, double roundto) {
    double x = input / roundto;
    return (round(x) * roundto);
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

// Calculates point perpendicular to line (p0->p1) at offset away from pStart
Vec2 PointPerpendicularToLine (const Vec2& p0, const Vec2& p1, double offset, const Vec2& pStart) 
{
    return pStart + (offset / hypot(p0.y-p1.y, p1.x-p0.x)) * Vec2(p0.y-p1.y, p1.x-p0.x);
}

Vec2 ArcCentre(const Vec2& p0, const Vec2& p1, double r, Direction direction)
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

// Finds the closest centre point to pC (at the same perpendicular distance as pC is from line (p0, p1))
Vec2 ArcCentre(const Vec2& p0, const Vec2& p1, const Vec2& pC) {
    // Get distance between line and point
    double d = DistanceBetween(p0, p1, pC);
    
    double th = Polar(p1 - p0).th;
    double thPerpendicular = CleanAngle(th + M_PI_2);
    Polar newCentre = Polar(d, thPerpendicular);
    
    Vec2 pMid = (p0 + p1) / 2.0;
    int flipSide = LeftOfLine(p0, p1, pC) ? -1 : 1;
    return pMid + newCentre.Cartesian() * flipSide;
}

// Finds the centre point of an arc (given p0 & p1) which is tangent to line (pt -> p0)
// Will return empty if:
//      Any of the points are the same
//      The 3 points are on the same line
std::optional<Vec2> ArcCentreFromTangentLine(const Vec2& pt, const Vec2& p0, const Vec2& p1) 
{
    // We calculate the line perpendicular to line (l0 -> p0) at p0
    // We calculate the line perpendicular to line (p0 -> p1) at the midpoint of line
    // the centre point is the intersection point between these 2 lines.
    // l1 is equiv. to p0
    
    // error points aren't different
    if(pt == p0 || p0 == p1 || p1 == pt)    { return {}; }
    
    Vec2 line1 = p0 - pt;
    Vec2 line2 = p1 - p0;
    
    // inverted gradients to give perpendicular lines
    double m1 = -line1.x / line1.y;
    double m2 = -line2.x / line2.y;
    
    // error if gradients are the same (includes check for 0 == -0 and inf == -inf)
    if(m1 == m2 || (std::isinf(m1) && std::isinf(m2))) { return {}; }
    
    // at points
    Vec2 a = p0;                // end of line1 & beginning of arc / line2
    Vec2 b = (p1 + p0) / 2.0;   // mid point of line2
    
    double c1 = a.y - a.x * m1;
    double c2 = b.y - b.x * m2;
        
    // Calculate x/y coordinates for centre point
    Vec2 pC;
    
    // line1 is horizontal
    if(line1.y == 0.0)      { pC.x = p0.x; }
    // line2 is horizontal
    else if(line2.y == 0.0) { pC.x = b.x; }
    // calculate x
    else                    { pC.x = (c2 - c1) / (m1 - m2); }
    // line1 is vertical
    if(line1.x == 0.0)      { pC.y = p0.y; }
    // line 2 is vertical
    else if(line2.x == 0.0) { pC.y = b.y; }
    // calculate y = m1 * x + c1, (use y = m2 * x + c2 if line is horizontal)
    else                    { pC.y = (line1.y == 0.0) ? pC.x*m2 + c2 : pC.x*m1 + c1; }
    
    // return centre point
    return pC;
}

// Finds the end point on a line which is tangent to an arc (p0(unused), p1, pC) with a length (d)
Vec2 ArcTangentLine(const Vec2& pC, const Vec2& p, Direction direction, double d) {
    
    // Angle of pC to p1
    double th = Polar(p - pC).th;
    // Angle of tangent line
    double thPerpendicular = CleanAngle(th + M_PI_2);
    
    Polar tangentLine = Polar(d, thPerpendicular);
    
    return p + tangentLine.Cartesian() * direction;
}

// returns tangent point of circle to point p0, (circle has centre pC, and radius r)     side: 1 is left, -1 is right      TODO: ***** CHECK DIRECITON IS CORRECT******
Vec2 CircleTangentPoint(const Vec2& p0, const Vec2& pC, double r, int side)
{
    Vec2 line = p0 - pC;
    Polar pol = Polar(line);
    // length of line
    double L = pol.r;
    // angle of line
    double th = asin(r / L);
    return p0 + Polar(sqrt(L*L-r*r), (pol.th - (double)side*th)).Cartesian();        
}


// returns true if point ios
bool LeftOfLine(const Vec2& p1, const Vec2& p2, const Vec2& pt)
{
    return ((p2.x-p1.x)*(pt.y-p1.y) - (p2.y-p1.y)*(pt.x-p1.x)) > 0.0;
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
