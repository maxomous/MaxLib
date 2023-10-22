#pragma once
/*
 * Name: 
 *    Geom
 * Description: 
 *    A math and geometry library, part of the Utils library
 */

#include <iostream>
#include <cmath>
#include <vector>
#include <optional>
#include <tuple>


namespace MaxLib {
namespace Geom {

// ************************************************************************************** //
// ********************************* Basic Functions ************************************ //


// Convert Degrees to Radians
inline double                Radians(double deg) { return (deg) * M_PI / 180.0; }
// Convert Radians to Degrees
inline double                Degrees(double rad) { return (rad) * 180.0 / M_PI; }
// Sin Cos & Tan of angle in Degrees                
inline double                Sin(double deg) { return sin(Radians(deg)); }
inline double                Cos(double deg) { return cos(Radians(deg)); }
inline double                Tan(double deg) { return tan(Radians(deg)); }


// ************************************************************************************** //
// *********************************** Geom Types *************************************** //

enum Direction { CW = 1, CCW = -1 };

// Vec2 (x, y)
class Vec2 {
public:
    // ensures the inhertited destructor is called
    virtual ~Vec2() {}
    // Variables
    double x;
    double y; 
    // Constructors
    Vec2()                   { x = 0; y = 0; }
    Vec2(double X, double Y) { x = X; y = Y; }
};
// Vec2 overload operators (+ - * / += -= == != <<)
static inline Vec2  operator+(const Vec2& a, const Vec2& b)  { return Vec2(a.x + b.x, a.y + b.y); }
static inline Vec2  operator-(const Vec2& a, const Vec2& b)  { return Vec2(a.x - b.x, a.y - b.y); }
static inline Vec2  operator*(const Vec2& a, const Vec2& b)  { return Vec2(a.x * b.x, a.y * b.y); }
static inline Vec2  operator/(const Vec2& a, const Vec2& b)  { return Vec2(a.x / b.x, a.y / b.y); }
static inline Vec2  operator+(const Vec2& a, const double b) { return Vec2(a.x + b, a.y + b); }
static inline Vec2  operator-(const Vec2& a, const double b) { return Vec2(a.x - b, a.y - b); }
static inline Vec2  operator*(const Vec2& a, const double b) { return Vec2(a.x * b, a.y * b); } 
static inline Vec2  operator/(const Vec2& a, const double b) { return Vec2(a.x / b, a.y / b); }
static inline Vec2  operator+(const double a, const Vec2& b) { return Vec2(a + b.x, a + b.y); }
static inline Vec2  operator-(const double a, const Vec2& b) { return Vec2(a - b.x, a - b.y); }
static inline Vec2  operator*(const double a, const Vec2& b) { return Vec2(a * b.x, a * b.y); } 
static inline Vec2  operator/(const double a, const Vec2& b) { return Vec2(a / b.x, a / b.y); }
static inline Vec2& operator+=(Vec2& a, const Vec2& b)       { a.x += b.x; a.y += b.y; return a;}
static inline Vec2& operator-=(Vec2& a, const Vec2& b)       { a.x -= b.x; a.y -= b.y; return a;}
static inline Vec2& operator+=(Vec2& a, const double b)      { a.x += b; a.y += b; return a;}
static inline Vec2& operator-=(Vec2& a, const double b)      { a.x -= b; a.y -= b; return a;}
static inline bool  operator==(const Vec2& a, const Vec2& b) { return (a.x == b.x && a.y == b.y); }
static inline bool  operator!=(const Vec2& a, const Vec2& b) { return (a.x != b.x || a.y != b.y); }
static inline bool  operator<(const Vec2& a, const Vec2& b)  { return (a.x < b.x && a.y < b.y); }
static inline bool  operator>(const Vec2& a, const Vec2& b)  { return (a.x > b.x && a.y > b.y); }
static inline bool  operator>=(const Vec2& a, const Vec2& b) { return (a.x >= b.x && a.y >= b.y); }
static inline bool  operator<=(const Vec2& a, const Vec2& b) { return (a.x <= b.x && a.y <= b.y); }
static inline std::ostream& operator<<(std::ostream& os, const Vec2& p) { os << "(" << p.x << ", " << p.y << ")"; return os; }

// Vec3 (x, y, x)
class Vec3 : public Vec2 {
public:
    // Variables
    double z;
    // Constructors
    Vec3()                             { x = 0; y = 0; z = 0; }
    Vec3(double X, double Y, double Z) { x = X; y = Y; z = Z; }
};
// Vec3 overload operators (+ - * / += -= == != <<)
static inline Vec3  operator+(const Vec3& a, const Vec3& b)  { return Vec3(a.x + b.x, a.y + b.y, a.z + b.z); }
static inline Vec3  operator-(const Vec3& a, const Vec3& b)  { return Vec3(a.x - b.x, a.y - b.y, a.z - b.z); }
static inline Vec3  operator*(const Vec3& a, const Vec3& b)  { return Vec3(a.x * b.x, a.y * b.y, a.z * b.z); }
static inline Vec3  operator/(const Vec3& a, const Vec3& b)  { return Vec3(a.x / b.x, a.y / b.y, a.z / b.z); }
static inline Vec3  operator+(const Vec3& a, const double b) { return Vec3(a.x + b, a.y + b, a.z + b); }
static inline Vec3  operator-(const Vec3& a, const double b) { return Vec3(a.x - b, a.y - b, a.z - b); }
static inline Vec3  operator*(const Vec3& a, const double b) { return Vec3(a.x * b, a.y * b, a.z * b); }
static inline Vec3  operator/(const Vec3& a, const double b) { return Vec3(a.x / b, a.y / b, a.z / b); }
static inline Vec3& operator+=(Vec3& a, const Vec3& b)       { a.x += b.x; a.y += b.y; a.z += b.z; return a;}
static inline Vec3& operator-=(Vec3& a, const Vec3& b)       { a.x -= b.x; a.y -= b.y; a.z -= b.z; return a;}
static inline Vec3& operator+=(Vec3& a, const double b)      { a.x += b; a.y += b; a.z += b; return a;}
static inline Vec3& operator-=(Vec3& a, const double b)      { a.x -= b; a.y -= b; a.z -= b; return a;}
static inline bool  operator==(const Vec3& a, const Vec3& b) { return (a.x == b.x && a.y == b.y && a.z == b.z); }
static inline bool  operator!=(const Vec3& a, const Vec3& b) { return (a.x != b.x || a.y != b.y || a.z != b.z); }
static inline bool  operator<(const Vec3& a, const Vec3& b)  { return (a.x < b.x && a.y < b.y && a.z < b.z); }
static inline bool  operator>(const Vec3& a, const Vec3& b)  { return (a.x > b.x && a.y > b.y && a.z > b.z); }
static inline bool  operator<=(const Vec3& a, const Vec3& b) { return (a.x <= b.x && a.y <= b.y && a.z <= b.z); }
static inline bool  operator>=(const Vec3& a, const Vec3& b) { return (a.x >= b.x && a.y >= b.y && a.z >= b.z); }
static inline std::ostream& operator<<(std::ostream& os, const Vec3& p) { os << "(" << p.x << ", " << p.y << ", " << p.z << ")"; return os; }

// Polar class
// 0 degs is at 12 o'clock
// Positive is clockwise
struct Polar {
    // Variables
    double r;    // Length
    double th;     // Angle
    // Constructor
    Polar(double R = 0.0, double Th = 0.0);
    Polar(const Vec2& p);
    Vec2 Cartesian();
};
// Overload polar operators for cout
static inline std::ostream& operator<<(std::ostream& os, const Polar& pol) { os << "(" << pol.r << ", " << Degrees(pol.th) << "degs)"; return os; }


// 2D Containers
// A container of (x, y) coords (Different names for ease of readability)
typedef std::vector<Vec2>  Geometry;
typedef std::vector<Vec2>  Points;
typedef std::vector<Vec2>  LineString; // A linestring is treated as a polypon when its first and last points are identical in some functions
// Polygon contains a shell and holes
struct Polygon
{
    Polygon() {}
    Polygon(LineString&& geometry) : shell(std::move(geometry)) {}
    Polygon(LineString&& Shell, std::vector<LineString>&& Holes) : shell(std::move(Shell)), holes(std::move(Holes)) {}
    LineString shell;
    std::vector<LineString> holes;
};
// A collection of any geometry type. Often used when we dont know what type of geometry will be returned
struct GeometryCollection
{
    std::vector<Vec2>       points;
    std::vector<LineString> lineStrings;
    std::vector<Polygon>    polygons;
    // Copy point, linestring or polygon to collection
    void AppendGeometry(const Vec2& pts)                { points.emplace_back(pts); }
    void AppendGeometry(const LineString& linestring)   { lineStrings.emplace_back(linestring); }
    void AppendGeometry(const Polygon& polygon)         { polygons.emplace_back(polygon); }
};

// ************************************************************************************** //
// ********************************* Geom Functions ************************************* //

// LineString Rendering
// Render element to linestring
LineString                  RenderLine(const Vec2& p0, const Vec2& p1); 
LineString                  RenderArc(const Vec2& pC, double radius, Direction direction, double th_Start, double th_End, double arcTolerance = 0.001);
LineString                  RenderArc(const Vec2& p0, const Vec2& p1, const Vec2& pC, Direction direction, double arcTolerance = 0.001);
LineString                  RenderCircle(const Vec2& pC, double radius, double arcTolerance = 0.001);
// Returns number of segments of a 90 deg arc which has a maximum deviation tolerance for a given radius
int                         ArcSegments(double radius, double tolerance);


// This takes a std::string of 3 values seperated by commas (,) and will return a 2D Point
// e.g. "4.000,0.000,0.000"
Vec2                        StringToVec2(const std::string& msg);
// This takes a std::string of 3 values seperated by commas (,) and will return a 3D Point
// e.g. "4.000,0.000,0.000"
Vec3                        StringToVec3(const std::string& msg);

// returns the sign (1. 0 or -1) of a number
// 0 can be set to 1 or -1 with zeroValue
template <typename T> 
int                         Sign(T val, int zeroValue = 0) { int a = (zeroValue == 1) ? (T(0) <= val) : (T(0) < val); int b = (zeroValue == -1) ? (val <= T(0)) : (val < T(0)); return a-b; }

// Rounds a number to the nearest decimal place (i.e. RoundTo(1.166, 0.01) = 1.17)
double                      RoundTo(double input, double dp);


// Calculates the dot product of 2x Vec2's
double                      DotProduct(const Vec2& v1, const Vec2& v2);

// Calculate determinant of matrix:  [a b]
//                                   [c d] 
double                      Determinant(double a, double b, double c, double d);

inline Vec2                 Abs(const Vec2& p)                                 { return { fabs(p.x), fabs(p.y) }; }
inline Vec3                 Abs(const Vec3& p)                                 { return { fabs(p.x), fabs(p.y), fabs(p.z) }; }
inline double               Hypot(const Vec2& p)                               { return sqrt(p.x*p.x + p.y*p.y); }
inline double               Hypot(const Vec3& p)                               { return sqrt(p.x*p.x + p.y*p.y + p.z*p.z); }


inline double               DistanceBetween(const Vec2& p0, const Vec2& p1)    { return Hypot(p1 - p0); }
inline double               DistanceBetween(const Vec3& p0, const Vec3& p1)    { return Hypot(p1 - p0); }
// Returns the minimum distance between Line l and Point p
double                      DistanceBetween(const Vec2& l0, const Vec2& l1, const Vec2& p);

// Returns Angle between 0 - 2PI
double                      CleanAngle(double angle);
// Modifies start and end angles to be in correct order based on direction 
// Makes both >= 0     
// Output will produce <= 4*PI difference between angles    (anticlockwise curve could be start:710degs to end:359degs
// If both angles are identical, they will be set to 2*PI out of phase
// Direction:  1 CW   -1 CCW
void                        CleanAngles(double& startAngle, double& endAngle, Direction direction);
    
// returns angle from positive x axis in CW direction based on centre point and end point
std::optional<double>       AngleBetween(const Vec2&  centre, const Vec2&  end);

// calculates angle between 3 points        p1 is start, p2 is centre, p3 is end
std::optional<double>       AngleBetween(const Vec2&  p1, const Vec2& p2, const Vec2& p3, Direction direction = Direction::CW);

// Calculates point perpendicular to line (p0->p1) at offset away from pStart
Vec2                        PointPerpendicularToLine(const Vec2& p0, const Vec2& p1, double offset, const Vec2& pStart);

// calculates centre from radius, start & end points (-r will return the second possible arc)
Vec2                        ArcCentre(const Vec2& p0, const Vec2& p1, double r, Direction direction);

// Finds the closest centre point to pC (at the same perpendicular distance as pC is from line (p0, p1))
Vec2                        ArcCentre(const Vec2& p0, const Vec2& p1, const Vec2& pC);

// Finds the centre point of an arc (given p0 & p1) which is tangent to line (pt -> p0)
// Will return empty if:
//      Any of the points are the same
//      The 3 points are on the same line
std::optional<Vec2>         ArcCentreFromTangentLine(const Vec2& l0, const Vec2& p0, const Vec2& p1);

// Finds the end point on a line of length d which is tangent to an arc with centre pC and tangent point p, the lines travels in direction 
Vec2                        ArcTangentLine(const Vec2& pC, const Vec2& p, Direction direction, double d);


// returns tangent point of circle to point p0, (circle has centre pC, and radius r)     side: 1 is left, -1 is right
Vec2                        CircleTangentPoint(const Vec2& p0, const Vec2& pC, double r, int side);


// returns true if point is left of line
bool                        LeftOfLine(const Vec2& p1, const Vec2& p2, const Vec2& pt);

// {} if no intersect,  p if intersect point
typedef std::optional<Vec2> Intersect;
// ({}, {}) if no intersect,  (p, {}) if 1 intersect point(tangent)   and  (p1, p2) if 2 intersect points
typedef std::pair<Intersect, Intersect> IntersectPair;

// calculates the intersection points between 2 circles and passes back in the 2 return pointers
IntersectPair               IntersectTwoCircles(const Vec2& c1, double r1, const Vec2& c2, double r2);
// calculates the intersection points between a line and a circle and passes back in the 2 return pointers
IntersectPair               IntersectLineCircle(const Vec2& p1, const Vec2& p2, const Vec2& c, double r);
// Calculate intersection of two lines.
Intersect                   IntersectLines(const Vec2& p1, const Vec2& p2, const Vec2& p3, const Vec2& p4);
// returns whether there is an intersect or not (faster alternative to IntersectLines() but does not return location)
bool                        IntersectLinesFast(const Vec2& p1, const Vec2& p2, const Vec2& p3, const Vec2& p4);


} // end namespace Geom
} // end namespace MaxLib
