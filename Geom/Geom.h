#pragma once
/*
 * Name: 
 *    Geom
 * Description: 
 *    A math and geometry library, part of the Utils library
 */

#include <iostream>
#include <cmath>
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
    Polar(Vec2 p);
    Vec2 Cartesian();
};
// Overload polar operators for cout
static inline std::ostream& operator<<(std::ostream& os, const Polar& pol) { os << "(" << pol.r << ", " << Degrees(pol.th) << "degs)"; return os; }


// ************************************************************************************** //
// ********************************* Geom Functions ************************************* //

// returns the sign (1. 0 or -1) of a number
// 0 can be set to 1 or -1 with zeroValue
template <typename T> 
int                          Sign(T val, int zeroValue = 0) { int a = (zeroValue == 1) ? (T(0) <= val) : (T(0) < val); int b = (zeroValue == -1) ? (val <= T(0)) : (val < T(0)); return a-b; }

// Rounds a number to the nearest decimal place (i.e. RoundTo(1.166, 0.01) = 1.17)
double                       RoundTo(double input, double dp);

inline Vec2                  Abs(const Vec2& p)   { return { fabs(p.x), fabs(p.y) }; }
inline Vec3                  Abs(const Vec3& p)   { return { fabs(p.x), fabs(p.y), fabs(p.z) }; }
inline double                Hypot(const Vec2& p) { return sqrt(p.x*p.x + p.y*p.y); }
inline double                Hypot(const Vec3& p) { return sqrt(p.x*p.x + p.y*p.y + p.z*p.z); }


// Returns Angle between 0 - 2PI
double                       CleanAngle(double angle);
// Modifies start and end angles to be in correct order based on direction 
// Makes both >= 0     
// Output will produce <= 4*PI difference between angles    (anticlockwise curve could be start:710degs to end:359degs
// If both angles are identical, they will be set to 2*PI out of phase
// Direction:  1 CW   -1 CCW
void                         CleanAngles(double& startAngle, double& endAngle, Direction direction);
    
// returns angle from positive x axis in CW direction based on centre point and end point
std::optional<double>        AngleBetween(Vec2 centre, Vec2 end);

// calculates angle between 3 points        p1 is start, p2 is centre, p3 is end
std::optional<double>        AngleBetween(Vec2 p1, Vec2 p2, Vec2 p3, Direction direction = Direction::CW);


// calculates centre from radius, start & end points (-r will return the second possible arc)
Vec2                         ArcCentreFromRadius(const Vec2& p0, const Vec2& p1, double r, Direction direction);

 
// Calculate determinant of matrix:  [a b]
//                                   [c d] 
double                       Determinant(double a, double b, double c, double d);
// returns true if point ios
bool                         LeftOfLine(const Vec2& p1, const Vec2& p2, const Vec2& pt);
// returns tangent point of circle to point p0, (circle has centre pC, and radius r)     side: 1 is left, -1 is right
Vec2                         TangentOfCircle(const Vec2& p0, const Vec2& pC, double r, int side);


// {} if no intersect,  p if intersect point
typedef std::optional<Vec2> Intersect;
// ({}, {}) if no intersect,  (p, {}) if 1 intersect point(tangent)   and  (p1, p2) if 2 intersect points
typedef std::pair<Intersect, Intersect> IntersectPair;

// calculates the intersection points between 2 circles and passes back in the 2 return pointers
IntersectPair                IntersectTwoCircles(const Vec2& c1, double r1, const Vec2& c2, double r2);
// calculates the intersection points between a line and a circle and passes back in the 2 return pointers
IntersectPair                IntersectLineCircle(const Vec2& p1, const Vec2& p2, const Vec2& c, double r);
// Calculate intersection of two lines.
Intersect                    IntersectLines(const Vec2& p1, const Vec2& p2, const Vec2& p3, const Vec2& p4);
// returns whether there is an intersect or not (faster alternative to IntersectLines() but does not return location)
bool                         IntersectLinesFast(const Vec2& p1, const Vec2& p2, const Vec2& p3, const Vec2& p4);


} // end namespace Geom
} // end namespace MaxLib
