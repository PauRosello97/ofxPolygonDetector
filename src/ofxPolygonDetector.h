
#pragma once

#include <algorithm>
#include <cstdint>
#include <string>
#include <cmath>
#include "ofMain.h"
#include <map>
#include <unordered_map>
#include <set>
#include <list>
#include <vector>

static const float minPointDiff = 0.001;
static const float minPointDiffSq = minPointDiff * minPointDiff;

struct ofxPolyLine;
struct ofxPolygonDetector;

struct PointType
{
    float x, y;

    PointType() {}
    PointType(float a, float b) : x(a), y(b) {}
    PointType(float a) : x(a), y(a) {}
    PointType(float a, float b, float c) : x(a), y(b) {}

    PointType& sub(const PointType& v);
    PointType& add(const PointType& v);
    PointType& mul(const float& v);
    PointType& div(const float& v);
    float squaredlen() const;
    float squaredist(const PointType& v) const;
    static void line(const PointType& p, const PointType& q, float& a, float& b, float& c);
    static float lineDist(float a, float b, float c, const PointType& p);
    static float lineDist(const PointType& la, const PointType& lb, const PointType& p);
};

using CycleSet = set<int>;

bool doIntersect(const PointType& p1, const PointType& q1, const PointType& p2, const PointType& q2);
int orientation(const PointType& p, const PointType& q, const PointType& r);
bool onSegment(const PointType& p, const PointType& q, const PointType& r);
bool collinearVecs(const PointType& p, const PointType& q, const PointType& r);
bool between(const PointType& p, const PointType& a, const PointType& b);
bool pointsDiffer(const PointType& a, const PointType& b, bool aprox = true);
bool overlap(const ofxPolyLine& l1, const ofxPolyLine& l2);
int simplifiedLine(const ofxPolyLine& line_1, const ofxPolyLine& line_2, ofxPolyLine& ret);
int iComparePointOrder(const PointType& p1, const PointType& p2);
bool bComparePointOrder(const PointType& p1, const PointType& p2);


enum class RmLinesType : int
{
    TakenTwice = 0, // rm lines taken twice
    Collinear, // dissolve collinear lines
    NoPointNeigh, // rm lines with one point not having neighbors
    PointConsumed, // point has 2 lines connected all taken
};

struct PolyCycle
{
    CycleSet idx;
    int startIdx, lastIdx;
    bool isClosed;
    bool fine;
    bool canBeClosed(ofxPolygonDetector& pd, int idToAdd) const;
    bool contains(int idP) const;
    bool addLineId(ofxPolygonDetector& pd, int id);
    string idxToString() const;
    string toString() const;
    bool equals(const PolyCycle& p) const;
    bool convex(ofxPolygonDetector& pd) const;
    bool pointConsumed(ofxPolygonDetector& pd, int pid) const;
    bool accepted(ofxPolygonDetector& pd);
};

using PolyCycles = vector<PolyCycle>;

struct ofxPolyLine
{
    PointType a, b;

    ofxPolyLine() {}
    ofxPolyLine(const PointType& aP, const PointType& bP) : a(aP), b(bP) {}

    vector<int> intersections;
    set<int> intersectedLines;
    PointType center;
    int id = 0;
    int test0 = 0, test1 = 0;
    int32_t origLine = -1;
    int32_t attr0 = -1; // used for testing
    bool ignore = false;
    int took = 0;
    int processed = 0;
    int aIdx = 0, bIdx = 0;
    int lastDissolveStep = 0;
    void calcCenter();
    bool hasCommonIdxPoints(const ofxPolyLine& line) const;
    void sortIntersectionsList(ofxPolygonDetector& pd);
    bool IntersectionPoint(const ofxPolyLine& line, PointType& pos) const;
    bool lineLineIntersectionPoint(const ofxPolyLine& line, PointType& pos) const;
    void calculateFirstAndLastPoint();
    static bool bCompareLineOrder(const ofxPolyLine& l1, ofxPolyLine& l2);
    static int iCompareLineOrder(const ofxPolyLine& l1, ofxPolyLine& l2);
    string neighToString(ofxPolygonDetector& pd, int* retNNeigh = nullptr) const;
    string toString(ofxPolygonDetector& pd) const;
    int numNeigh(ofxPolygonDetector& pd) const;
    int numIntersections(ofxPolygonDetector& pd) const;
    int canBeRemoved(ofxPolygonDetector& pd, RmLinesType type) const;
    ofxPolyLine& mul(float m);
    ofxPolyLine& add(const PointType& p);
    void setIgnore(ofxPolygonDetector& pd, const char* msg);
    bool contains(const PointType& point) const;
    bool contains(const ofxPolyLine& line) const;
    bool collinear(const ofxPolyLine& l) const;
    int minPid() const;
    int maxPid() const;
    int32_t commonPid(const ofxPolyLine& l) const;
    int otherPid(int pid) const;
    bool compareNeigh(ofxPolygonDetector& pd, int nid1, int nid2) const;
    bool sortNeigh(ofxPolygonDetector& pd) const;
    int& incTook(ofxPolygonDetector& pd);
    bool betweenNeighbors(ofxPolygonDetector& pd, const ofxPolyLine& l1, const ofxPolyLine& l2) const;
};

struct ofxPolyPol
{
    vector<PointType> p;
    PointType c;
    size_t firstIdx = 0;
    ofColor color;
    PolyCycle cycle;
    int id;
    int dissolveStep = 0;
    double _area;
    int getCount() const;
    void calculateFirstAndLastPoint();
    bool minus(const ofxPolyPol& other);
    void addLine(const ofxPolyLine& l);
    int roundArea();
    bool addPointChecked(const PointType& v);
    void setColor(ofColor color);
    void draw();
    PointType center();
    double triangleArea(ofxPolygonDetector& pd);
};

struct ofxPolygonDetector
{
    using LineVector = vector<ofxPolyLine>;
    using LineList = list<ofxPolyLine>;
    using PolyVector = vector<ofxPolyPol>;

    PolyCycles _cycles;

    LineVector
        origLines, // from user
        lines; // active lines
    PolyVector polys;

    map<int, vector<int>> _neighbors; // key: lid, val: vec(neighborLine)
    map<int, vector<int>> collinearLineMap; // key: lid, val: vec(collinearLine)

    vector<PointType> intersectionPoints;
    map<int, vector<int>> pointToLines; // key: pid, val: list of nids
    unordered_map<int, int32_t> lineIdToIdx; // key: lid, val: index in lines

    int dissolveCount = 0;

    void reset();
    void addLine(const ofxPolyLine& line);
    bool detectPolygons(LineVector lineVector);
    PolyVector& getPolys() { return polys; }

    bool addPointToLine(int pid, int lid);

    ofxPolyLine* findLine(int id, bool useIgnore = true);
    ofxPolyLine* findLine(int pidA, int pidB, bool useIgnore = true);
    ofxPolyLine* findOrigLine(int id);

    int getPolyCount() const { return (int)polys.size(); };
    void sortLines();

    // The order of operations
    void removeZeroLengthLines();
    void removeOverlappings(); 
    int detectAllIntersections();
    bool createLines();
    bool findPolys();
    void simplifyPolys(double smaller_polygon_length);

    bool buildCycle(int id, PolyCycle cycle); // as value!

    ofxPolyLine newLine(int i, int j, ofxPolyLine& origLine);

    // collinearity
    void setCollinear(int l1, int l2);
    bool collinearIdx(int l1, int l2);
    bool collinearIdx(const ofxPolyLine& l1, const ofxPolyLine& l2);

    bool rmLines(RmLinesType type);
    bool rmEarPoints();

    // dissolve
    bool dissolveCollinearLine(ofxPolyLine& l);
    bool dissolveCollinear(ofxPolyLine& l1, ofxPolyLine& l2);
    bool dissolve();

    void dumpLines(const char* msg, bool useIgnore = false);
};
