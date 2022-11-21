#ifndef EULER_CLASSES_H
#define EULER_CLASSES_H

#include <vector>
#include <array>
#include <string>
#include <math.h>
#include <iostream>
#include <map>


typedef unsigned int uint;
typedef double real;


struct vec2 {
    real x;
    real y;
    vec2() {x=0.; y=0;}
    vec2(real a) {x=a; y=a;}
    vec2(real a, real b) {x=a; y=b;}
    vec2 operator + (const vec2 &v) const {
        return vec2(x + v.x, y + v.y);
    }
    vec2& operator += (const vec2 &v) {
        x += v.x;
        y += v.y;
        return *this;
    }
    vec2 operator - (const vec2 &v) const {
        return vec2(x - v.x, y - v.y);
    }
    vec2& operator -= (const vec2 &v) {
        x -= v.x;
        y -= v.y;
        return *this;
    }
    vec2 operator * (const real &v) const {
        return vec2(x*v, y*v);
    }
    vec2& operator *= (const real &v) {
        x *= v;
        y *= v;
        return *this;
    }
    vec2 operator / (const real &v) const {
        return vec2(x/v, y/v);
    }
    vec2& operator /= (const real &v) {
        x /= v;
        y /= v;
        return *this;
    }
};

real norm(const vec2 &p) {
    return sqrt(p.x*p.x + p.y*p.y);
}

real dot(const vec2& a, const vec2& b) {
    return a.x*b.x + a.y*b.y;
}

real cross(const vec2& a, const vec2& b) {
    return a.x*b.y - a.y*b.x;
}

std::ostream& operator << (std::ostream& o, const vec2& a)
{
    o << a.x << " " << a.y;
    return o;
}


struct inputVar {
    real rho, u, v, p;
    inputVar() {rho=1.; u=0.; v=0.; p=0.;}
    inputVar(real a, real b, real c, real d) {
        rho=a;
        u=b;
        v=c;
        p=d;
    }
    int size() const {
        return 4;
    }
    real& operator[](uint idx) {
        if (idx == 0) { return rho; }
        else if (idx == 1) { return u; }
        else if (idx == 2) { return v; }
        else if (idx == 3) { return p; }
    }
    const real& operator[](uint idx) const {
        if (idx == 0) { return rho; }
        else if (idx == 1) { return u; }
        else if (idx == 2) { return v; }
        else if (idx == 3) { return p; }
    }
};

struct var {
    real rho, rhou, rhov, rhoe;
    var() {rho=1.; rhou=0.; rhov=0.; rhoe=0.;}
    var(real a) {
        rho=a;
        rhou=a;
        rhov=a;
        rhoe=a;
    }
    var(real a, real b, real c, real d) {
        rho=a;
        rhou=b;
        rhov=c;
        rhoe=d;
    }
    var& operator = (const real &v) {
        rho = v;
        rhou = v;
        rhov = v;
        rhoe = v;
        return *this;
    }
    var operator + (const var &v) const {
        return var(rho+v.rho, rhou+v.rhou, rhov+v.rhov, rhoe+v.rhoe);
    }
    var operator + (const real &v) const {
        return var(rho+v, rhou+v, rhov+v, rhoe+v);
    }
    var& operator += (const var &v) {
        rho += v.rho;
        rhou += v.rhou;
        rhov += v.rhov;
        rhoe += v.rhoe;
        return *this;
    }
    var operator - (const var &v) const {
        return var(rho-v.rho, rhou-v.rhou, rhov-v.rhov, rhoe-v.rhoe);
    }
    var& operator -= (const var &v) {
        rho -= v.rho;
        rhou -= v.rhou;
        rhov -= v.rhov;
        rhoe -= v.rhoe;
        return *this;
    }
    var operator * (const real &s) const {
        return var(rho*s, rhou*s, rhov*s, rhoe*s);
    }
    var operator * (const var &s) const {
        return var(rho*s.rho, rhou*s.rhou, rhov*s.rhov, rhoe*s.rhoe);
    }
    var& operator *= (const real &v) {
        rho *= v;
        rhou *= v;
        rhov *= v;
        rhoe *= v;
        return *this;
    }
    var operator / (const real &s) const {
        return var(rho/s, rhou/s, rhov/s, rhoe/s);
    }
    var operator / (const var &s) const {
        return var(rho/s.rho, rhou/s.rhou, rhov/s.rhov, rhoe/s.rhoe);
    }
    var& operator /= (const real &v) {
        rho /= v;
        rhou /= v;
        rhov /= v;
        rhoe /= v;
        return *this;
    }
    var& operator /= (const var &v) {
        rho /= v.rho;
        rhou /= v.rhou;
        rhov /= v.rhov;
        rhoe /= v.rhoe;
        return *this;
    }

    int size() const {
        return 4;
    }
    real& operator[](uint idx) {
        if (idx == 0) { return rho; }
        else if (idx == 1) { return rhou; }
        else if (idx == 2) { return rhov; }
        else if (idx == 3) { return rhoe; }
        else { return rho; }
    }
    const real& operator[](uint idx) const {
        if (idx == 0) { return rho; }
        else if (idx == 1) { return rhou; }
        else if (idx == 2) { return rhov; }
        else if (idx == 3) { return rhoe; }
        else { return rho; }
    }
};

bool is_var_nan(var& v) {
    if ((std::isnan(v.rho)|std::isnan(v.rhou))|(std::isnan(v.rhov)|std::isnan(v.rhoe))) {
        return true;
    } else {
        return false;
    }
}

var sqrt(const var& v) {
    return var(
        sqrt(v.rho),
        sqrt(v.rhou),
        sqrt(v.rhov),
        sqrt(v.rhoe)
    );
}


struct varVec2 {
    var x;
    var y;

    varVec2(var a, var b) {x = a; y = b;}
    varVec2(var a) {x = a; y = a;}
    varVec2(real a) {x = a; y = a;}
    varVec2() {x = 0.; y = 0.;}

    varVec2 operator + (const real &v) const {
        return varVec2(x+v, y+v);
    }
    varVec2& operator += (const varVec2 &v) {
        x += v.x;
        y += v.y;
        return *this;
    }
    varVec2 operator + (const varVec2 &v) const {
        return varVec2(x+v.x, y+v.y);
    }
    varVec2 operator - (const varVec2 &v) const {
        return varVec2(x-v.x, y-v.y);
    }
    varVec2& operator -= (const varVec2 &v) {
        x -= v.x;
        y -= v.y;
        return *this;
    }
    varVec2 operator / (const varVec2 &v) const {
        return varVec2(x/v.x, y/v.y);
    }
    varVec2& operator /= (const real &v) {
        x /= v;
        y /= v;
        return *this;
    }
    varVec2 operator * (const varVec2 &v) const {
        return varVec2(x*v.x, y*v.y);
    }
    varVec2 operator * (const var &v) const {
        return varVec2(x*v, y*v);
    }
    varVec2 operator * (const real &v) const {
        return varVec2(x*v, y*v);
    }
    varVec2 operator * (const vec2 &b) const {
        return varVec2(
            x*b.x,
            y*b.y
        );
    }
};

var dot(const varVec2& a, const vec2& b) {
    return a.x*b.x + a.y*b.y;
}

varVec2 prodVarVec2(const var& v, const vec2 &s) {
    return varVec2(v*s.x, v*s.y);
}


namespace std {
    var abs(const var& v) {
        return var(
            std::abs(v.rho),
            std::abs(v.rhou),
            std::abs(v.rhov),
            std::abs(v.rhoe)
        );
    }
    real max(const var& r) {
        return max(r.rho, max(r.rhou, max(r.rhov, r.rhoe)));
    }
    var max(const var& a, const var& b) {
        return var(
            max(a.rho, b.rho),
            max(a.rhou, b.rhou),
            max(a.rhov, b.rhov),
            max(a.rhoe, b.rhoe)
        );
    }
    var min(const var& a, const var& b) {
        return var(
            min(a.rho, b.rho),
            min(a.rhou, b.rhou),
            min(a.rhov, b.rhov),
            min(a.rhoe, b.rhoe)
        );
    }
    inputVar max(const inputVar& a, const inputVar& b) {
        return inputVar(
            max(a.rho, b.rho),
            max(a.u, b.u),
            max(a.v, b.v),
            max(a.p, b.p)
        );
    }
    inputVar min(const inputVar& a, const inputVar& b) {
        return inputVar(
            min(a.rho, b.rho),
            min(a.u, b.u),
            min(a.v, b.v),
            min(a.p, b.p)
        );
    }
    varVec2 max(const varVec2& a, const varVec2& b) {
        return varVec2(
            max(a.x, b.x),
            max(a.y, b.y)
        );
    }
    varVec2 min(const varVec2& a, const varVec2& b) {
        return varVec2(
            min(a.x, b.x),
            min(a.y, b.y)
        );
    }
    real sign(const real& x) {
        if (x > 0) { return 1.;}
        else if (x < 0) { return -1.;}
        else { return 0.; }
    }
    var sign(const var& v) {
        return var(
            sign(v.rho),
            sign(v.rhou),
            sign(v.rhov),
            sign(v.rhoe)
        );
    }
}

std::ostream& operator << (std::ostream& o, const var& a)
{
    o << a.rho << " " << a.rhou << " " << a.rhov << " " << a.rhoe;
    return o;
}


struct constants {
    real gamma;
    constants() {}
    constants(const real gin) {gamma = gin;}
};


var convertVariable(const inputVar& q, const constants& consts) {
    return var(
        q.rho,
        q.u*q.rho,
        q.v*q.rho,
        q.p/(consts.gamma - 1.) + 0.5*q.rho*(q.u*q.u + q.v*q.v)
    );
}

inputVar convertVariable(const var& q, const constants& consts) {
    return inputVar(
        q.rho,
        q.rhou/q.rho,
        q.rhov/q.rho,
        (consts.gamma - 1.)*(q.rhoe - 0.5/q.rho*(q.rhou*q.rhou + q.rhov*q.rhov))
    );
}


// Mesh struct
// Supports only triangles and quads
// Simple vector form
struct mesh {
    std::vector<vec2> nodes;
    std::vector<real> nodesGradWeights;

    std::vector<vec2> edgesNormals;
    std::vector<vec2> edgesCenters;
    std::vector<uint> edgesFaces0;
    std::vector<uint> edgesFaces1;
    std::vector<uint> edgesNodes0;
    std::vector<uint> edgesNodes1;
    std::vector<real> edgesLengths;

    std::map<std::tuple<uint, uint>, uint> edgesRef;

    std::vector<vec2> boundsNormals;
    std::vector<vec2> boundsCenters;
    std::vector<uint> boundsFaces;
    std::vector<uint> boundsNodes0;
    std::vector<uint> boundsNodes1;
    std::vector<real> boundsLengths;
    std::vector<std::string> boundsTags;
    std::vector<var (*)(const var&, const var&, const varVec2&, const vec2&, const vec2&, const constants&)> boundsFuncs;

    std::map<std::tuple<uint, uint>, uint> boundsRef;

    std::vector<uint> facesNodes0;
    std::vector<uint> facesNodes1;
    std::vector<uint> facesNodes2;
    std::vector<uint> facesNodes3;
    std::vector<bool> faceIsTriangle;
    std::vector<real> facesAreas;
    std::vector<vec2> facesCenters;

    std::vector<std::string> uniqueBoundsTags;

    real import_duration;

    mesh() {};
    mesh(const std::string filename, const bool verb);

    void compute_mesh();
    void add_last_quad_edges();
    void add_last_tri_edges();
    void read_su2_nodes(std::vector<std::string>& content);
    void read_su2_elements(std::vector<std::string>& content);
    void read_su2_boundaries(std::vector<std::string>& content);
    void read_raw_su2(const std::string filename, const bool verb);

    void convert_node_face_info();
    bool find_if_edge_in_mesh(const uint n0, const uint n1);

    uint find_edge_with_nodes(const uint n0, const uint n1);
    uint find_bound_with_nodes(const uint n0, const uint n1);

    int size() { return facesAreas.size(); }
    std::tuple<real, real> get_xy(const int i) {
        return std::make_tuple(facesCenters[i].x, facesCenters[i].y);
    }
};


#endif
