#include "demsi_polygon.h"

#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>

namespace DEMSI {

// Polygon constructor with separate x and y vertex position vectors
Polygon::Polygon(const std::vector<double> xVertices, const std::vector<double> yVertices) {

  if (xVertices.size() == yVertices.size()) {

    for (int iVertex = 0 ; iVertex < xVertices.size() ; iVertex++) {

      vertices.push_back({xVertices[iVertex],yVertices[iVertex]});

    } // iVertex

  }

}

// Polygon constructor with array of Vector2D structs
Polygon::Polygon(const std::vector<Vector2D> verticesIn) {

  vertices = verticesIn;

}

// add a vertex to the polygon vertices
void Polygon::add_vertex(const double x, const double y) {

  vertices.push_back({x,y});

}

// cross product of two 2D vectors
double cross_product(Vector2D a, Vector2D b) {

  return a.x * b.y - a.y * b.x;

}

// test if a point is in a half plane
bool in_half_plane(const Vector2D x, const Vector2D pdot, const Vector2D pn) {

  Vector2D vec;
  vec.x = x.x - pn.x;
  vec.y = x.y - pn.y;

  return (cross_product(pdot,vec) >= 0.0);

}

// solve a 2 by 2 matrix system
void solve_2by2_matrix_eqn(const double a, const double b, const double c, const double d, const double e, const double f, double& x, double& y) {

  double det = a*d - b*c;

  x = (d*e - b*f) / det;
  y = (a*f - c*e) / det;

}

// determine if two line segments overlap and if they do calculate intersection point
bool edge_intersection(const Vector2D p0, const Vector2D p1, const Vector2D q0, const Vector2D q1, Vector2D& i) {

  double alpha, beta;

  solve_2by2_matrix_eqn(p1.x - p0.x,
			q0.x - q1.x,
			p1.y - p0.y,
			q0.y - q1.y,
			q0.x - p0.x,
			q0.y - p0.y,
			alpha, beta);

  i.x = p0.x + alpha * (p1.x - p0.x);
  i.y = p0.y + alpha * (p1.y - p0.y);

  return (alpha >= 0.0 and alpha <= 1.0 and
	  beta  >= 0.0 and beta  <= 1.0);

}

// determine if a point is within a polygon
bool point_in_polygon_old(const Vector2D point, const std::vector<Vector2D> polygon) {

  // this function determines the sum of angles subtended by the edges in a polygon
  // to a point. If the point is inside the polygon this sum is 2 pi, whereas it is
  // zero if the point is outside the polygon
  double angleSum = 0.0;

  for (int i0 = 0 ; i0 < polygon.size() ; i0++) {
    int i1 = (i0 + 1) % polygon.size();

    Vector2D a;
    a.x = polygon[i0].x - point.x;
    a.y = polygon[i0].y - point.y;
    double amag = std::sqrt(a.x*a.x + a.y*a.y);

    Vector2D b;
    b.x = polygon[i1].x - point.x;
    b.y = polygon[i1].y - point.y;
    double bmag = std::sqrt(b.x*b.x + b.y*b.y);

    angleSum += std::asin(cross_product(a,b) / (amag * bmag));

  }

  return (angleSum > 3.14);

}

double is_left(Vector2D P0, Vector2D P1, Vector2D P2) {
  return ( (P1.x - P0.x) * (P2.y - P0.y) - (P2.x - P0.x) * (P1.y - P0.y) );
}

bool point_in_polygon(const Vector2D point, const std::vector<Vector2D> polygon) {

  // http://geomalgorithms.com/a03-_inclusion.html

  int windingNumber = 0;

  for (int i0 = 0 ; i0 < polygon.size() ; i0++) {
    int i1 = (i0 + 1) % polygon.size();

    if (polygon[i0].y <= point.y) {
      if (polygon[i1].y > point.y and is_left(polygon[i0], polygon[i1], point) > 0.0) {
	++windingNumber;
      }
    } else {
      if (polygon[i1].y <= point.y and is_left(polygon[i0], polygon[i1], point) < 0.0) {
	--windingNumber;
      }
    }

  }

  return (windingNumber != 0);

}

// calculate the intersection of two convex polygons
std::vector<Vector2D> polygon_intersection(const std::vector<Vector2D> P, const std::vector<Vector2D> Q) {

  // O'Rourke, J., C. Chien, T. Olson, and D. Naddor (1982), A New Linear Algorithm for Intersecting Convex Polygons, Computer Graphics and Image Processing, 19, 384-391

  std::vector<Vector2D> PnQ;

  // number of vertices on polygon
  int L = P.size();
  int M = Q.size();

  int ip = 0, ipn = L-1, ipp = 0; // index of P polygon vertex
  int iq = 0, iqn = M-1, iqp = 0; // index of Q polygon vertex

  // edge vectors
  Vector2D pdot, qdot;

  pdot.x = P[ipp].x - P[ipn].x;
  pdot.y = P[ipp].y - P[ipn].y;

  qdot.x = Q[iqp].x - Q[iqn].x;
  qdot.y = Q[iqp].y - Q[iqn].y;

  // logical whether p or q is inside
  // set to false initially so that we dont start recording vertices until we have found first intersection
  bool pInside = false;
  bool qInside = false;

  // store intersections so we know if are back to first found intersection
  std::vector<Vector2D> intersections;

  // loop through maximum number of advancements we need to find polygon intersection
  for (int k = 0 ; k < 2 * (L + M); k++) {

    // check if found intersection
    Vector2D intersection;
    if (edge_intersection(P[ipn], P[ipp], Q[iqn], Q[iqp], intersection)) {

      // check if got back to initial intersection
      if (intersections.size() > 0 and intersection.x == intersections[0].x and intersection.y == intersections[0].y) {

	// back to original intersection to have found intersection polygon
	break;

      } else {

	// add intersection to output polygon
	PnQ.push_back(intersection);

	// check if p in q halfPlane
	pInside = in_half_plane(P[ipp], qdot, Q[iqn]);
	qInside = not pInside;

      }

      // add new intersection to store of intersections
      intersections.push_back(intersection);

    }

    // advance p or q
    if (cross_product(qdot,pdot) >= 0.0) {

      if (in_half_plane(P[ipp], qdot, Q[iqn])) {

	// advance q
	if (qInside) PnQ.push_back(Q[iqp]);
	iq++;
	iqn = (iq-1) % M;
	iqp = iq % M;

	qdot.x = Q[iqp].x - Q[iqn].x;
	qdot.y = Q[iqp].y - Q[iqn].y;

      } else {

	// advance p
	if (pInside) PnQ.push_back(P[ipp]);
	ip++;
	ipn = (ip-1) % L;
	ipp = ip % L;

	pdot.x = P[ipp].x - P[ipn].x;
	pdot.y = P[ipp].y - P[ipn].y;

      }

    } else {

      if (in_half_plane(Q[iqp], pdot, P[ipn])) {

	// advance p
	if (pInside) PnQ.push_back(P[ipp]);
	ip++;
	ipn = (ip-1) % L;
	ipp = ip % L;

	pdot.x = P[ipp].x - P[ipn].x;
	pdot.y = P[ipp].y - P[ipn].y;

      } else {

	// advance q
	if (qInside) PnQ.push_back(Q[iqp]);
	iq++;
	iqn = (iq-1) % M;
	iqp = iq % M;

	qdot.x = Q[iqp].x - Q[iqn].x;
	qdot.y = Q[iqp].y - Q[iqn].y;

      }

    }

  } // loop while advancing

  // check if didnt find intersection
  if (PnQ.size() == 0) {

    // see if P in Q
    if (point_in_polygon(P[0], Q)) {

      for (int i = 0 ; i < P.size() ; i++) {
	PnQ.push_back(P[i]);
      } // i

    } else {

      // P not in Q; see if Q in P
      if (point_in_polygon(Q[0], P)) {

	for (int i = 0 ; i < Q.size() ; i++) {
	  PnQ.push_back(Q[i]);
	} // i

      } // Q in P

    } // P in Q

  } // P and Q intersect

  return PnQ;

}

// Find the intersection polygon between this Polygon and another
Polygon Polygon::intersection(const Polygon other) const {

  std::vector<Vector2D> poly = polygon_intersection(this->vertices, other.vertices);
  return Polygon(poly);

}

// Return the number of vertices the polygon has.
int Polygon::nvertices(void) const {
  return this->vertices.size();
}

// Return the area of the polygon.
double Polygon::area(void) const {

  int nVertices = this->nvertices();
  double areaOut = 0; // Accumulates area in the loop
  int j = nVertices-1; // The last vertex is the 'previous' one to the first

  for (int i = 0 ; i < nVertices ; i++) {
    areaOut += (this->vertices[i].x + this->vertices[j].x) * (this->vertices[i].y - this->vertices[j].y);
    j = i; // j is previous vertex to i
  }
  return areaOut / 2.0;

}

// get x value of given vertex
double Polygon::get_x(const int iVertex) const {

  return vertices[iVertex].x;

}

// get y value of given vertex
double Polygon::get_y(const int iVertex) const {

  return vertices[iVertex].y;

}

// get the maximum size (vertex to vertex) of the polygon
double Polygon::max_size(void) const {

  double maxSize = 0.0;

  for (int iVertex = 1 ; iVertex < this->vertices.size() ; iVertex++) {
    maxSize = std::max(maxSize,
		       std::sqrt(std::pow(this->vertices[iVertex].x-this->vertices[0].x,2) +
				 std::pow(this->vertices[iVertex].y-this->vertices[0].y,2)));
  } // iVertex

  return maxSize;

}


// write the polygon vertex positions to a text file.
void Polygon::write(const std::string filename) const {

  std::ofstream myfile;
  myfile.open(filename.c_str());
  myfile << std::setprecision(12);
  for (int i = 0 ; i < this->vertices.size() ; i++) {
    myfile << this->vertices[i].x << " " << this->vertices[i].y << std::endl;
  } // i
  if (this->vertices.size() > 0) myfile << this->vertices[0].x << " " << this->vertices[0].y << std::endl;
  myfile.close();

}

// over load << operator
std::ostream & operator<<(std::ostream & os, const Polygon & polygon)
{
  os << "Polygon: {nVertices: " << polygon.vertices.size() << "}";
  return os;
}

} // namespace DEMSI
