#ifndef DEMSI_POLYGON_H_
#define DEMSI_POLYGON_H_

#include <vector>
#include <string>

/*! \file   demsi_polygon.h
    \brief  Header file for the DEMSI::Polygon class
*/

namespace DEMSI {

//------------------------------------------------------------------------------
/*! \struct Vector2D
    \brief This struct defines a two dimensional vector with x and y components.
*/
//------------------------------------------------------------------------------
struct Vector2D {
  /*! x component of two dimensional vector. */
  double x;
  /*! y component of two dimensional vector. */
  double y;
};

//------------------------------------------------------------------------------
/*! \class Polygon
    \brief This class defines a two dimensional polygon.

    This class defines a two dimensional convex polygon as a vector of vertices
    consisting of Vector2D structs. The class defines method to find the
    intersection of two polygons (also a polygon), and the area of the polygon.
    Polygon vertex order is counter-clockwise and repeated first/last vertex is
    not used.

    Intersection method:
    O'Rourke, J., C. Chien, T. Olson, and D. Naddor (1982), A New Linear
    Algorithm for Intersecting Convex Polygons, Computer Graphics and Image
    Processing, 19, 384-391

    Previous intersection method:
    http://rosettacode.org/wiki/Sutherland-Hodgman_polygon_clipping#C
    https://en.wikipedia.org/wiki/Sutherland-Hodgman_algorithm
    Ivan Sutherland, Gary W. Hodgman: Reentrant Polygon Clipping.
    Communications of the ACM, vol. 17, pp. 32–42, 1974

    Inclusion of point in polygon:
    https://en.wikipedia.org/wiki/Point_in_polygon#cite_ref-6
    Dan Sunday (2001)
    http://geomalgorithms.com/a03-_inclusion.html
*/
//------------------------------------------------------------------------------
class Polygon {

public:

  /*! \brief Polygon constructor
   */
  Polygon(void) {};

  /*! \brief Polygon constructor
      \param xVertices x position of polygon vertices
      \param yVertices y position of polygon vertices
  */
  Polygon(const std::vector<double> xVertices, const std::vector<double> yVertices);

  /*! \brief Polygon constructor
      \param verticesIn Array of polygon vertex coordinates
  */
  Polygon(const std::vector<Vector2D> verticesIn);

  /*! \brief Default destructor.*/
  ~Polygon() = default;

  /*! \brief Add a vertex to the polygon at the end of the vertex list
      \param x x position of the new vertex
      \param y y position of the new vertex
  */
  void add_vertex(const double x, const double y);

  /*! \brief Find the intersection of this polygon and another
      \param other The other polygon to find the intersection with.
      \return The intersection polygon of this and the other polygon.
  */
  DEMSI::Polygon intersection(const DEMSI::Polygon other) const;

  /*! \brief Return the number of vertices the polygon has.
      \return The number of vertices the polygon has.
  */
  int nvertices(void) const;

  /*! \brief Return the area of the polygon.
      \return The area of the polygon.
  */
  double area(void) const;

  /*! \brief Get the x position of a given vertex
      \param iVertex vertex index to get position of
      \return The x position of a given vertex
  */
  double get_x(const int iVertex) const;

  /*! \brief Get the y position of a given vertex
      \param iVertex vertex index to get position of
      \return The y position of a given vertex
  */
  double get_y(const int iVertex) const;

  /*! \brief Get the maximum size (vertex to vertex) of the polygon
      \return The maximum size (vertex to vertex) of the polygon
  */
  double max_size(void) const;

  /*! \brief Write the polygon vertex locations to a text file.
      \param filename The name of the text file to write to.
  */
  void write(const std::string filename) const;

  friend std::ostream & operator<<(std::ostream & os, const Polygon & polygon);

private:

  /*! Array of vertex positions for the polygon */
  std::vector<Vector2D> vertices;

};

} // namespace DEMSI

#endif // DEMSI_POLYGON_H_
