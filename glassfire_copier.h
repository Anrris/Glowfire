#ifndef GLASSFIRE_COPIER_H
#define GLASSFIRE_COPIER_H

#include <boost/geometry.hpp>                   /* Boost Library: For RTree*/
#include <boost/geometry/geometries/point.hpp>  /* Boost Library */
#include <boost/geometry/geometries/box.hpp>    /* Boost Library */
#include <boost/geometry/index/rtree.hpp>       /* Boost Library */

namespace bg  = boost::geometry;
namespace bgm = boost::geometry::model;
namespace bgi = boost::geometry::index;

template<typename AxisType, size_t Dimension>
class point_copier{
    typedef std::vector<AxisType> Feature;
    typedef bgm::point<AxisType, Dimension, bg::cs::cartesian>  RtreePoint;

    template<size_t _x_, size_t to>
    struct t_static_for
    {
        typedef std::vector<AxisType> Feature;

        AxisType m_val_min = std::numeric_limits<AxisType>::max();
        AxisType m_val_max = std::numeric_limits<AxisType>::min();

        void checkDimension(size_t dim){
            if(dim != Dimension){
                throw std::length_error("Input data does not match required dimension.");
            }
        }

        void copy_to(Feature &coord, RtreePoint &pt) {
            checkDimension(coord.size());

            m_val_min = std::min(m_val_min, coord.at(_x_));
            m_val_max = std::max(m_val_max, coord.at(_x_));

            pt.template set<_x_>(coord.at(_x_));
            t_static_for<_x_+1,to>().copy_to(coord, pt);
        }

        void create_box_corner(AxisType dist, Feature &coord, RtreePoint &pt1, RtreePoint &pt2) {
            checkDimension(coord.size());

            pt1.template set<_x_>(coord.at(_x_)-dist);
            pt2.template set<_x_>(coord.at(_x_)+dist);
            t_static_for<_x_+1,to>().create_box_corner(dist, coord, pt1, pt2);
        }
    };

    template<size_t to>
    struct t_static_for<to, to>
    {
        void copy_to(Feature &arr, RtreePoint &pt) {}
        void create_box_corner(AxisType dist, Feature &coord, RtreePoint &pt1, RtreePoint &pt2) {}
    };

    typedef t_static_for<0, Dimension>    static_for;

public:
    void copy_to(Feature &arr, RtreePoint &p){
        static_for  m_static_for;
        m_static_for.copy_to(arr, p);
    }

    RtreePoint to_point(Feature arr){
        RtreePoint p;
        static_for  m_static_for;
        m_static_for.copy_to(arr, p);
        return p;
    }

    void create_box_corner(AxisType dist, Feature &coord, RtreePoint &pt1, RtreePoint &pt2) {
        static_for  m_static_for;
        m_static_for.create_box_corner(dist, coord, pt1, pt2);
    }
};


#endif