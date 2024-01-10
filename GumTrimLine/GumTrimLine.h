#ifndef GUMTRIMLINE_H
#define GUMTRIMLINE_H

#include <array>
#include <functional>
#include <string>
#include <vector>
#include <iterator>
#include "../Polyhedron.h"
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Weighted_point_3.h>
#include "tinycolormap.hpp"
#include "../EasyOBJ.h"
#include "../Ortho.h"

namespace internal
{
    template <typename Curve>
    class CurveIterator;

    template <typename Kernel>
    class Curve
    {
    public:
        using Point_3 = Kernel::Point_3;
        using Vector_3 = Kernel::Vector_3;
        using K = Kernel;
        using iterator = CurveIterator<Curve<Kernel>>;
        using const_iterator = CurveIterator<const Curve<Kernel>>;

        void AddPoint(Point_3 p, int label)
        {
            _points.push_back(p);
            _labels.push_back(label);
        }
        Point_3& operator[](size_t i)
        {
            return _points[i];
        }
        const Point_3& operator[](size_t i) const
        {
            return _points[i];
        }
        int& Label(size_t i)
        {
            return _labels[i];
        }
        const int& Label(size_t i) const
        {
            return _labels[i];
        }
        void UpdateData()
        {
            std::unordered_map<int, std::vector<Point_3>> labelwise_points;
            for(size_t i = 0; i < _points.size(); i++)
            {
                labelwise_points[_labels[i]].push_back(_points[i]);
            }
            for(const auto& [label, points] : labelwise_points)
            {
                typename Kernel::Point_3 c;
                typename Kernel::Plane_3 plane;
                CGAL::linear_least_squares_fitting_3(points.begin(), points.end(), plane, c, CGAL::Dimension_tag<0>());
                _centroids.insert({label, c});
                typename Kernel::Vector_3 dir = plane.perpendicular_line(c).direction().to_vector();
                _upwards.insert({label, dir});
            }

            // unify all directions. Is this robust?
            for(auto& [label, dir] : _upwards)
            {
                auto p0 = _centroids[label];
                auto p1 = p0 + dir;
                double signed_v = 0.0;
                for(size_t i = 0; i < _points.size(); i++)
                {
                    auto p2 = _points[(i + 1) % _points.size()];
                    auto p3 = _points[i];
                    signed_v += CGAL::determinant(p1 - p0, p2 - p0, p3 - p0);
                }
                if(signed_v < 0)
                {
                    dir = -dir;
                }
            }

            // collect all appearing labels
            std::array<bool, 50> label_map;
            std::fill(label_map.begin(), label_map.end(), false);
            for(int label : _labels)
            {
                if(label >= 0 && label <= 49)
                {
                    label_map[label] = true;
                }
            }
            {
                std::vector<int> labels0;
                std::vector<int> labels1;
                for(int label = 0; label < 50; label++)
                {
                    if(label_map[label])
                    {
                        if(label >= 11 && label <= 18 || label >= 31 && label <= 38)
                        {
                            labels0.push_back(label);
                        }
                        else if(label >= 21 && label <= 28 || label >= 41 && label <= 48)
                        {
                            labels1.push_back(label);
                        }
                    }
                }
                std::sort(labels0.begin(), labels0.end());
                std::reverse(labels0.begin(), labels0.end());
                std::sort(labels1.begin(), labels1.end());
                _ordered_labels.insert(_ordered_labels.end(), labels0.begin(), labels0.end());
                _ordered_labels.insert(_ordered_labels.end(), labels1.begin(), labels1.end());
            }
        }
        void LoadCrownFrame( const CrownFrames<Kernel>& frames )
        {
            for(auto& [label, up] : _upwards)
            {
                up = frames.GetFrame(label).up;
            }
            for(auto& [label, c] : _centroids)
            {
                c = frames.GetFrame(label).pos;
            }
        }
        int MaxLabel() const
        {
            return _ordered_labels.back();
        }
        int MinLabel() const
        {
            return _ordered_labels.front();
        }
        size_t size() const { return _points.size(); }
        std::vector<Point_3> GetPoints() const { return _points; }
        Point_3 GetCentroidOfLabel(int label) const { return _centroids.at(label); }
        Vector_3 GetUpwardOfLabel(int label) const { return _upwards.at(label); }
        Curve<Kernel> GetSubCurve(size_t first, size_t last) const
        {
            if(first >= _points.size() || last >= _points.size())
            {
                throw std::range_error("out of curve range");
            }
            Curve<Kernel> subcurve;
            size_t idx = first;
            subcurve.AddPoint(_points[idx], _labels[idx]);
            do
            {
                idx++;
                if(idx == _points.size())
                {
                    idx = 0;
                }
                subcurve.AddPoint(_points[idx], _labels[idx]);
            } while (idx != last);
            return subcurve;
        }
        Curve<Kernel> GetSubCurve(const const_iterator& first, const const_iterator& last) const
        {
            return GetSubCurve(first.Idx(), last.Idx());
        }
        Curve<Kernel> GetSubCurve(const iterator& first, const iterator& last) const
        {
            return GetSubCurve(first.Idx(), last.Idx());
        }
        void InsertAt(const Curve<Kernel>& curve, size_t pos)
        {
            if(pos > _points.size())
            {
                throw std::range_error("out of curve range");
            }
            _points.insert(_points.begin() + pos, curve._points.begin(), curve._points.end());
            _labels.insert(_labels.begin() + pos, curve._labels.begin(), curve._labels.end());
        }
        std::vector<Point_3>::iterator begin() { return _points.begin(); }
        std::vector<Point_3>::iterator end() { return _points.end(); }
        std::vector<Point_3>::const_iterator begin() const { return _points.begin(); }
        std::vector<Point_3>::const_iterator end() const { return _points.end(); }
        void WriteOBJ(const std::string& path) const
        {
            static const std::array<std::array<float, 3>, 10> COLORS = {
                std::array<float, 3>{142.0f / 255, 207.0f / 255, 201.0f / 255},
                std::array<float, 3>{255.0f / 255, 190.0f / 255, 122.0f / 255},
                std::array<float, 3>{250.0f / 255, 127.0f / 255, 111.0f / 255},
                std::array<float, 3>{130.0f / 255, 176.0f / 255, 210.0f / 255},
                std::array<float, 3>{190.0f / 255, 184.0f / 255, 220.0f / 255},
                std::array<float, 3>{40.0f / 255, 120.0f / 255, 181.0f / 255},
                std::array<float, 3>{248.0f / 255, 172.0f / 255, 140.0f / 255},
                std::array<float, 3>{255.0f / 255, 136.0f / 255, 132.0f / 255},
                std::array<float, 3>{84.0f / 255, 179.0f / 255, 69.0f / 255},
                std::array<float, 3>{137.0f / 255, 131.0f / 255, 191.0f / 255}
            };
            std::ofstream ofs(path);
            for(size_t i = 0; i < _points.size(); i++)
            {
                auto c = COLORS[_labels[i] % COLORS.size()];
                if(_labels[i] < 0 || _labels[i] >= 49)
                {
                    c[0] = 1.0f;
                    c[1] = 0.f;
                    c[2] = 0.f;
                }
                ofs << "v " << _points[i].x() << ' ' << _points[i].y() << ' ' << _points[i].z() << ' ' << c[0] << ' ' << c[1] << ' ' << c[2] << "\n";
            }
            for(auto& [label, c] : _centroids)
            {
                ofs << "v " << c.x() << ' ' << c.y() << ' ' << c.z() << " 1 0 0\n";
            }
            for(auto& [label, up] : _upwards)
            {
                for(int i = 0; i < 100; i++)
                {
                    Point_3 p = _centroids.at(label) + up * 0.1f * (i + 1);
                    ofs << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << " 0 0 1\n";
                }
            }
            ofs << std::endl;
        }
        void ForEachSegment(const std::function<void(const Point_3& p0, const Point_3& p1, int l0, int l1, size_t idx)>& func) const
        {
            for(size_t i = 0; i < _points.size(); i++)
            {
                func(_points[i], _points[(i + 1) % _points.size()], _labels[i], _labels[(i + 1) % _points.size()], i);
            }
        }
        CurveIterator<Curve<Kernel>> CreateIterator(size_t pos)
        {
            return CurveIterator<Curve<Kernel>>(this, pos);
        }
        CurveIterator<const Curve<Kernel>> CreateIterator(size_t pos) const
        {
            return CurveIterator<const Curve<Kernel>>(this, pos);
        }
        CurveIterator<const Curve<Kernel>> CreateConstIterator(size_t pos) const
        {
            return CurveIterator<const Curve<Kernel>>(this, pos);
        }
        template <typename AABBTree>
        void FixShape(int label, const AABBTree& guide_mesh, double fix_factor)
        {
            if(_points.size() <= 2)
            {
                return;
            }
            std::vector<std::vector<iterator>> segs;
            iterator start_pos = CreateIterator(0);
            int cnt = 0;
            while(start_pos.Label() != label)
            {
                start_pos++;
                if(++cnt >= _points.size())
                {
                    return;
                }
            }
            cnt = 0;
            while(start_pos.Label() == label)
            {
                start_pos--;
                if(++cnt >= _points.size())
                {
                    return;
                }
            }
            start_pos++;
            if(start_pos.Label() == label)
            {
                segs.emplace_back();
            }
            else
            {
                return;
            }
            iterator end_pos = start_pos;
            iterator it = start_pos;
            do
            {
                if(it.Label() == label)
                {
                    segs.back().push_back(it);
                }
                else
                {
                    if(!segs.back().empty())
                    {
                        segs.emplace_back();
                    }
                }
                it++;
            }while(it != end_pos);
            if(segs.back().empty())
            {
                segs.pop_back();
            }
            if(segs.size() == 1)
            {
                return;
            }
            auto tooth_dir = _upwards.at(label);
            typename Kernel::Line_3 tooth_axis(_centroids[label], _upwards[label]);
            //EasyOBJ file("fixshape" + std::to_string(label) + ".obj");
            for(auto& seg : segs)
            {
                iterator mid = *std::min_element(seg.begin(), seg.end(), [&](auto& lh, auto& rh) {
                    double dotl = (tooth_axis.projection(lh.Point()) - _centroids[label]) * tooth_dir;
                    double dotr = (tooth_axis.projection(rh.Point()) - _centroids[label]) * tooth_dir;
                    return dotl < dotr;
                });
                typename Kernel::Line_3 line1(seg.front().Point(), mid.Point());
                typename Kernel::Line_3 line2(mid.Point(), seg.back().Point());
                typename Kernel::Plane_3 plane1(mid.Point(), seg.front().Point(), mid.Point() + CGAL::cross_product(line1.to_vector(), tooth_dir));
                typename Kernel::Plane_3 plane2(seg.back().Point(), mid.Point(), mid.Point() + CGAL::cross_product(line2.to_vector(), tooth_dir));
                double len1 = CGAL::squared_distance(mid.Point(), seg.front().Point());
                double len2 = CGAL::squared_distance(mid.Point(), seg.back().Point());
                auto it = seg.front();
                //file.AddV(mid.Point(), 0, 1, 0);
                while(it != mid)
                {
                    bool distance_test = CGAL::squared_distance(line1, it.Point()) >= len1 / 25.0;
                    bool angle_test = CGAL::approximate_angle(line1.to_vector(), (it.Point() - (it - 1).Point())) > 30.0;
                    bool plane_test = plane1.has_on_positive_side(it.Point()) && CGAL::approximate_angle(it.Point() - (it - 1).Point(), plane1.projection(it.Point()) - plane1.projection((it - 1).Point())) > 30.0;
                    if(plane_test || angle_test || distance_test)
                    {
                        auto break_start = it - 1;
                        auto break_end = it;
                        while((plane_test || angle_test || distance_test) && break_end != mid)
                        {
                            break_end++;
                            distance_test = CGAL::squared_distance(line1, break_end.Point()) >= len1 / 25.0;
                            angle_test = CGAL::approximate_angle(line1.to_vector(), (break_end.Point() - break_start.Point())) > 30.0;
                            plane_test = plane1.has_on_positive_side(it.Point()) && CGAL::approximate_angle(break_end.Point() - break_start.Point(), plane1.projection(break_end.Point()) - plane1.projection(break_start.Point())) > 30.0;
                        }
                        double diff = (double)(break_end - break_start);
                        while(it != break_end)
                        {
                            it.Point() = break_start.Point() + (double)(it - break_start) / diff * (break_end.Point() - break_start.Point());
                            it++;
                        }
                        // do some smooth
                        for(int j = 0; j < 5; j++)
                        {
                            it = break_start + 1;
                            while(it != break_end)
                            {
                                it.Point() = CGAL::midpoint((it - 1).Point(), (it + 1).Point());
                                it++;
                            }
                        }
                        it = break_start + 1;
                        while(it != break_end)
                        {
                            it.Point() = guide_mesh.closest_point(it.Point());
                            //file.AddV(it.Point(), 1, 0, 0);
                            it++;
                        }
                    }
                    else
                    {
                        //file.AddV(it.Point(), 0, 0, 1);
                        it++;
                    }
                }
                it = seg.back();
                while(it != mid)
                {
                    bool distance_test = CGAL::squared_distance(line2, it.Point()) >= len2 / 25.0;
                    bool angle_test = CGAL::approximate_angle(-line2.to_vector(), (it.Point() - (it - 1).Point())) > 30.0;
                    bool plane_test = plane1.has_on_positive_side(it.Point()) && CGAL::approximate_angle(it.Point() - (it + 1).Point(), plane2.projection(it.Point()) - plane2.projection((it + 1).Point())) > 30.0;
                    if(plane_test || angle_test || distance_test)
                    {
                        auto break_start = it + 1;
                        auto break_end = it;
                        while((plane_test || angle_test || distance_test) && break_end != mid)
                        {
                            break_end--;
                            distance_test = CGAL::squared_distance(line2, break_end.Point()) >= len2 / 25.0;
                            angle_test = CGAL::approximate_angle(-line2.to_vector(), (break_end.Point() - break_start.Point())) > 30.0;
                            plane_test = plane1.has_on_positive_side(it.Point()) && CGAL::approximate_angle(break_end.Point() - break_start.Point(), plane2.projection(break_end.Point()) - plane2.projection(break_start.Point())) > 30.0;
                        }
                        double diff = (double)(break_end - break_start);
                        while(it != break_end)
                        {
                            it.Point() = break_start.Point() + (double)(it - break_start) / diff * (break_end.Point() - break_start.Point());
                            it--;
                        }
                        // do some smooth
                        for(int j = 0; j < 5; j++)
                        {
                            it = break_start - 1;
                            while(it != break_end)
                            {
                                it.Point() = CGAL::midpoint((it - 1).Point(), (it + 1).Point());
                                it--;
                            }
                        }
                        it = break_start - 1;
                        while(it != break_end)
                        {
                            it.Point() = guide_mesh.closest_point(it.Point());
                            //file.AddV(it.Point(), 1, 0, 0);
                            it--;
                        }
                    }
                    else
                    {
                        //file.AddV(it.Point(), 0, 0, 1);
                        it--;
                    }
                }

                for(int j = 0; j < 10; j++)
                {
                    for(int k = 1; k < seg.size() - 1; k++)
                    {
                        if(seg[k] == mid)
                        {
                            continue;
                        }
                        seg[k].Point() = CGAL::midpoint((seg[k] - 1).Point(), (seg[k] + 1).Point());
                    }
                }
                for(int k = 1; k < seg.size() - 1; k++)
                {
                    seg[k].Point() = guide_mesh.closest_point(seg[k].Point());
                }
            }
        }
        template <typename AABBTree>
        void FixAllCurve(const AABBTree& guide_mesh, double fix_factor)
        {
            printf("Adjusting curve...\n");
            std::unordered_set<int> all_labels(_labels.begin(), _labels.end());
            for(int label : all_labels)
            {
                if(label % 10 <= 4)
                {
                    FixShape(label, guide_mesh, fix_factor);
                }
            }
        }
    protected:
        std::vector<Point_3> _points;
        std::vector<int> _labels;
        std::unordered_map<int, Point_3> _centroids;
        std::unordered_map<int, Vector_3> _upwards;
        std::vector<int> _ordered_labels;
    };

    template <typename Curve>
    class CurveIterator
    {
    public:
        using Kernel = typename Curve::K;

        CurveIterator(Curve* curve) : _curve(curve), _pos(0u){}
        CurveIterator(Curve* curve, size_t pos) : _curve(curve), _pos(pos) {}
        CurveIterator(const CurveIterator<Curve>& rh) : _curve(rh._curve), _pos(rh._pos) {}
        CurveIterator& operator=(const CurveIterator<Curve>& rh)
        {
            _curve = rh._curve;
            _pos = rh._pos;
            return *this;
        }
        CurveIterator& operator++()
        {
            _pos = (_pos == _curve->size() - 1) ? 0 : _pos + 1;
            return *this;
        }
        CurveIterator operator++(int)
        {
            CurveIterator ret = *this;
            _pos = (_pos == _curve->size() - 1) ? 0 : _pos + 1;
            return ret;
        }
        CurveIterator& operator--()
        {
            _pos = _pos == 0 ? _curve->size() - 1 : _pos - 1;
            return *this;
        }
        CurveIterator operator--(int)
        {
            CurveIterator ret = *this;
            _pos = _pos == 0 ? _curve->size() - 1 : _pos - 1;
            return ret;
        }
        CurveIterator& operator+=(int diff)
        {
            if(diff >= 0)
            {
                _pos = (_pos + diff) % _curve->size();
            }
            else
            {
                _pos = (_pos + _curve->size() + diff) % _curve->size();
            }
            return *this;
        }
        bool operator==(const CurveIterator<Curve>& it) const
        {
            return _curve == it._curve && _pos == it._pos;
        }
        bool operator!=(const CurveIterator<Curve>& it) const
        {
            return _curve != it._curve || _pos != it._pos;
        }
        size_t Idx() const
        {
            return _pos;
        }
        std::conditional_t<std::is_const_v<Curve>, const typename Kernel::Point_3&, typename Kernel::Point_3&> Point()
        {
            return (*_curve)[_pos];
        }
        int Label() const
        {
            return _curve->Label(_pos);
        }
        typename Kernel::Segment_3 Segment() const
        {
            size_t i0 = _pos;
            size_t i1 = (_pos == _curve->size() - 1) ? 0 : _pos + 1;
            return typename Kernel::Segment_3((*_curve)[i0], (*_curve)[i1]);
        }
        friend CurveIterator operator+(const CurveIterator& it, int diff)
        {
            CurveIterator ret = it;
            if(diff >= 0)
            {
                ret._pos = (ret._pos + diff) % ret._curve->size();
            }
            else
            {
                ret._pos = (ret._pos + ret._curve->size() + diff) % ret._curve->size();
            }
            return ret;
        }
        friend CurveIterator operator+(int diff, const CurveIterator& it)
        {
            CurveIterator ret = it;
            if(diff >= 0)
            {
                ret._pos = (ret._pos + diff) % ret._curve->size();
            }
            else
            {
                ret._pos = (ret._pos + ret._curve->size() + diff) % ret._curve->size();
            }
            return ret;
        }
        friend CurveIterator operator-(const CurveIterator& it, int diff)
        {
            return it + (-diff);
        }
        friend int operator-(const CurveIterator& lh, const CurveIterator& rh)
        {
            if(lh.Idx() >= rh.Idx())
            {
                return static_cast<int>(lh.Idx()) - static_cast<int>(rh.Idx());
            }
            else
            {
                return (static_cast<int>(lh.Idx()) + static_cast<int>(lh._curve->size()) - static_cast<int>(rh.Idx())) % static_cast<int>(lh._curve->size());
            }
        }

    protected:
        Curve* _curve;
        size_t _pos;
    };

    template <typename Kernel, typename AABBTree>
    Curve<Kernel> Merge(const Curve<Kernel> &curve0, const Curve<Kernel> &curve1, const AABBTree& guide_mesh )
    {
        // if(curve0.MaxLabel() > curve1.MinLabel())
        // {
        //     throw AlgError("Currently the algorithm do not allow overlapping labels."
        //      "This may happen because of bad segmentation or teeth of unusual shape."
        //      "labels of curve 0: (" + std::to_string(curve0.MinLabel()) + ", " + std::to_string(curve0.MaxLabel()) + ")"
        //       ", labels of curve 1: (" + std::to_string(curve1.MinLabel()) + ", " + std::to_string(curve1.MaxLabel()) + ")");
        // }
        printf("Merging curve0 of label (%d-%d) and curve1 of label (%d-%d)\n", curve0.MinLabel(), curve0.MaxLabel(), curve1.MinLabel(), curve1.MaxLabel());
        typename Kernel::Point_3 c0 = curve0.GetCentroidOfLabel(curve0.MaxLabel());
        typename Kernel::Point_3 c1 = curve1.GetCentroidOfLabel(curve1.MinLabel());
        typename Kernel::Vector_3 up0 = curve0.GetUpwardOfLabel(curve0.MaxLabel());
        typename Kernel::Vector_3 up1 = curve1.GetUpwardOfLabel(curve1.MinLabel());
        typename Kernel::Vector_3 mid_up = (up0 + up1) * 0.5;
        typename Kernel::Plane_3 connecting_plane(c0, c1, c0 + mid_up);
        typename Kernel::Line_3 line(c0, c1);
        typename Kernel::Vector_3 line_dir = c1 - c0;

        std::pair<size_t, size_t> closest_pair;
        // {
        //     double min_dist = std::numeric_limits<double>::max();
        //     curve0.ForEachSegment([&](const typename Kernel::Point_3& p0, const typename Kernel::Point_3& p1, int l0, int l1, size_t idx)
        //     {
        //         if(l0 == curve0.MaxLabel() && CGAL::angle(p1 - p0, connecting_plane.orthogonal_vector()) == CGAL::Angle::ACUTE)
        //         {
        //             double d = CGAL::squared_distance(connecting_plane.projection(p0), p0);
        //             if(d < min_dist)
        //             {
        //                 closest_pair.first = idx;
        //                 min_dist = d;
        //             }
        //         }
        //     });
        // }
        // {
        //     double min_dist = std::numeric_limits<double>::max();
        //     curve1.ForEachSegment([&](const typename Kernel::Point_3& p0, const typename Kernel::Point_3& p1, int l0, int l1, size_t idx)
        //     {
        //         if(l0 == curve1.MinLabel() && CGAL::angle(p1 - p0, connecting_plane.orthogonal_vector()) == CGAL::Angle::OBTUSE)
        //         {
        //             double d = CGAL::squared_distance(connecting_plane.projection(p0), p0);
        //             if(d < min_dist)
        //             {
        //                 closest_pair.second = idx;
        //                 min_dist = d;
        //             }
        //         }
        //     });
        // }

        bool found_fail1 = true;
        bool found_fail2 = true;
        auto closest_its = std::make_pair(curve0.CreateIterator(0), curve1.CreateIterator(0));
        {
            double min_dist = std::numeric_limits<double>::max();
            auto it = curve0.CreateIterator(0);
            auto end = it;
            do
            {
                if(it.Label() == curve0.MaxLabel())
                {
                    found_fail1 = false;
                    double d = CGAL::squared_distance(c1, it.Point());
                    if(d < min_dist)
                    {
                        closest_its.first = it;
                        min_dist = d;
                    }
                }
                it++;
            } while (it != end);
        }
        {
            double min_dist = std::numeric_limits<double>::max();
            auto it = curve1.CreateIterator(0);
            auto end = it;
            do
            {
                if(it.Label() == curve1.MinLabel())
                {
                    found_fail2 = false;
                    double d = CGAL::squared_distance(c0, it.Point());
                    if(d < min_dist)
                    {
                        closest_its.second = it;
                        min_dist = d;
                    }
                }
                it++;
            } while (it != end);
        }
        
        if(found_fail1 || found_fail2)
        {
            printf("Warning: failed to found closest merge pair. It may caused by bad geom or unknown bug. A trival searching method will be applied.\n");
            auto it0 = curve0.CreateIterator(0);
            auto end0 = it0;
            double min_dist = std::numeric_limits<double>::max();
            do
            {
                auto it1 = curve1.CreateIterator(0);
                auto end1 = it1;
                if(it0.Label() != curve0.MaxLabel())
                {
                    continue;
                }
                do
                {
                    if(it1.Label() == curve1.MinLabel())
                    {
                        double curr_dist = CGAL::squared_distance(it0.Point(), it1.Point());
                        if(curr_dist < min_dist)
                        {
                            closest_its.first = it0;
                            closest_its.second = it1;
                            min_dist = curr_dist;
                        }
                    }
                    it1++;
                } while (it1 != end1);
                it0++;
            } while (it0 != end0);
        }
        
        closest_pair.first = closest_its.first.Idx();
        closest_pair.second = closest_its.second.Idx();
        
        // TODO: Change following code to using iterators.
        double w0 = 1.0f;
        double w1 = 0.0f;

        auto end_it = closest_its;
        auto start_it = closest_its;
        double max_score = -999.0;
        auto it = closest_its.first;
        for (size_t i = 0; i < curve0.size(); i++, it++)
        {
            if(it.Label() != closest_its.first.Label())
            {
                break;
            }
            typename Kernel::Vector_3 vec = it.Segment().to_vector();
            double distance_score = std::sqrt(CGAL::squared_distance(connecting_plane.projection(it.Point()), it.Point()));
            double direction_score = CGAL::approximate_angle(line_dir, vec);
            double score = w0 * distance_score + w1 * direction_score;
            if(connecting_plane.has_on_positive_side(it.Point()) && score >= max_score)
            {
                max_score = score;
                end_it.first = it;
            }
            else
            {
                //break;
            }
            it++;
        }
        max_score = -999.0;
        it = closest_its.second;
        for (size_t i = 0; i < curve1.size(); i++)
        {
            if(it.Label() != closest_its.second.Label())
            {
                break;
            }
            typename Kernel::Vector_3 vec = it.Segment().to_vector();
            double distance_score = std::sqrt(CGAL::squared_distance(connecting_plane.projection(it.Point()), it.Point()));
            double direction_score = 180.0 - CGAL::approximate_angle(line_dir, vec);
            double score = w0 * distance_score + w1 * direction_score;
            if(connecting_plane.has_on_negative_side(it.Point()) && score >= max_score)
            {
                max_score = score;
                end_it.second = it;
            }
            else
            {
                //break;
            }
            it++;
        }
        max_score = -999.0;
        it = closest_its.first;
        for (size_t i = 0; i < curve0.size(); i++)
        {
            if(it.Label() != closest_its.first.Label())
            {
                break;
            }
            typename Kernel::Vector_3 vec = it.Segment().to_vector();
            double distance_score = std::sqrt(CGAL::squared_distance(connecting_plane.projection(it.Point()), it.Point()));
            double direction_score = CGAL::approximate_angle(line_dir, vec);
            double score = w0 * distance_score + w1 * distance_score;
            if(connecting_plane.has_on_negative_side(it.Point()) && score >= max_score)
            {
                max_score = score;
                start_it.first = it;
            }
            else
            {
                //break;
            }
            it--;
        }

        max_score = -999.0;
        it = closest_its.second;
        for (size_t i = 0; i < curve1.size(); i++)
        {
            if(it.Label() != closest_its.second.Label())
            {
                break;
            }
            typename Kernel::Vector_3 vec = it.Segment().to_vector();
            double direction_score = 180.0 - CGAL::approximate_angle(line_dir, vec);
            double distance_score = std::sqrt(CGAL::squared_distance(connecting_plane.projection(it.Point()), it.Point()));
            double score = w0 * distance_score + w1 * distance_score;
            if(connecting_plane.has_on_positive_side(it.Point()) && score >= max_score)
            {
                max_score = score;
                start_it.second = it;
            }
            else
            {
                //break;
            }
            it--;
        }

        // if(start_pair.first < closest_pair.first)
        // {
        //     start_pair.first = static_cast<size_t>(start_pair.first * 0.7 + closest_pair.first * 0.3);
        // }
        // else
        // {
        //     start_pair.first = static_cast<size_t>(start_pair.first * 0.7 + (closest_pair.first + curve0.size()) * 0.3) % curve0.size();
        // }
        // if(end_pair.first > closest_pair.first)
        // {
        //     end_pair.first = static_cast<size_t>(end_pair.first * 0.7 + closest_pair.first * 0.3);
        // }
        // else
        // {
        //     end_pair.first = static_cast<size_t>(closest_pair.first * 0.3 + (curve0.size() + end_pair.first) * 0.7) % curve0.size();
        // }
        // if(start_pair.second < closest_pair.second)
        // {
        //     start_pair.second = static_cast<size_t>(start_pair.second * 0.7 + closest_pair.second * 0.3);
        // }
        // else
        // {
        //     start_pair.second = static_cast<size_t>(start_pair.second * 0.7 + (closest_pair.second + curve1.size()) * 0.3) % curve1.size();
        // }
        // if(end_pair.second > closest_pair.second)
        // {
        //     end_pair.second = static_cast<size_t>(end_pair.second * 0.7 + closest_pair.second * 0.3);
        // }
        // else
        // {
        //     end_pair.second = static_cast<size_t>(closest_pair.second * 0.3 + (curve1.size() + end_pair.second) * 0.7) % curve1.size();
        // }

        // std::ofstream ofs("./merge" + std::to_string(curve0.size()) + ".obj");
        // for (size_t i = 0; i < curve0.size(); i++)
        // {
        //     float r = 0;
        //     float g = 0;
        //     float b = 0;
        //     if(i == closest_its.first.Idx())
        //     {
        //         r = 1.f;
        //     }
        //     if(i == start_it.first.Idx())
        //     {
        //         g = 1.f;
        //     }
        //     if(i == end_it.first.Idx())
        //     {
        //         b = 1.f;
        //     }

        //     ofs << "v " << curve0[i].x() << ' ' << curve0[i].y() << ' ' << curve0[i].z() << ' ' << r << ' ' << g << ' ' << b << '\n';
        // }
        // for(size_t i = 0; i < curve1.size(); i++)
        // {
        //     float r = 0;
        //     float g = 0;
        //     float b = 0;
        //     if(i == closest_its.second.Idx())
        //     {
        //         r = 1.f;
        //     }
        //     if(i == start_it.second.Idx())
        //     {
        //         g = 1.f;
        //     }
        //     if(i == end_it.second.Idx())
        //     {
        //         b = 1.f;
        //     }
        //     ofs << "v " << curve1[i].x() << ' ' << curve1[i].y() << ' ' << curve1[i].z() << ' ' << r << ' ' << g << ' ' << b << '\n';
        // }
        // ofs.close();
        
        Curve<Kernel> midcurve1;
        {
            double error = 0.0;
            int max_try = std::min(closest_its.first - start_it.first, end_it.second - closest_its.second) - 1;
            do
            {
                typename Kernel::Point_3 pos0 = start_it.first.Point();
                typename Kernel::Point_3 pos1 = end_it.second.Point();
                typename Kernel::Point_3 mid = CGAL::midpoint(pos0, pos1);
                mid = CGAL::midpoint(connecting_plane.projection(mid), mid);
                for(int i = 1; i < 40; i++)
                {
                    double t = (double)i / 40.0;
                    typename Kernel::Vector_3 p = (1.0 - t) * (1.0 - t) * (pos0 - CGAL::ORIGIN) + 2.0 * t * (1.0 - t) * (mid - CGAL::ORIGIN) + t * t * (pos1 - CGAL::ORIGIN);
                    midcurve1.AddPoint(CGAL::ORIGIN + p, start_it.first.Label());
                }
                error = 0.0;
                for(size_t i = 0; i < midcurve1.size(); i++)
                {
                    if(guide_mesh.closest_point_and_primitive(midcurve1[i]).second->_label != 0)
                    {
                        error += 1.0;
                    }
                }
                error /= midcurve1.size();
                if(error < 0.3)
                {
                    break;
                }
                else
                {
                    //printf("Merge error: %f\n", error);
                    start_it.first++;
                    end_it.second--;
                    max_try--;
                    midcurve1 = Curve<Kernel>();
                }
            } while(max_try > 0);
            
        }
        Curve<Kernel> midcurve2;
        {
            double error = 0.0;
            int max_try = std::min(closest_its.second - start_it.second, end_it.first - closest_its.first) - 1;
            do
            {
                typename Kernel::Point_3 pos0 = start_it.second.Point();
                typename Kernel::Point_3 pos1 = end_it.first.Point();
                typename Kernel::Point_3 mid = CGAL::midpoint(pos0, pos1);
                mid = CGAL::midpoint(connecting_plane.projection(mid), mid);
                for(int i = 1; i < 40; i++)
                {
                    double t = (double)i / 40.0;
                    typename Kernel::Vector_3 p = (1.0 - t) * (1.0 - t) * (pos0 - CGAL::ORIGIN) + 2.0 * t * (1.0 - t) * (mid - CGAL::ORIGIN) + t * t * (pos1 - CGAL::ORIGIN);
                    midcurve2.AddPoint(CGAL::ORIGIN + p, end_it.first.Label());
                }
                error = 0.0;
                for(size_t i = 0; i < midcurve2.size(); i++)
                {
                    if(guide_mesh.closest_point_and_primitive(midcurve2[i]).second->_label != 0)
                    {
                        error += 1.0;
                    }
                }
                error /= midcurve2.size();
                if(error < 0.3)
                {
                    break;
                }
                else
                {
                    printf("Merge error: %f\n", error);
                    start_it.second++;
                    end_it.first--;
                    max_try--;
                    midcurve2 = Curve<Kernel>();
                }
            }while(max_try > 0);
        }
        Curve<Kernel> subcurve1 = curve0.GetSubCurve(end_it.first, start_it.first);
        Curve<Kernel> subcurve2 = curve1.GetSubCurve(end_it.second, start_it.second);

        if(midcurve1.size() >= 3)
        {
            for(int j = 0; j < 0; j++)
            {
                midcurve1[0] = CGAL::midpoint(midcurve1[1], subcurve1[subcurve1.size() - 1]);
                midcurve1[midcurve1.size() - 1] = CGAL::midpoint(midcurve1[midcurve1.size() - 2], subcurve2[0]);
                for(int k = 1; k < midcurve1.size() - 1; k++)
                {
                    midcurve1[k] = CGAL::midpoint(midcurve1[k - 1], midcurve1[k + 1]);
                }
            }
        }
        if(midcurve2.size() >= 3)
        {
            for(int j = 0; j < 0; j++)
            {
                midcurve2[0] = CGAL::midpoint(midcurve2[1], subcurve2[subcurve2.size() - 1]);
                midcurve2[midcurve2.size() - 1] = CGAL::midpoint(midcurve2[midcurve2.size() - 2], subcurve1[0]);
                for(int k = 1; k < midcurve1.size() - 1; k++)
                {
                    midcurve2[k] = CGAL::midpoint(midcurve2[k - 1], midcurve2[k + 1]);
                }
            }
        }

        Curve<Kernel> result_curve = subcurve1;
        result_curve.InsertAt(midcurve1, result_curve.size());
        result_curve.InsertAt(subcurve2, result_curve.size());
        result_curve.InsertAt(midcurve2, result_curve.size());

        //std::ofstream ofs("./merge.obj");
        static const std::array<std::array<float, 3>, 10> COLORS = {
            std::array<float, 3>{142.0f / 255, 207.0f / 255, 201.0f / 255},
            std::array<float, 3>{255.0f / 255, 190.0f / 255, 122.0f / 255},
            std::array<float, 3>{250.0f / 255, 127.0f / 255, 111.0f / 255},
            std::array<float, 3>{130.0f / 255, 176.0f / 255, 210.0f / 255},
            std::array<float, 3>{190.0f / 255, 184.0f / 255, 220.0f / 255},
            std::array<float, 3>{40.0f / 255, 120.0f / 255, 181.0f / 255},
            std::array<float, 3>{248.0f / 255, 172.0f / 255, 140.0f / 255},
            std::array<float, 3>{255.0f / 255, 136.0f / 255, 132.0f / 255},
            std::array<float, 3>{84.0f / 255, 179.0f / 255, 69.0f / 255},
            std::array<float, 3>{137.0f / 255, 131.0f / 255, 191.0f / 255}
        };
        // for (size_t i = 0; i < result_curve.size(); i++)
        // {
        //     //std::array<float, 3> c = COLORS[result_curve.Label(i) % COLORS.size()];
        //     auto p0 = i == 0 ? result_curve[result_curve.size() - 1] : result_curve[i - 1];
        //     auto p1 = result_curve[i];
        //     auto p2 = result_curve[(i + 1) % result_curve.size()];
        //     double dx = std::sqrt(CGAL::squared_distance(CGAL::midpoint(p0, p2), p1));
        //     auto c = tinycolormap::GetColor(dx * 100);
        //     ofs << "v " << result_curve[i].x() << ' ' << result_curve[i].y() << ' ' << result_curve[i].z() << ' ' << c[0] << ' ' << c[1] << ' ' << c[2] << '\n';
        // }
        result_curve.UpdateData();
        return result_curve;
    }
}


bool GumTrimLine( std::string input_file, std::string label_file, std::string frame_file, std::string output_file, int smooth, double fix_factor );

#endif