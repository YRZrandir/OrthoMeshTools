#ifndef EASY_OBJ_H
#define EASY_OBJ_H
#include <fstream>
#include <string>

class EasyOBJ
{
public:
    EasyOBJ(std::string filename)
        :_ofs(filename)
    {
    }

    size_t AddV(double x, double y, double z, double r, double g, double b)
    {
        _ofs << "v " << x << ' ' << y << ' ' << z << ' ' << r << ' ' << g << ' ' << b << '\n';
        _vcount++;
        return _vcount;
    }

    size_t AddV(double x, double y, double z)
    {
        _ofs << "v " << x << ' ' << y << ' ' << z << '\n';
        _vcount++;
        return _vcount;
    }

    template <typename Point>
    size_t AddV(const Point& p)
    {
        _ofs << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
        _vcount++;
        return _vcount;
    }

    template <typename Point>
    size_t AddV(const Point& p, double r, double g, double b)
    {
        _ofs << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << ' ' << r << ' ' << g << ' ' << b << '\n';
        _vcount++;
        return _vcount;
    }

    void AddL(size_t i0, size_t i1)
    {
        _ofs << "l " << i0 << ' ' << i1 << '\n';
    }

    template <typename Point>
    void AddL(const Point& p0, const Point& p1)
    {
        size_t i0 = AddV(p0);
        size_t i1 = AddV(p1);
        AddL(i0, i1);
    }

    void AddF(size_t i0, size_t i1, size_t i2)
    {
        _ofs << "f " << i0 << ' ' << i1 << ' ' << i2 << '\n';
    }

protected:
    std::ofstream _ofs;
    size_t _vcount = 0;
};

#endif