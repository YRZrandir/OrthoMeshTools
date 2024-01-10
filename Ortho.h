#ifndef ORTHO_H
#define ORTHO_H
#include <array>
#include <exception>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <nlohmann/json.hpp>

class IOError : public std::runtime_error
{
public:
    IOError( const std::string& what_arg ) : std::runtime_error(what_arg) {}
    IOError( const char* what_arg ) : std::runtime_error(what_arg) {}
};

class MeshError : public std::runtime_error
{
public:
    MeshError( const std::string& what_arg ) : std::runtime_error(what_arg) {}
    MeshError( const char* what_arg ) : std::runtime_error(what_arg) {}
};

class AlgError : public std::runtime_error
{
public:
    AlgError( const std::string& what_arg ) : std::runtime_error(what_arg) {}
    AlgError( const char* what_arg ) : std::runtime_error(what_arg) {}
};

template <typename Kernel>
struct Frame
{
    typename Kernel::Vector_3 right;
    typename Kernel::Vector_3 front;
    typename Kernel::Vector_3 up;
    typename Kernel::Point_3 pos;
};

template <typename Kernel>
class CrownFrames
{
public:
    CrownFrames(const std::string& path )
    {
        using namespace nlohmann;
        std::ifstream label_ifs( path );
        if(label_ifs.fail())
        {
            throw IOError("Cannot open file: " + path);
        }
        json js = json::parse( label_ifs );
        for(int label = 11; label < 50; label++)
        {
            if(js.find(std::to_string(label)) != js.end())
            {
                std::vector<std::vector<double>> frame_data = js[std::to_string(label)].get<std::vector<std::vector<double>>>();
                if(frame_data.size() != 4 || frame_data[0].size() != 3 || frame_data[1].size() != 3 || frame_data[2].size() != 3)
                {
                    throw IOError("crown frame json format is invalid.\n");
                }
                Frame<Kernel> frame;
                frame.right = typename Kernel::Vector_3(frame_data[0][0], frame_data[0][1], frame_data[0][2]);
                frame.front = typename Kernel::Vector_3(frame_data[1][0], frame_data[1][1], frame_data[1][2]);
                frame.up = typename Kernel::Vector_3(frame_data[2][0], frame_data[2][1], frame_data[2][2]);
                frame.pos = typename Kernel::Point_3(frame_data[3][0], frame_data[3][1], frame_data[3][2]);
                _frames[label] = frame;
            }
        }
    }

    const Frame<Kernel>& GetFrame(int label) const
    {
        if(_frames.count(label) == 0)
        {
            throw AlgError("Cannot find frame data of label " + std::to_string(label));
        }
        return _frames.at(label);
    }

    std::unordered_map<int, Frame<Kernel>>::const_iterator Begin() const
    {
        return _frames.begin();
    }

    std::unordered_map<int, Frame<Kernel>>::const_iterator End() const
    {
        return _frames.end();
    }

    const std::unordered_map<int, Frame<Kernel>>& Frames() const
    {
        return _frames;
    }

protected:
    std::unordered_map<int, Frame<Kernel>> _frames;
};

constexpr std::array<float, 3> LabelColorMap(int label)
{
    constexpr std::array<std::array<float, 3>, 10> COLORS = {
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
    return COLORS[label % COLORS.size()];
}
#endif