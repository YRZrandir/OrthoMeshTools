#ifndef ORTHO_H
#define ORTHO_H
#include <array>
#include <exception>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <nlohmann/json.hpp>
#include <Eigen/Eigen>

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

    Frame() = default;
    Frame(const typename Kernel::Aff_transformation_3& t)
    {
        right = typename Kernel::Vector_3(t.m(0, 0), t.m(1, 0), t.m(2, 0));
        front = typename Kernel::Vector_3(t.m(0, 1), t.m(1, 1), t.m(2, 1));
        up = typename Kernel::Vector_3(t.m(0, 2), t.m(1, 2), t.m(2, 2));
        pos = typename Kernel::Point_3(t.m(0, 3), t.m(1, 3), t.m(2, 3));
    }

    typename Kernel::Aff_transformation_3 WorldToLocal() const
    {
        return LocalToWorld().inverse();
    }

    typename Kernel::Aff_transformation_3 LocalToWorld() const
    {
        return typename Kernel::Aff_transformation_3(
            right.x(), front.x(), up.x(), pos.x(),
            right.y(), front.y(), up.y(), pos.y(),
            right.z(), front.z(), up.z(), pos.z(),
            1.0
        );
    }
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

    CrownFrames() = default;

    void Insert(int label, const Frame<Kernel>& frame)
    {
        _frames.insert({label, frame});
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

template <typename Scalar>
class CBCTRegis
{
public:
    CBCTRegis( const std::string& path )
    {
        nlohmann::json json = nlohmann::json::parse(std::ifstream(path));
        for(int label = 10; label < 50; label++)
        {
            std::string ios_to_cbct_name = "ios_to_cbct_matrix_" + std::to_string(label);
            std::string cbct_to_ios_name = "cbct_to_ios_matrix_" + std::to_string(label);
            for(int i = 0; i < json.size(); i++)
            {
                if(json[i].find(ios_to_cbct_name) != json[i].end())
                {
                    Eigen::Matrix<Scalar, 4, 4> ios_to_cbct;
                    std::vector<std::vector<double>> data = json[i][ios_to_cbct_name].get<std::vector<std::vector<double>>>();
                    for(int r = 0; r < 4; r++)
                        for(int c = 0; c < 4; c++)
                            ios_to_cbct(r, c) = data[r][c];
                    _mats[label].first = Eigen::Transform<Scalar, 3, Eigen::Affine>(ios_to_cbct);
                }
                if(json[i].find(cbct_to_ios_name) != json[i].end())
                {
                    Eigen::Matrix<Scalar, 4, 4> cbct_to_ios;
                    std::vector<std::vector<double>> data = json[i][cbct_to_ios_name].get<std::vector<std::vector<double>>>();
                    for(int r = 0; r < 4; r++)
                        for(int c = 0; c < 4; c++)
                            cbct_to_ios(r, c) = data[r][c];
                    _mats[label].second = Eigen::Transform<Scalar, 3, Eigen::Affine>(cbct_to_ios);
                }
            }
        }
        printf("CBCT: %zd\n", _mats.size());
    }

    Eigen::Transform<Scalar, 3, Eigen::Affine> IOS_to_CBCT(int label) const
    {
        if(_mats.count(label) != 0)
        {
            return _mats.at(label).first;
        }
        else
        {
            throw AlgError("No cbct info: " + std::to_string(label));
            //return Eigen::Transform<Scalar, 3, Eigen::Affine>::Identity();
        }
    }
    Eigen::Transform<Scalar, 3, Eigen::Affine> CBCT_to_IOS(int label) const
    {
        if(_mats.count(label) != 0)
        {
            return _mats.at(label).second;
        }
        else
        {
            throw AlgError("No cbct info: " + std::to_string(label));
            //return Eigen::Transform<Scalar, 3, Eigen::Affine>::Identity();
        }
    }
protected:
    std::unordered_map<int, std::pair<Eigen::Transform<Scalar, 3, Eigen::Affine>, Eigen::Transform<Scalar, 3, Eigen::Affine>>> _mats;
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