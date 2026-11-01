#include "grating.h"
#include "crossgrating2.h"
#include "gcross.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace SCATMECH;

struct Input {
    std::string model;
    double wavelength;
    int resolution;
    std::vector<std::pair<std::string, std::string> > params;
};

void number(std::ostream& out, double value) {
    if (!std::isfinite(value)) throw std::runtime_error("Geometry contains a nonfinite value");
    out << value;
}

void component(std::ostream& out, COMPLEX value) {
    out << '[';
    number(out, real(value));
    out << ',';
    number(out, imag(value));
    out << ']';
}

void isotropic(std::ostream& out, const dielectric_function& material, double wavelength) {
    COMPLEX value = material.epsilon(wavelength);
    out << '[';
    for (int c = 0; c < 6; ++c) {
        if (c) out << ',';
        component(out, c < 3 ? value : COMPLEX(1, 0));
    }
    out << ']';
}

void configure(Model& model, const Input& input) {
    model.init();
    for (size_t i = 0; i < input.params.size(); ++i)
        model.set_parameter(input.params[i].first, input.params[i].second);
}

void edges(std::ostream& out, const std::vector<double>& values) {
    out << '[';
    for (size_t i = 0; i < values.size(); ++i) {
        if (i) out << ',';
        number(out, values[i]);
    }
    out << ']';
}

std::vector<int> display_bounds(int count, int resolution) {
    int size = std::min(count, resolution);
    std::vector<int> result(size + 1);
    for (int i = 0; i <= size; ++i) result[i] = i * count / size;
    return result;
}

std::vector<double> fractions(const std::vector<int>& bounds, int count) {
    std::vector<double> result;
    for (size_t i = 0; i < bounds.size(); ++i) result.push_back(double(bounds[i]) / count - 0.5);
    return result;
}

void export_one(Grating& model, const Input& input, std::ostream& out) {
    configure(model, input);
    model.set_lambda(input.wavelength);
    int levels = model.get_levels();
    double period = model.get_period();
    bool continuous = input.model == "Sinusoidal_Volume_Grating";
    out << "{\"model\":\"" << input.model << "\",\"dimension\":1,\"lattice\":[[";
    number(out, period);
    out << ",0],[0,";
    number(out, period);
    out << "]],\"medium_i\":";
    isotropic(out, model.get_medium_i(), input.wavelength);
    out << ",\"medium_t\":";
    isotropic(out, model.get_medium_t(), input.wavelength);
    out << ",\"layers\":[";
    for (int level = 0; level < levels; ++level) {
        if (level) out << ',';
        std::vector<double> x;
        x.push_back(-0.5);
        if (continuous) {
            for (int i = 1; i < input.resolution; ++i) x.push_back(double(i) / input.resolution - 0.5);
        } else {
            const std::vector<double>& native = model.get_position()[level];
            for (size_t i = 0; i < native.size(); ++i) {
                double u = native[i] / period;
                u -= std::floor(u + 0.5);
                if (u > -0.5 + 1e-14 && u < 0.5 - 1e-14) x.push_back(u);
            }
        }
        x.push_back(0.5);
        std::sort(x.begin(), x.end());
        x.erase(std::unique(x.begin(), x.end(), [](double a, double b) { return std::abs(a - b) < 1e-12; }), x.end());
        out << "{\"thickness\":";
        number(out, model.get_thickness(level));
        out << ",\"u_edges\":";
        edges(out, x);
        out << ",\"v_edges\":[-0.5,0.5],\"continuous\":" << (continuous ? "true" : "false") << ",\"tensors\":[[";
        for (size_t i = 0; i + 1 < x.size(); ++i) {
            if (i) out << ',';
            double midpoint = (x[i] + x[i + 1]) * 0.5 * period;
            out << '[';
            for (int c = 0; c < 6; ++c) {
                if (c) out << ',';
                component(out, c < 3 ? model.eps(midpoint, level, c) : model.mu(midpoint, level, c - 3));
            }
            out << ']';
        }
        out << "]]}";
    }
    out << "],\"notes\":[\"";
    if (continuous) out << "Continuous native permittivity sampled at " << input.resolution << " display cells per period";
    else out << "Exact native boundary segments retained";
    out << "; all " << levels << " native layers retained\"]}";
}

template<class T> class Inspectable : public T {
public:
    COMPLEX value(int x, int y, int layer, int c) {
        if (c < 3) {
            if (this->isotropic) return this->eps(x, y, layer);
            return (c == 0 ? this->eps1 : c == 1 ? this->eps2 : this->eps3)(x, y, layer);
        }
        if (this->nonmagnetic) return COMPLEX(1, 0);
        if (this->isotropic) return this->mu(x, y, layer);
        return (c == 3 ? this->mu1 : c == 4 ? this->mu2 : this->mu3)(x, y, layer);
    }
};

template<class T> void radial_contour(Inspectable<T>&, const Input&, int, int, std::ostream&) {}

void radial_contour(Inspectable<Cylinder_CrossGrating>& model, const Input& input, int level, int levels, std::ostream& out) {
    double pi = std::acos(-1.);
    double zeta = model.get_zeta() * pi / 180.;
    int count = std::max(128, std::min(512, input.resolution * 4));
    std::vector<double> angles;
    for (int index = 0; index <= count; ++index) angles.push_back(index == count ? pi : -pi + 2. * pi * index / count);
    for (int u = -1; u <= 1; u += 2) {
        for (int v = -1; v <= 1; v += 2) {
            double x = u * model.get_d1() + v * model.get_d2() * std::sin(zeta);
            double y = v * model.get_d2() * std::cos(zeta);
            double angle = std::atan2(y, x);
            angles.push_back(angle >= pi ? angle - 2. * pi : angle);
        }
    }
    std::sort(angles.begin(), angles.end());
    angles.erase(std::unique(angles.begin(), angles.end(), [](double a, double b) { return std::abs(a - b) < 1e-12; }), angles.end());
    std::vector<double> radii;
    for (size_t index = 0; index < angles.size(); ++index) {
        double degrees = angles[index] * 180. / pi;
        double bottom = model.get_rbottom().value(degrees);
        double top = model.get_rtop().value(degrees);
        radii.push_back(bottom + (level - 0.5) / levels * (top - bottom));
    }
    out << ",\"radial\":{\"angles\":";
    edges(out, angles);
    out << ",\"radii\":";
    edges(out, radii);
    out << ",\"inside\":";
    isotropic(out, model.get_inside(), input.wavelength);
    out << ",\"outside\":";
    isotropic(out, model.get_outside(), input.wavelength);
    out << '}';
}

void cross_header(CrossGrating& model, const Input& input, std::ostream& out) {
    double angle = model.get_zeta() * std::acos(-1.) / 180.;
    out << "{\"model\":\"" << input.model << "\",\"dimension\":2,\"lattice\":[[";
    number(out, model.get_d1());
    out << ",0],[";
    number(out, model.get_d2() * std::sin(angle));
    out << ',';
    number(out, model.get_d2() * std::cos(angle));
    out << "]],\"medium_i\":";
    isotropic(out, model.get_medium_i(), input.wavelength);
    out << ",\"medium_t\":";
    isotropic(out, model.get_medium_t(), input.wavelength);
    out << ",\"layers\":[";
}

template<class T> void export_cross(const Input& input, std::ostream& out) {
    Inspectable<T> model;
    configure(model, input);
    model.set_lambda(input.wavelength);
    model.set_order1(0);
    model.set_order2(0);
    int levels = model.get_levels();
    int nx = model.get_grid1();
    int ny = model.get_grid2();
    std::vector<int> bx = display_bounds(nx, input.resolution);
    std::vector<int> by = display_bounds(ny, input.resolution);
    cross_header(model, input, out);
    for (int level = levels; level >= 1; --level) {
        if (level != levels) out << ',';
        out << "{\"thickness\":";
        number(out, model.get_thick(level));
        out << ",\"u_edges\":";
        edges(out, fractions(bx, nx));
        out << ",\"v_edges\":";
        edges(out, fractions(by, ny));
        out << ",\"continuous\":false,\"tensors\":[";
        for (size_t j = 0; j + 1 < by.size(); ++j) {
            if (j) out << ',';
            int y = (by[j] + by[j + 1]) / 2 + 1;
            out << '[';
            for (size_t i = 0; i + 1 < bx.size(); ++i) {
                if (i) out << ',';
                int x = (bx[i] + bx[i + 1]) / 2 + 1;
                out << '[';
                for (int c = 0; c < 6; ++c) {
                    if (c) out << ',';
                    component(out, model.value(x, y, level, c));
                }
                out << ']';
            }
            out << ']';
        }
        out << ']';
        radial_contour(model, input, level, levels, out);
        out << '}';
    }
    out << "],\"notes\":[\"Native grid " << nx << " x " << ny << "; display grid "
        << bx.size() - 1 << " x " << by.size() - 1 << "; all " << levels << " native layers retained\"";
    if (input.model == "Cylinder_CrossGrating") out << ",\"Smooth cylinder contours evaluated from native radius formulas and tables at each solver layer; solver grid unchanged\"";
    out << "]}";
}

int main(int argc, char** argv) {
    try {
        if (argc != 3) throw std::runtime_error("Expected request and output filenames");
        std::ifstream request(argv[1]);
        Input input;
        std::string line;
        if (!std::getline(request, input.model)) throw std::runtime_error("Missing model");
        if (!std::getline(request, line)) throw std::runtime_error("Missing wavelength");
        input.wavelength = std::stod(line);
        if (!std::getline(request, line)) throw std::runtime_error("Missing resolution");
        input.resolution = std::stoi(line);
        if (!(input.wavelength > 0) || input.resolution < 1) throw std::runtime_error("Invalid preview parameters");
        while (std::getline(request, line)) {
            size_t tab = line.find('\t');
            if (tab == std::string::npos) throw std::runtime_error("Invalid named parameter");
            input.params.push_back(std::make_pair(line.substr(0, tab), line.substr(tab + 1)));
        }
        std::ofstream out(argv[2]);
        if (!out) throw std::runtime_error("Cannot open geometry output");
        out << std::setprecision(17);
        std::unique_ptr<Grating> one;
        if (input.model == "Single_Line_Grating") one.reset(new Single_Line_Grating);
        else if (input.model == "Corner_Rounded_Grating") one.reset(new Corner_Rounded_Grating);
        else if (input.model == "Sinusoidal_Relief_Grating") one.reset(new Sinusoidal_Relief_Grating);
        else if (input.model == "Sinusoidal_Volume_Grating") one.reset(new Sinusoidal_Volume_Grating);
        else if (input.model == "Triangular_Grating") one.reset(new Triangular_Grating);
        else if (input.model == "Generic_Grating") one.reset(new Generic_Grating);
        if (one) export_one(*one, input, out);
        else if (input.model == "Cylinder_CrossGrating") export_cross<Cylinder_CrossGrating>(input, out);
        else if (input.model == "Rectangle_CrossGrating") export_cross<Rectangle_CrossGrating>(input, out);
        else if (input.model == "Sphere_CrossGrating") export_cross<Sphere_CrossGrating>(input, out);
        else if (input.model == "Pyramidal_Pit_CrossGrating") export_cross<Pyramidal_Pit_CrossGrating>(input, out);
        else if (input.model == "Generic_CrossGrating") export_cross<Generic_CrossGrating>(input, out);
        else if (input.model == "Null_CrossGrating") {
            Null_CrossGrating model;
            configure(model, input);
            model.set_lambda(input.wavelength);
            model.set_order1(0);
            model.set_order2(0);
            model.get_levels();
            cross_header(model, input, out);
            out << "],\"notes\":[\"Flat native interface\"]}";
        } else throw std::runtime_error("Unsupported native leaf model: " + input.model);
        out.close();
        if (!out) throw std::runtime_error("Cannot write geometry output");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << std::endl;
        return 1;
    }
}
