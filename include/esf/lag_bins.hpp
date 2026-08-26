#pragma once

#include <cstddef>
#include <vector>

namespace agnsf {
namespace esf {

class LagBins {
public:
    enum class GridType {
        Custom,
        Linear,
        Logarithmic
    };

    explicit LagBins(
        std::vector<double> edges
    );

    static LagBins linear(
        double min,
        double max,
        double step
    );

    static LagBins logarithmic(
        double min,
        double max,
        double step
    );

    GridType grid_type() const noexcept;

    double min() const noexcept;
    double max() const noexcept;

    std::size_t size() const noexcept;

    const std::vector<double>& edges() const noexcept;

    bool contains(double lag) const noexcept;

    bool try_index(
        double lag,
        std::size_t& index
    ) const noexcept;

    std::size_t index(double lag) const;

private:
    LagBins(
        GridType grid_type,
        std::vector<double> edges
    );

    GridType grid_type_;

    std::vector<double> edges_;
};

} // namespace esf
} // namespace agnsf