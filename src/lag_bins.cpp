#include <esf/lag_bins.hpp>

#include <algorithm>
#include <cmath>
#include <iterator>
#include <stdexcept>
#include <utility>

namespace esf {

LagBins::LagBins(
    std::vector<double> edges
)
    : grid_type_(GridType::Custom),
      edges_(std::move(edges))
{
    if (edges_.size() < 2) {
        throw std::invalid_argument(
            "LagBins requires at least two edges."
        );
    }

    for (double edge : edges_) {
        if (!std::isfinite(edge)) {
            throw std::invalid_argument(
                "Bin edges must be finite."
            );
        }

        if (edge < 0.0) {
            throw std::invalid_argument(
                "Bin edges must be non-negative."
            );
        }
    }

    for (std::size_t i = 1; i < edges_.size(); ++i) {
        if (edges_[i] <= edges_[i - 1]) {
            throw std::invalid_argument(
                "Bin edges must be strictly increasing."
            );
        }
    }
}


LagBins::LagBins(
    GridType grid_type,
    std::vector<double> edges
)
    : grid_type_(grid_type),
      edges_(std::move(edges))
{
}


LagBins LagBins::linear(
    double min,
    double max,
    double step
)
{
    if (!std::isfinite(min) ||
        !std::isfinite(max) ||
        !std::isfinite(step)) {
        throw std::invalid_argument(
            "Linear lag-bin parameters must be finite."
        );
    }

    if (min < 0.0) {
        throw std::invalid_argument(
            "Linear lag-bin minimum must be non-negative."
        );
    }

    if (max <= min) {
        throw std::invalid_argument(
            "Linear lag-bin maximum must be greater than minimum."
        );
    }

    if (step <= 0.0) {
        throw std::invalid_argument(
            "Linear lag-bin step must be positive."
        );
    }

    std::vector<double> edges;

    edges.push_back(min);

    double edge = min;

    while (edge + step < max) {
        edge += step;
        edges.push_back(edge);
    }

    if (edges.back() != max) {
        edges.push_back(max);
    }

    return LagBins(
        GridType::Linear,
        std::move(edges)
    );
}


LagBins LagBins::logarithmic(
    double min,
    double max,
    double step
)
{
    if (!std::isfinite(min) ||
        !std::isfinite(max) ||
        !std::isfinite(step)) {
        throw std::invalid_argument(
            "Logarithmic lag-bin parameters must be finite."
        );
    }

    if (min <= 0.0) {
        throw std::invalid_argument(
            "Logarithmic lag bins require a positive minimum."
        );
    }

    if (max <= min) {
        throw std::invalid_argument(
            "Logarithmic lag-bin maximum must be greater than minimum."
        );
    }

    if (step <= 0.0) {
        throw std::invalid_argument(
            "Logarithmic lag-bin step must be positive."
        );
    }

    const double log_min = std::log10(min);
    const double log_max = std::log10(max);

    std::vector<double> edges;

    edges.push_back(min);

    double log_edge = log_min + step;

    while (log_edge < log_max) {
        edges.push_back(
            std::pow(10.0, log_edge)
        );

        log_edge += step;
    }

    edges.push_back(max);

    return LagBins(
        GridType::Logarithmic,
        std::move(edges)
    );
}


LagBins::GridType
LagBins::grid_type() const noexcept
{
    return grid_type_;
}


double LagBins::min() const noexcept
{
    return edges_.front();
}


double LagBins::max() const noexcept
{
    return edges_.back();
}


std::size_t LagBins::size() const noexcept
{
    return edges_.size() - 1;
}


const std::vector<double>&
LagBins::edges() const noexcept
{
    return edges_;
}


bool LagBins::contains(double lag) const noexcept
{
    std::size_t index_value;

    return try_index(lag, index_value);
}


bool LagBins::try_index(
    double lag,
    std::size_t& index_value
) const noexcept
{
    if (!std::isfinite(lag)) {
        return false;
    }

    if (lag < min() || lag >= max()) {
        return false;
    }

    auto it = std::upper_bound(
        edges_.begin(),
        edges_.end(),
        lag
    );

    index_value =
        static_cast<std::size_t>(
            std::distance(
                edges_.begin(),
                it
            )
        ) - 1;

    return true;
}


std::size_t LagBins::index(double lag) const
{
    std::size_t index_value;

    if (!try_index(lag, index_value)) {
        throw std::out_of_range(
            "Lag is outside the bin range."
        );
    }

    return index_value;
}

} // namespace esf