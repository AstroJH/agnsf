#pragma once

#include <cstddef>
#include <limits>
#include <vector>

#include <esf/sf_uncertainty.hpp>

namespace agnsf {
namespace esf {

struct SFBinResult {
    // Number of pairs contributing to this bin.
    std::size_t count = 0;

    // Noise-corrected second-order structure function.
    double sf_squared =
        std::numeric_limits<double>::quiet_NaN();

    // Structure function.
    double sf =
        std::numeric_limits<double>::quiet_NaN();

    // Measurement uncertainty on SF or ESF (NaN when not estimated).
    SFUncertainty measurement;

    // Source-to-source sampling uncertainty on ESF.
    //
    // Single light curves never estimate sampling uncertainty; this
    // stays NaN by design.
    SFUncertainty sampling;
};


class SFResult {
public:
    SFResult() = default;

    explicit SFResult(
        std::vector<SFBinResult> bins
    );

    std::size_t size() const noexcept;

    const SFBinResult& bin(
        std::size_t index
    ) const;

    const std::vector<SFBinResult>& bins() const noexcept;

private:
    std::vector<SFBinResult> bins_;
};

} // namespace esf
} // namespace agnsf