#pragma once

#include <cstddef>
#include <limits>
#include <vector>

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