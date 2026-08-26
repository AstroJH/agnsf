#include <esf/sf_result.hpp>

#include <stdexcept>
#include <utility>

namespace agnsf {
namespace esf {

SFResult::SFResult(
    std::vector<SFBinResult> bins
)
    : bins_(std::move(bins))
{
}


std::size_t SFResult::size() const noexcept
{
    return bins_.size();
}


const SFBinResult& SFResult::bin(
    std::size_t index
) const
{
    if (index >= bins_.size()) {
        throw std::out_of_range(
            "SF bin index is out of range."
        );
    }

    return bins_[index];
}


const std::vector<SFBinResult>&
SFResult::bins() const noexcept
{
    return bins_;
}

} // namespace esf
} // namespace agnsf