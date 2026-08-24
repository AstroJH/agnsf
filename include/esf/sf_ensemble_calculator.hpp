#pragma once

#include <vector>

#include <esf/lag_bins.hpp>
#include <esf/light_curve.hpp>
#include <esf/sf_result.hpp>

namespace esf {

/**
 * Ensemble SF obtained by first computing an SF for each
 * light curve and then averaging the individual SF values.
 *
 * For each lag bin, only light curves with a finite SF
 * contribute to the mean.
 */
class SFEnsembleCalculator {
public:
    SFResult calculate(
        const std::vector<LightCurve>& data,
        const LagBins& bins
    ) const;
};

} // namespace esf