#include <esf/sf_calculator.hpp>

#include <esf/bin_accumulator.hpp>
#include <esf/sf_uncertainty_estimator.hpp>

namespace agnsf {
namespace esf {

namespace {

/**
 * Validate the uncertainty configuration for a single SF.
 *
 * Measurement supports Analytic (Monte Carlo may be added later);
 * sampling is an ESF concept and must be Off here.
 */
void validate_config(
    const UncertaintyConfig& config
)
{
    if (config.measurement != UncertaintyMethod::Off &&
        config.measurement != UncertaintyMethod::Analytic) {

        throw std::invalid_argument(
            "SFCalculator: measurement uncertainty supports "
            "only Off or Analytic"
        );
    }

    if (config.sampling != UncertaintyMethod::Off) {

        throw std::invalid_argument(
            "SFCalculator: sampling uncertainty is not defined "
            "for a single light curve"
        );
    }
}


SFResult calculate_impl(
    const double* time,
    const double* value,
    const double* error,
    std::size_t size,
    const LagBins& bins,
    SFMethod method,
    const UncertaintyConfig& config,
    double redshift
)
{
    validate_config(config);

    if (redshift <= -1.0) {
        throw std::invalid_argument(
            "redshift must be > -1"
        );
    }

    // Rest-frame lag scale
    const double lag_scale = 1.0 / (1.0 + redshift);

    std::vector<BinAccumulator> accumulators(
        bins.size()
    );

    const double min_lag = bins.min();
    const double max_lag = bins.max();

    for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = i + 1; j < size; ++j) {
            const double dt = time[j] - time[i];
            const double lag = dt * lag_scale; // to rest frame

            if (lag >= max_lag) {
                break;
            }

            if (lag < min_lag) {
                continue;
            }

            std::size_t bin;

            if (!bins.try_index(lag, bin)) {
                continue;
            }

            const double delta =
                value[j] - value[i];

            accumulators[bin].add(
                delta,
                error[i],
                error[j]
            );
        }
    }

    // Measurement-uncertainty estimator, created once per call.
    std::unique_ptr<SFMeasurementUncertaintyEstimator>
        measurement_estimator;

    if (config.measurement == UncertaintyMethod::Analytic) {
        measurement_estimator =
            make_sf_measurement_uncertainty_estimator(method);
    }

    std::vector<SFBinResult> results;
    results.reserve(bins.size());

    for (const auto& accumulator : accumulators) {

        SFBinResult result;

        result.count = accumulator.count();
        result.sf_squared = accumulator.sf_squared(method);
        result.sf = accumulator.sf(method);

        if (measurement_estimator) {
            result.measurement =
                measurement_estimator->estimate(accumulator);
        }

        results.push_back(result);
    }

    return SFResult(
        std::move(results)
    );
}

} // namespace


SFResult SFCalculator::calculate(
    const agnsf::LightCurve& data,
    const LagBins& bins,
    SFMethod method,
    const UncertaintyConfig& config,
    double redshift
) const
{
    return calculate_impl(
        data.time_data(),
        data.value_data(),
        data.error_data(),
        data.size(),
        bins,
        method,
        config,
        redshift
    );
}


SFResult SFCalculator::calculate(
    const agnsf::LightCurveView& data,
    const LagBins& bins,
    SFMethod method,
    const UncertaintyConfig& config,
    double redshift
) const
{
    return calculate_impl(
        data.time_data(),
        data.value_data(),
        data.error_data(),
        data.size(),
        bins,
        method,
        config,
        redshift
    );
}

} // namespace esf
} // namespace agnsf
