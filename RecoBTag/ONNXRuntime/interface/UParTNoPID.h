#ifndef RecoBTag_ONNXRuntime_UParTNoPID_h
#define RecoBTag_ONNXRuntime_UParTNoPID_h

#include <cassert>
#include <cstddef>
#include <vector>

#include "RecoBTag/ONNXRuntime/interface/tensor_configs.h"

namespace UparTNoPID {
  inline constexpr unsigned n_charged_features = 22;

  // Project the release's charged-candidate layout onto the custom NoPID model:
  // omit HadFrac, CaloFrac and pdgID (columns 19, 20, 21), retaining hit counts.
  inline void removeChargedPIDFeatures(std::vector<float>& values) {
    constexpr auto width = UparT::N_InputFeatures.at(UparT::kChargedCandidates);
    static_assert(width == n_charged_features + 3, "Review the NoPID feature mapping after layout changes");
    assert(values.size() % width == 0);
    const auto candidates = values.size() / width;
    std::size_t output = 0;
    for (std::size_t candidate = 0; candidate < candidates; ++candidate) {
      for (unsigned feature = 0; feature < width; ++feature) {
        if (feature < 19 || feature > 21)
          values[output++] = values[candidate * width + feature];
      }
    }
    values.resize(output);
  }
}  // namespace UparTNoPID

#endif
