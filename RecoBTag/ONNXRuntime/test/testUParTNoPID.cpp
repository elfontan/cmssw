#include <numeric>
#include <stdexcept>
#include <vector>

#include "RecoBTag/ONNXRuntime/interface/UParTNoPID.h"

int main() {
  // Distinct rows catch both feature-order mistakes and in-place row corruption.
  std::vector<float> values(75, 0.f);
  std::iota(values.begin(), values.begin() + 25, 0.f);
  std::iota(values.begin() + 25, values.begin() + 50, 100.f);
  const std::vector<float> expected{
      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 22, 23, 24,
      100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 122, 123, 124,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  UparTNoPID::removeChargedPIDFeatures(values);
  if (values != expected)
    throw std::runtime_error("NoPID projection changed feature order, hit counts, row boundaries or padding");
  values.clear();
  UparTNoPID::removeChargedPIDFeatures(values);
  if (!values.empty())
    throw std::runtime_error("NoPID projection changed an empty input");
}
