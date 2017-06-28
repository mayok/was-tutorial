#include <complex>
#include <vector>
#include <math.h>

class MFCC {
private:
  const double PI = 4 * atan(1.0);

private:
  inline double hz2mel(double f) {
    return 2595 * std::log10(1 + f/700);
  }

  inline double mel2hz(double m) {
    return 700 * (std::pow(10, m/2595) - 1);
  }

public:
  // constructor
  MFCC() {

  }

  // main process
  // TODO: should be return some value
  void process() {

  }
}
