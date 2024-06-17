// benchPreds.hpp

#ifndef GKMC_BenchPreds_H
#define GKMC_BenchPreds_H

namespace gkmc {

template <class Sys, class Logger> struct StepsPredicate {
  StepsPredicate(size_t maxSteps, Logger logger)
    : _maxSteps{maxSteps}, _logger{logger} {}
  auto operator()(Sys &sys, double t) {
    _logger(sys, t, ++_steps);
    return _steps < _maxSteps;
  }
  size_t steps() const { return _steps; }
private:
  size_t _maxSteps;
  Logger _logger;
  size_t _steps = 0;
};

template <class Sys, class Logger> struct TimePredicate {
  TimePredicate(double maxTime, Logger logger)
    : _maxTime{maxTime}, _logger{logger} {}
  auto operator()(Sys &sys, double t) {
    _logger(sys, t, _steps++);
    return t < _maxTime;
  }
  size_t steps() const { return _steps; }
private:
  double _maxTime;
  Logger _logger;
  size_t _steps = 0;
};

} // namespace gkmc

#endif