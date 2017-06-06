#include <Eigen/Dense>

#include "python_shared.h"

% for f in functions:
#include <IGsim/${f}.h>
% endfor

void python_export_sim(py::module &m)
{
% for f in functions:
  #include "py_sim/py_${f}.cpp"
% endfor
}
