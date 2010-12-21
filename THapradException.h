#include <stdexcept>

using std::runtime_error;


class TKinematicException : public runtime_error {
public:
    TKinematicException()
        : runtime_error("Wrong kinematics! Skip the point!") {};
};
