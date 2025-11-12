#ifndef OPENMM_DIHEDRALFORCEIMPL_H_
#define OPENMM_DIHEDRALFORCEIMPL_H_


#include "DihedralForce.h"
#include "openmm/internal/ForceImpl.h"
#include "openmm/Kernel.h"
#include <utility>
#include <set>
#include <string>

namespace DihedralPlugin {

class System;

/**
 * This is the internal implementation of ExampleForce.
 */

class OPENMM_EXPORT_EXAMPLE DihedralForceImpl : public OpenMM::ForceImpl {
public:
    DihedralForceImpl(const DihedralForce& owner);
    ~DihedralForceImpl();
    void initialize(OpenMM::ContextImpl& context);
    const DihedralForce& getOwner() const {
        return owner;
    }
    void updateContextState(OpenMM::ContextImpl& context, bool& forcesInvalid) {
        // This force field doesn't update the state directly.
    }
    double calcForcesAndEnergy(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy, int groups);
    std::map<std::string, double> getDefaultParameters() {
        return std::map<std::string, double>(); // This force field doesn't define any parameters.
    }
    std::vector<std::string> getKernelNames();
    std::vector<std::tuple<int, int, int, int> > getDihedralParticles() const;
    void updateParametersInContext(OpenMM::ContextImpl& context);
private:
    const DihedralForce& owner;
    OpenMM::Kernel kernel;
};

} // namespace DihedralPlugin

#endif /*OPENMM_DIHEDRALFORCEIMPL_H_*/
