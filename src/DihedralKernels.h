#ifndef DIHEDRAL_KERNELS_H_
#define DIHEDRAL_KERNELS_H_


#include "DihedralForce.h"
#include "openmm/KernelImpl.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include <string>

namespace DihedralPlugin {

/**
 * This kernel is invoked by DihedralForce to calculate the forces acting on the system and the energy of the system.
 */
class CalcDihedralForceKernel : public OpenMM::KernelImpl {
public:
    static std::string Name() {
        return "CalcDihedralForce";
    }
    CalcDihedralForceKernel(std::string name, const OpenMM::Platform& platform) : OpenMM::KernelImpl(name, platform) {
    }
    /**
     * Initialize the kernel.
     */
    virtual void initialize(const OpenMM::System& system, const DihedralForce& force) = 0;
    /**
    
     */
    virtual double execute(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     */
    virtual void copyParametersToContext(OpenMM::ContextImpl& context, const DihedralForce& force) = 0;
};

} // namespace DihedralPlugin

#endif /*DIHEDRAL_KERNELS_H_*/
