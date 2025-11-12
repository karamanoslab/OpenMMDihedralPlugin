#ifndef DIHEDRAL_KERNELS_H_
#define DIHEDRAL_KERNELS_H_


#include "DihedralForce.h"
#include "openmm/KernelImpl.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include <string>

namespace DihedralPlugin {

/**
 * This kernel is invoked by ExampleForce to calculate the forces acting on the system and the energy of the system.
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
     * 
     * @param system     the System this kernel will be applied to
     * @param force      the ExampleForce this kernel will be used for
     */
    virtual void initialize(const OpenMM::System& system, const DihedralForce& force) = 0;
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    virtual double execute(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy) = 0;
    /**
     * Copy changed parameters over to a context.
     *
     * @param context    the context to copy parameters to
     * @param force      the ExampleForce to copy the parameters from
     */
    virtual void copyParametersToContext(OpenMM::ContextImpl& context, const DihedralForce& force) = 0;
};

} // namespace DihedralPlugin

#endif /*DIHEDRAL_KERNELS_H_*/
