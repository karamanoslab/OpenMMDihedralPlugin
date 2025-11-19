#ifndef REFERENCE_DIHEDRAL_KERNELS_H_
#define REFERENCE_DIHEDRAL_KERNELS_H_


#include "DihedralKernels.h"
#include "openmm/Platform.h"
#include <vector>

namespace DihedralPlugin {

/**
 * This kernel is invoked by DihedralForce to calculate the forces acting on the system and the energy of the system.
 */
class ReferenceCalcDihedralForceKernel : public CalcDihedralForceKernel {
public:
    ReferenceCalcDihedralForceKernel(std::string name, const OpenMM::Platform& platform) : CalcDihedralForceKernel(name, platform) {
    }
    /**
     * Initialize the kernel.
     */
    void initialize(const OpenMM::System& system, const DihedralForce& force);
    /**
     * Execute the kernel to calculate the forces and/or energy.
     *
     * @param context        the context in which to execute this kernel
     * @param includeForces  true if forces should be calculated
     * @param includeEnergy  true if the energy should be calculated
     * @return the potential energy due to the force
     */
    double execute(OpenMM::ContextImpl& context, bool includeForces, bool includeEnergy);
    /**
     * Copy changed parameters over to a context.
     */
    void copyParametersToContext(OpenMM::ContextImpl& context, const DihedralForce& force);
private:
    int numDihedrals;
    std::vector<int> particle1, particle2, particle3, particle4;
    std::vector<double> ka1, phia1,  epsd,  ka2,  phia2,  eps0,  kb1,  phib1, eps1,  kb2,  phib2,  eps2;
    bool periodic;
};

} // namespace DihedralPlugin

#endif /*REFERENCE_DIHEDRAL_KERNELS_H_*/
