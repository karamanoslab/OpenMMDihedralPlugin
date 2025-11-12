#ifndef OPENMM_REFERENCEDIHEDRALKERNELFACTORY_H_
#define OPENMM_REFERENCEDIHEDRALKERNELFACTORY_H_



#include "openmm/KernelFactory.h"

namespace OpenMM {

/**
 * This KernelFactory creates kernels for the reference implementation of the Example plugin.
 */

class ReferenceDihedralKernelFactory : public KernelFactory {
public:
    KernelImpl* createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const;
};

} // namespace OpenMM

#endif /*OPENMM_REFERENCEDIHEDRALKERNELFACTORY_H_*/
