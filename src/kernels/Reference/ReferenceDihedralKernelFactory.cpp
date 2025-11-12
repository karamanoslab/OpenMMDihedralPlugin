
#include "ReferenceDihedralKernelFactory.h"
#include "ReferenceDihedralKernels.h"
#include "openmm/reference/ReferencePlatform.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace DihedralPlugin;
using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    for (int i = 0; i < Platform::getNumPlatforms(); i++) {
        Platform& platform = Platform::getPlatform(i);
        if (dynamic_cast<ReferencePlatform*>(&platform) != NULL) {
            ReferenceDihedralKernelFactory* factory = new ReferenceDihedralKernelFactory();
            platform.registerKernelFactory(CalcDihedralForceKernel::Name(), factory);
        }
    }
}

extern "C" OPENMM_EXPORT void registerDihedralReferenceKernelFactories() {
    registerKernelFactories();
}

KernelImpl* ReferenceDihedralKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    ReferencePlatform::PlatformData& data = *static_cast<ReferencePlatform::PlatformData*>(context.getPlatformData());
    if (name == CalcDihedralForceKernel::Name())
        return new ReferenceCalcDihedralForceKernel(name, platform);
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}
