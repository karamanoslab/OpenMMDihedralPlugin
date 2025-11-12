from openmm.app import *
from openmm import *
from openmm.unit import *
import numpy as np

Platform.loadPluginsFromDirectory("/localdata/tkaraman/Projects/openmm_tests/torsion_plugin/build")

import dihedralplugin

def compute_energy(context):
    state = context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)


def run_test_dihedral_force():
    epsilon = 1e-4  # nm
    tolerance = 1e-3  # kJ/mol/nm
    

    # Create system
    system = System()
    for _ in range(4):
        system.addParticle(12.0)

    #set box for periodic boundary conditions
    x=Quantity(np.array([20, 0,  0]), nanometers)
    y=Quantity(np.array([0,  20, 0]), nanometers)
    z=Quantity(np.array([0,  0, 20]), nanometers)
    system.setDefaultPeriodicBoxVectors(x, y, z)

    # Create and add DihedralForce
    force = dihedralplugin.DihedralForce()
    force.addDihedral(0, 1, 2, 3,
                   11.4, 0.9, -2, 0.15, 1.02, #double ka1, double phia1, double epsd, double ka2, double phia2
                   0.27, 1.8, -1.55, 0.14,   #double eps0, double kb1, double phib1, double eps1,
                   0.65, -2.5, 0.26 )        #double kb2, double phib2, double eps2

    system.addForce(force)
    force.setUsesPeriodicBoundaryConditions(True)
    print('Periodic conditions ', force.usesPeriodicBoundaryConditions())


    integrator = VerletIntegrator(0.001)
    platform = Platform.getPlatformByName("Reference")
    context = Context(system, integrator, platform)

    # Initial positions
    positions = [
        Vec3(1, 0, 0),
        Vec3(0, 0, 0),
        Vec3(0, 1, 0),
        Vec3(1, 1, 1)
    ] * nanometers
    context.setPositions(positions)

    # Analytical forces
    state = context.getState(getEnergy=True, getForces=True)
    analytical_forces = state.getForces(asNumpy=True)
    energy = state.getPotentialEnergy()
    print("Analytical Energy =", energy)

    print("Analytical Forces:")
    for i in range(4):
        print(f"{i}: {analytical_forces[i]}")

    print("\nNumerical Forces:")
    for i in range(4):
        for d in range(3):
            displaced_pos = positions.value_in_unit(nanometers).copy()
            vec = displaced_pos[i]
            if d == 0:
                vec = Vec3(vec[0] + epsilon, vec[1], vec[2])
            elif d == 1:
                vec = Vec3(vec[0], vec[1] + epsilon, vec[2])
            else:
                vec = Vec3(vec[0], vec[1], vec[2] + epsilon)
            displaced_pos[i] = vec
            context.setPositions(displaced_pos * nanometers)
            e_plus = compute_energy(context)

            # Backward displacement
            displaced_pos = positions.value_in_unit(nanometers).copy()
            vec = displaced_pos[i]
            if d == 0:
                vec = Vec3(vec[0] - epsilon, vec[1], vec[2])
            elif d == 1:
                vec = Vec3(vec[0], vec[1] - epsilon, vec[2])
            else:
                vec = Vec3(vec[0], vec[1], vec[2] - epsilon)
            displaced_pos[i] = vec
            context.setPositions(displaced_pos * nanometers)
            e_minus = compute_energy(context)


            f_num = -(e_plus - e_minus) / (2 * epsilon)
            f_ana = analytical_forces[i][d].value_in_unit(kilojoule_per_mole/nanometer)
            diff = abs(f_num - f_ana)
            print(f"Atom {i}, dim {d}: f_num = {f_num:.5f}, f_ana = {f_ana:.5f}, diff = {diff:.5g}")

            if diff > tolerance:
                raise RuntimeError("Force mismatch: check derivatives.")

    print("\nTest passed: forces match within tolerance.")

if __name__ == "__main__":
    run_test_dihedral_force()