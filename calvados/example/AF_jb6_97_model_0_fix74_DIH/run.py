from calvados import sim
from argparse import ArgumentParser
import openmm



if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--path',nargs='?', default='.', const='.', type=str)
    parser.add_argument('--config',nargs='?', default='config.yaml', const='config.yaml', type=str)
    parser.add_argument('--components',nargs='?', default='components.yaml', const='components.yaml', type=str)

    args = parser.parse_args()

    path = args.path
    fconfig = args.config
    fcomponents = args.components

    sim.run(path=path,fconfig=fconfig,fcomponents=fcomponents)


from calvados.analysis import save_rg

save_rg("/localdata/tkaraman/Programs/CALVADOS/JDGFtest/AF_jb6_97_model_0_fix74_DIH","AF_jb6_97_model_0","/localdata/tkaraman/Programs/CALVADOS/JDGFtest/input/residues.csv",".",10)
