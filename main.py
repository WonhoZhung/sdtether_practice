import os


AVAIL_IN_FORMAT_LIST = ["sdf", "pdb", "mol2", "xyz"]


def generate_prm(
        ref_poc_fn,
        ref_lig_fn,
        out_prm_fn,
        title="",
        max_trans=1.0,
        max_rot=10.0,
        max_dihedral=10.0
        ):
    assert ref_poc_fn.split('.')[-1] == "mol2", "Wrong receptor format" 
    assert ref_lig_fn.split('.')[-1] == "sdf", "Wrong ligand format"
    assert out_prm_fn.split('.')[-1] == "prm", "Wrong param format"
    
    lines = f"""RBT_PARAMETER_FILE_V1.00
TITLE {title}

RECEPTOR_FILE {ref_poc_fn}

SECTION LIGAND
    TRANS_MODE TETHERED
    ROT_MODE TETHERED 
    DIHEDRAL_MODE TETHERED 
    MAX_TRANS {max_trans} 
    MAX_ROT {max_rot}
    MAX_DIHEDRAL {max_dihedral}
END_SECTION

##################################################################
### CAVITY DEFINITION: REFERENCE LIGAND METHOD
##################################################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL {ref_lig_fn}
    RADIUS 6.0
    SMALL_SPHERE 1.0
    MIN_VOLUME 100
    MAX_CAVITIES 1
    VOL_INCR 0.0
   GRIDSTEP 0.5
END_SECTION

#################################
#CAVITY RESTRAINT PENALTY
#################################
SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION"""

    with open(out_prm_fn, 'w') as w:
        w.writelines(lines)
    return

def prepare_pocket(
        in_poc_fn,
        out_poc_fn, # should be .mol2 format
        verbose=False
        ):
    #if os.path.exists(out_poc_fn):
    #    return 0
    assert in_poc_fn.split('.')[-1] in AVAIL_IN_FORMAT_LIST, "Wrong input format" 
    assert out_poc_fn.split('.')[-1] == "mol2", "Wrong output format"
    command = f"obabel {in_poc_fn} -O {out_poc_fn} "
    if not verbose:
        command += ">& /dev/null"
    else:
        print(command)
    tag = os.system(command)
    return tag

def prepare_ligand(
        ref_lig_fn,
        in_lig_fn,
        out_lig_fn, 
        smarts,
        verbose=False
        ):
    #if os.path.exists(out_lig_fn):
    #    return 0
    assert in_lig_fn.split('.')[-1] in AVAIL_IN_FORMAT_LIST, "Wrong input format" 
    assert ref_lig_fn.split('.')[-1] in AVAIL_IN_FORMAT_LIST, "Wrong ref format" 
    assert out_lig_fn.split('.')[-1] == "sdf", "Wrong output format"
    command = f"sdtether {ref_lig_fn} {in_lig_fn} {out_lig_fn} '{smarts}' "
    if not verbose:
        command += ">& /dev/null"
    else:
        print(command)
    tag = os.system(command)
    return tag

def make_grid(
        in_prm_fn,
        verbose=False
        ):
    assert in_prm_fn.split('.')[-1] == "prm", "Wrong param format"
    command = f"rbcavity -W -d -r {in_prm_fn} "
    if not verbose:
        command += ">& /dev/null"
    else:
        print(command)
    tag = os.system(command)
    return tag

def run_docking(
        in_lig_fn,
        out_lig_prefix,
        in_prm_fn,
        n,
        seed=None,
        verbose=False
        ):
    assert in_lig_fn.split('.')[-1] == "sdf", "Wrong input format"
    assert in_prm_fn.split('.')[-1] == "prm", "Wrong param format"
    command = f"rbdock -i {in_lig_fn} -o {out_lig_prefix} -r {in_prm_fn} -p dock.prm -n {n} "
    if seed is not None:
        command += "--seed {seed} "
    if not verbose:
        command += ">& /dev/null"
    else:
        print(command)
    tag = os.system(command)
    return tag


if __name__ == "__main__":

    import argparse
    import tempfile
    from time import time
    
    parser = argparse.ArgumentParser()

    parser.add_argument("-p", help="protein in PDB format", type=str, required=True)
    parser.add_argument("-l", help="input ligand in SDF format", type=str, required=True)
    parser.add_argument("-r", help="reference ligand in SDF format", type=str, required=True)
    parser.add_argument("-s", help="SMILES to be tethered", type=str, required=True)
    parser.add_argument("-o", help="output ligand prefix", type=str, required=True)
    parser.add_argument("-n", help="num docking per ligand", type=int, default=1)
    parser.add_argument("--seed", help="seed", type=int, default=None)
    parser.add_argument("-v", help="verbose", action="store_true")

    args = parser.parse_args()

    fd, path = tempfile.mkstemp(suffix=".prm", prefix="SD_tmp_", dir="/tmp")

    try: 
        st = time()
        prepare_ligand(args.r, args.l, args.l, args.s, verbose=args.v)
        prepare_pocket(args.p, args.p[:-3]+"mol2", verbose=args.v)
        if args.v:
            print(f"Prepare inputs done --> Duration: {time() - st:.2f}(s)")

        st = time()
        generate_prm(args.p[:-3]+"mol2", args.r, path, args.s)
        if args.v:
            print(f"Generate param file done --> Duration: {time() - st:.2f}(s)")

        st = time()
        make_grid(path, verbose=args.v)
        if args.v:
            print(f"Making grid done --> Duration: {time() - st:.2f}(s)")

        st = time()
        run_docking(args.l, args.o, path, n=args.n, seed=args.seed, verbose=args.v)
        if args.v:
            print(f"Running docking done --> Duration: {time() - st:.2f}(s)")
    except Exception as e:
        print(e)

    os.close(fd)
    os.unlink(path)
    #os.system(f"rm /tmp/SD_tmp_*")


