import os


AVAIL_IN_FORMAT_LIST = ["sdf", "pdb", "mol2", "sd", "xyz"]


def generate_prm(
        ref_poc_fn,
        ref_lig_fn,
        out_prm_fn,
        title=""
        ):
    assert ref_poc_fn.split('.')[-1] == "mol2", "Wrong receptor format" 
    assert ref_lig_fn.split('.')[-1] == "sd", "Wrong ligand format"
    assert out_prm_fn.split('.')[-1] == "prm", "Wrong param format"
    
    lines = f"""RBT_PARAMETER_FILE_V1.00
TITLE {title}

RECEPTOR_FILE {ref_poc_fn}

SECTION LIGAND
    TRANS_MODE TETHERED
    ROT_MODE TETHERED 
    DIHEDRAL_MODE TETHERED 
    MAX_TRANS 1.0 
    MAX_ROT 10.0 
    MAX_DIHEDRAL 10.0
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
        out_poc_fn # should be .mol2 format
        ):
    if os.path.exists(out_poc_fn):
        return 0
    assert in_poc_fn.split('.')[-1] in AVAIL_IN_FORMAT_LIST, "Wrong input format" 
    assert out_poc_fn.split('.')[-1] == "mol2", "Wrong output format"
    command = f"obabel {in_poc_fn} -O {out_poc_fn}"
    print(command)
    tag = os.system(command)
    return tag

def prepare_ligand(
        ref_lig_fn,
        in_lig_fn,
        out_lig_fn, # should be .sd format
        smarts
        ):
    if os.path.exists(out_lig_fn):
        return 0
    assert in_lig_fn.split('.')[-1] in AVAIL_IN_FORMAT_LIST, "Wrong input format" 
    assert ref_lig_fn.split('.')[-1] in AVAIL_IN_FORMAT_LIST, "Wrong ref format" 
    assert out_lig_fn.split('.')[-1] == "sd", "Wrong output format"
    command = f"sdtether {ref_lig_fn} {in_lig_fn} {out_lig_fn} '{smarts}'"
    print(command)
    tag = os.system(command)
    return tag

def make_grid(
        in_prm_fn
        ):
    assert in_prm_fn.split('.')[-1] == "prm", "Wrong param format"
    command = f"rbcavity -W -d -r {in_prm_fn}"
    print(command)
    tag = os.system(command)
    return tag

def run_docking(
        in_lig_fn,
        out_lig_prefix,
        in_prm_fn,
        n=1,
        seed=0,
        ):
    assert in_lig_fn.split('.')[-1] == "sd", "Wrong input format"
    assert in_prm_fn.split('.')[-1] == "prm", "Wrong param format"
    command = f"rbdock -i {in_lig_fn} -o {out_lig_prefix} -r {in_prm_fn} -p dock.prm -n {n} --seed {seed}"
    print(command)
    tag = os.system(command)
    return tag


if __name__ == "__main__":

    import argparse
    from time import time
    
    parser = argparse.ArgumentParser()

    parser.add_argument("-p", help="protein_fn in PDB format", type=str)
    parser.add_argument("-l", help="ligand_fn in SDF format", type=str)
    parser.add_argument("-c", help="complex_fn in PDB format", type=str)
    parser.add_argument("-t", help="traj_fn in DCD format", type=str)
    parser.add_argument("-o", help="log_fn", type=str)
    parser.add_argument("-v", help="verbose", action="store_true")

    args = parser.parse_args()

    KEY = "1jzs"
    N = 1
    lig_ori_fn = f"files/{KEY}_ligand_0.sdf"
    lig_ref_fn = f"files/{KEY}_ligand.sd"
    poc_ori_fn = f"files/{KEY}_protein_cleaned.pdb"
    scaff_smi = "C1COCC(CC2CO2)C1"

    lig_in_fn = f"files/{KEY}_ligand_input.sd"
    lig_out_prefix = f"files/{KEY}_ligand_dock"
    poc_in_fn = f"files/{KEY}_protein_cleaned.mol2"

    prm_fn = f"files/{KEY}.prm"

    st = time()
    prepare_ligand(lig_ref_fn, lig_ori_fn, lig_in_fn, scaff_smi)
    prepare_pocket(poc_ori_fn, poc_in_fn)
    print(f"{time() - st:.2f}(s)")

    generate_prm(poc_in_fn, lig_ref_fn, prm_fn, "test")
    print(f"{time() - st:.2f}(s)")

    make_grid(prm_fn)
    print(f"{time() - st:.2f}(s)")

    run_docking(lig_in_fn, lig_out_prefix, prm_fn, n=N)
    print(f"{time() - st:.2f}(s)")



