from openff.toolkit import Topology
from openff.toolkit.utils import get_data_file_path
import contextlib

def successfully_loaded(top):
    match_info = [atom.metadata["match_info"] for atom in top.atoms]
    return all([bool(match) for match in match_info])

#____________SIMPLE PROTEIN TESTS (should already be implemented)______________
top = Topology.from_pdb(get_data_file_path("proteins/TwoMol_SER_CYS.pdb"))
assert successfully_loaded(top)

top = Topology.from_pdb(get_data_file_path("proteins/T4-protein.pdb"))
assert successfully_loaded(top)
# #____________________SIMPLE POLYMERS (custom substructures)____________________
rubber_substructs = {
                        "naturalrubber": "[#6D4+0:16](-[#6D3+0:17](=[#6D3+0:18](-[#6D4+0:19](-[#1D1+0:20])(-[#1D1+0:21])-[*:22])-[#6D4+0:23](-[#1D1+0:24])(-[#1D1+0:25])-[#1D1+0:26])-[#1D1+0:27])(-[#1D1+0:28])(-[#1D1+0:29])-[*:30]",
                        "naturalrubber_TERM1": "[#6D4+0:16](-[#6D3+0:17](=[#6D3+0:18](-[#6D4+0:19](-[#1D1+0:20])(-[#1D1+0:21])-[#1D1+0:22])-[#6D4+0:23](-[#1D1+0:24])(-[#1D1+0:25])-[#1D1+0:26])-[#1D1+0:27])(-[#1D1+0:28])(-[#1D1+0:29])-[*:30]",
                        "naturalrubber_TERM2": "[#6D4+0:16](-[#6D3+0:17](=[#6D3+0:18](-[#6D4+0:19](-[#1D1+0:20])(-[#1D1+0:21])-[*:22])-[#6D4+0:23](-[#1D1+0:24])(-[#1D1+0:25])-[#1D1+0:26])-[#1D1+0:27])(-[#1D1+0:28])(-[#1D1+0:29])-[#1D1+0:30]"
                    }
top = Topology.from_pdb("naturalrubber.pdb", _custom_substructures=rubber_substructs)
assert successfully_loaded(top)

PE_substructs = {
                    "PE": "[#6D4+0:9](-[#1D1+0:10])(-[#1D1+0:11])(-[#6D4+0:12](-[#1D1+0:13])(-[#1D1+0:14])-[*:15])-[*:16]",
                    "PE_TERM1": "[#6D4+0:9](-[#1D1+0:10])(-[#1D1+0:11])(-[#6D4+0:12](-[#1D1+0:13])(-[#1D1+0:14])-[#1D1+0:15])-[*:16]"
                }
top = Topology.from_pdb("polyethylene.pdb", _custom_substructures=PE_substructs)
assert successfully_loaded(top)

pnipam_substructs = {
                        "pnipam": "[#6D4+0:22](-[#1D1+0:23])(-[#1D1+0:24])(-[#6D4+0:25](-[#1D1+0:26])(-[#6D3+0:27](=[#8D1+0:28])-[#7D3+0:29](-[#1D1+0:30])-[#6D4+0:31](-[#1D1+0:32])(-[#6D4+0:33](-[#1D1+0:34])(-[#1D1+0:35])-[#1D1+0:36])-[#6D4+0:37](-[#1D1+0:38])(-[#1D1+0:39])-[#1D1+0:40])-[*:41])-[*:42]",
                        "pnipam_TERM2": "[#6D4+0:22](-[#1D1+0:23])(-[#1D1+0:24])(-[#6D4+0:25](-[#1D1+0:26])(-[#6D3+0:27](=[#8D1+0:28])-[#7D3+0:29](-[#1D1+0:30])-[#6D4+0:31](-[#1D1+0:32])(-[#6D4+0:33](-[#1D1+0:34])(-[#1D1+0:35])-[#1D1+0:36])-[#6D4+0:37](-[#1D1+0:38])(-[#1D1+0:39])-[#1D1+0:40])-[#1D1+0:41])-[*:42]",
                        "pnipam_TERM3": "[#6D4+0:22](-[#1D1+0:23])(-[#1D1+0:24])(-[#6D4+0:25](-[#1D1+0:26])(-[#6D3+0:27](=[#8D1+0:28])-[#7D3+0:29](-[#1D1+0:30])-[#6D4+0:31](-[#1D1+0:32])(-[#6D4+0:33](-[#1D1+0:34])(-[#1D1+0:35])-[#1D1+0:36])-[#6D4+0:37](-[#1D1+0:38])(-[#1D1+0:39])-[#1D1+0:40])-[*:41])-[#1D1+0:42]"
                    }
top = Topology.from_pdb("pnipam_modified.pdb", _custom_substructures=pnipam_substructs)
assert successfully_loaded(top)

#______________FUNCTIONALIZED PROTEIN (existing + custom substructures)________

five_xg9_substructs = {
                        "peg": "[#6D4+0:13](-[#1D1+0:14])(-[#1D1+0:15])(-[#6D4+0:16](-[#1D1+0:17])(-[#1D1+0:18])-[#8D2+0:19]-[*:20])-[*:21]",
                        "peg_TERM1": "[#6D4+0:16](-[#1D1+0:17])(-[#1D1+0:18])(-[#6D4+0:19](-[#1D1+0:20])(-[#1D1+0:21])-[#8D2+0:22]-[#1D1+0:23])-[*:27]",
                        "peg_TERM2": "[#6D4+0:16](-[#1D1+0:17])(-[#1D1+0:18])(-[#6D4+0:19](-[#1D1+0:20])(-[#1D1+0:21])-[#8D2+0:22]-[*:23])-[#1D1+0:27]",
                        "H2O2": "[#8D2+0:1](-[#1D1+0:2])-[#8D2+0:1](-[#1D1+0:3])",
                        "TRP_neg": "[#1D1+0:1]-[#6D3+0:2]1-[#7D2-1:3]-[#6D3+0:4]2-[#6D3+0:5](-[#1D1+0:6])=[#6D3+0:7](-[#1D1+0:8])-[#6D3+0:9](-[#1D1+0:10])=[#6D3+0:11](-[#1D1+0:12])-[#6D3+0:13]=2-[#6D3+0:14]=1-[#6D4+0:15](-[#1D1+0:16])(-[#1D1+0:17])-[#6D4+0:18]([#1D1+0:19])([#7D3+0:20]([#1D1+0:24])([*:25]))([#6D3+0:21](=[#8D1+0:22])([*:23]))",
                        "LEU_Ald_TERM": "[#7D3+0:1]([#1D1+0:2])([*:3])-[#6D4+0:4]([#1D1+0:5])([#6D4+0:6]([#1D1+0:7])([#1D1+0:8])([#1D1+0:9]))-[#6D3+0:10](-[#1D1+0:11])=[#8D1+0:12]"
                      }
top = Topology.from_pdb("5xg9_hydrogenated_fixed_bond.pdb", _custom_substructures=five_xg9_substructs)

assert successfully_loaded(top)
print("All files run successfully. Beginning Tests:")

#________________________________ERROR TESTING_________________________________
@contextlib.contextmanager
def report_run(run_name):
    print(f"Starting Run: {run_name}")
    print("-"*50)
    try:
        yield
    except Exception as e:
        print(e)
    finally:
        print("-"*50, end="\n\n")

with report_run("Incorrect Atomic Number Spec"):
    PE_substructs = {"PE": "[CD4+0:9](-[#1D1+0:10])(-[#1D1+0:11])(-[#6D4+0:12](-[#1D1+0:13])(-[#1D1+0:14])-[*:15])-[*:16]"}
    top = Topology.from_pdb("polyethylene.pdb", _custom_substructures=PE_substructs)

with report_run("Incorrect Connectivity"):
    PE_substructs = {"PE": "[#6D3+0:9](-[#1D1+0:10])(-[#1D1+0:11])(-[#6D4+0:12](-[#1D1+0:13])(-[#1D1+0:14])-[*:15])-[*:16]"}
    top = Topology.from_pdb("polyethylene.pdb", _custom_substructures=PE_substructs)

with report_run("No Formal Charge"):
    PE_substructs = {"PE": "[#6D4+0:9](-[#1D1+0:10])(-[#1D1+0:11])(-[#6D4:12](-[#1D1+0:13])(-[#1D1+0:14])-[*:15])-[*:16]"}
    top = Topology.from_pdb("polyethylene.pdb", _custom_substructures=PE_substructs)

with report_run("Unspecified Bond Type"):
    PE_substructs = {"PE": "[#6D4+0:9](~[#1D1+0:10])(-[#1D1+0:11])(-[#6D4+0:12](-[#1D1+0:13])(-[#1D1+0:14])-[*:15])-[*:16]"}
    top = Topology.from_pdb("polyethylene.pdb", _custom_substructures=PE_substructs)

with report_run("Duplicate Entry Names"):
    PE_substructs = {"LEU": "[#6D4+0:9](-[#1D1+0:10])(-[#1D1+0:11])(-[#6D4+0:12](-[#1D1+0:13])(-[#1D1+0:14])-[*:15])-[*:16]"}
    top = Topology.from_pdb("polyethylene.pdb", _custom_substructures=PE_substructs)

with report_run("Ambiguous Bond Info Assigned"):
    PE_substructs = {"PE": "[#6D4+0:9](-[#1D1+0:10])(-[#1D1+0:11])(-[#6D4+0:12](-[#1D1+0:13])(-[#1D1+0:14])-[*:15])-[*:16]",
                     "PE_wrong": "[#6D4+0:9](-[#1D1+0:10])(-[#1D1+0:11])(=[#6D4+0:12](-[#1D1+0:13])(-[#1D1+0:14])-[*:15])-[*:16]"
                     }
    top = Topology.from_pdb("polyethylene.pdb", _custom_substructures=PE_substructs)   

with report_run("Ambiguous Bond Info Assigned"):
    PE_substructs = {"PE": "[#6D4+0:9](-[#1D1+0:10])(-[#1D1+0:11])(-[#6D4+0:12](-[#1D1+0:13])(-[#1D1+0:14])-[*:15])-[*:16]",
                     "PE_wrong": "[#6D4+0:9](-[#1D1+0:10])(-[#1D1+0:11])(-[#6D4+1:12](-[#1D1+0:13])(-[#1D1+0:14])-[*:15])-[*:16]"
                     }
    top = Topology.from_pdb("polyethylene.pdb", _custom_substructures=PE_substructs)  