"""Re-verify the RDKit ETKDGv3 baseline numbers quoted in the design (section 1):
   total 86.5 s / mean 9.8 ms / slowest 8.11 s / 46 failures / 1 C++ exception.
Single process, seed 0xf00d, AddHs, as measure_params.py does.
"""
import sys
import time

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

CORPUS = "/Users/tom/Projects/momega/omgkit/.claude/worktrees/agent-a92e349554d652a6b/harness/corpus/large.smi"

smis = []
for line in open(CORPUS):
    line = line.strip()
    if line and not line.startswith("#"):
        smis.append(line.split("\t")[0])
print(f"corpus lines (non-comment): {len(smis)}", flush=True)

fail = 0
exc = 0
parse_fail = 0
times = []
worst = []
t0 = time.time()
for i, smi in enumerate(smis):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        parse_fail += 1
        continue
    mh = Chem.AddHs(m)
    p = AllChem.ETKDGv3()
    p.randomSeed = 0xF00D
    t = time.time()
    try:
        rc = AllChem.EmbedMolecule(mh, p)
        if rc < 0:
            fail += 1
    except Exception as e:  # noqa: BLE001
        exc += 1
        fail += 1
        print(f"  EXCEPTION on #{i}: {smi[:70]}\n     {str(e)[:200]}", flush=True)
    dt = time.time() - t
    times.append(dt)
    worst.append((dt, smi))
tot = time.time() - t0
times.sort()
print(f"total wall: {tot:.1f} s   embeds: {len(times)}   "
      f"mean {1000*sum(times)/len(times):.1f} ms   median {1000*times[len(times)//2]:.1f} ms")
print(f"p99 {1000*times[int(0.99*len(times))]:.0f} ms   max {times[-1]:.2f} s")
print(f"parse failures: {parse_fail}   embed failures: {fail} "
      f"({100.0*fail/len(times):.2f}%)   C++ exceptions: {exc}")
worst.sort(key=lambda w: -w[0])
print("slowest 5:")
for dt, smi in worst[:5]:
    print(f"   {dt:.2f} s  {smi[:80]}")
