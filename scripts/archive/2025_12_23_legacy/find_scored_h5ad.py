# Without running out of memory, it helps you extract the files from hundreds or thousands of .h5ad files 
# that contain key analytical data such as fate_prob or pseudotime, 
# and sorts them in order of information richness.
import os, glob
import scanpy as sc

root = r"E:\scgpt"
keys = ["p_iPSC","p_target","score_target","fate_prob","wot_pseudotime","pseudotime","latent_time"]

files = glob.glob(os.path.join(root, "**", "*.h5ad"), recursive=True)

out = []
for f in files:
    try:
        a = sc.read_h5ad(f, backed="r")
        cols = set(a.obs_keys())
        hit = [k for k in keys if k in cols]
        if hit:
            out.append((len(hit), hit, f))
    except Exception:
        pass

out = sorted(out, key=lambda x: (-x[0], x[2]))
print("Found candidates:", len(out))
for n, hit, path in out[:30]:
    print(n, hit, path)
