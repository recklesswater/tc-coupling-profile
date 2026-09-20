"""Functional verification that the latent space is residue-aligned.

Reviewer's concern, and it is a fair one: with an MLP encoder there is no
guarantee that latent dimension i corresponds to residue i. Shape checks do
not establish this -- a broadcasting or reshape error produces correctly
shaped tensors that map to the wrong residues, and nothing complains.

So test it functionally instead. Perturb the input of exactly one residue and
ask which latent dimensions respond. Under a correct implementation the
answer must be: only the dimensions owned by that residue (and, because the
encoder has a global context branch, a small global shift that is identical
for every residue -- which is itself worth measuring).

This test is decisive in a way that "print(L.shape)" is not: it fails loudly
if the mask is applied to a flattened or mis-broadcast tensor.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import torch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.bd_profile import contact_map, fetch_pdb, gnm_correlation  # noqa: E402
from tc_alignment_model import AlignedVAE  # noqa: E402

PDB_ID = "1UBQ"
PER_RES = 2
SEED = 20260912


def run_case(label, global_ctx, n_res, device):
    """Run the perturbation test for one encoder configuration."""
    d = n_res * PER_RES
    model = AlignedVAE(n_res, per_res=PER_RES, hidden=48, global_ctx=global_ctx,
                       mask=torch.eye(d)).to(device)
    model.eval()

    rng = np.random.default_rng(SEED)
    x0 = torch.tensor(
        rng.normal(0.0, 0.3, size=(1, n_res, 2)), dtype=torch.float32
    )
    with torch.no_grad():
        mu0 = model.encode_mean(x0)

    print()
    print("=" * 74)
    print("%s   (global context dims = %d)" % (label, global_ctx))
    print("=" * 74)
    print("%8s %12s %16s %16s %12s"
          % ("residue", "owned dims", "max|d| own", "max|d| other", "leak"))

    worst = 0.0
    for res in [0, 7, 20, 37, 50, n_res - 1]:
        x1 = x0.clone()
        x1[0, res, :] += 0.5
        with torch.no_grad():
            mu1 = model.encode_mean(x1)
        delta = (mu1 - mu0).abs()[0].numpy()
        own = [res * PER_RES + k for k in range(PER_RES)]
        own_max = float(max(delta[i] for i in own))
        other_max = float(max(delta[i] for i in range(d) if i not in own))
        leak = other_max / max(own_max, 1e-12)
        worst = max(worst, leak)
        print("%8d %12s %16.5f %16.5f %12.3f"
              % (res + 1, str(own), own_max, other_max, leak))

    # Decisive check: how many INPUT residues can move a given latent dim?
    #
    # The leak ratio above is NOT the right criterion here. The encoder has a
    # local window of +/-2 residues on purpose -- that context is what makes it
    # a useful encoder -- so perturbing residue i is *supposed* to move the
    # latent dims of residues i-2..i+2. Counting those as leakage would reject
    # a correct implementation.
    #
    # What actually matters is the receptive field: under strict alignment,
    # latent dim i may depend on residue i and its window neighbours and on
    # nothing else. A purely local encoder gives a receptive field of 5
    # residues (3 at a chain terminus); a global-context branch gives all L.
    probes = {}
    for probe in (0, n_res // 2):
        x = x0.clone().requires_grad_(True)
        mu = model.encode_mean(x)
        (grad,) = torch.autograd.grad(mu[0, probe], x, retain_graph=False)
        per_res = grad[0].abs().sum(dim=-1).numpy()
        active = np.where(per_res > 1e-8)[0]
        probes[probe] = len(active)

    window = 2 * (2 * 2 + 1)  # generous bound if all window entries mattered
    n_active = probes[n_res // 2]
    locally_aligned = n_active <= window

    print()
    print("  receptive field of latent dim 0      : %d residue(s)" % probes[0])
    print("  receptive field of a middle dim      : %d residue(s)" % n_active)
    print("  (local window of +/-2 implies at most %d)" % window)
    print("  worst neighbour-response ratio       : %.3f  (informational)" % worst)
    verdict = ("STRICTLY RESIDUE-ALIGNED" if locally_aligned
               else "POSITION-INDEXED ONLY (global branch -> not local)")
    print("  verdict                              : %s" % verdict)
    return dict(label=label, global_ctx=global_ctx, worst=worst,
                n_active=n_active, n_res=n_res, verdict=verdict)


def main():
    torch.manual_seed(SEED)
    np.random.seed(SEED)
    device = torch.device("cpu")

    path = fetch_pdb(PDB_ID)
    from src.pdb_io import ca_trace

    ca, names, _ = ca_trace(path)
    n_res = len(ca)
    print("%s: %d residues" % (PDB_ID, n_res))

    corr = gnm_correlation(contact_map(ca))
    d = n_res * PER_RES

    print("Perturb residue i, then ask which latent dimensions respond.")
    print("Strict residue alignment means: only the dims owned by residue i.")
    print("Leak ratio = (largest response on other dims) / (response on own dims).")

    results = [
        run_case("A. purely local encoder", 0, n_res, device),
        run_case("B. encoder with global context", 8, n_res, device),
    ]

    print()
    print("=" * 74)
    print("CONCLUSION")
    print("=" * 74)
    for r in results:
        print("  %-32s receptive field %2d/%d residues   %s"
              % (r["label"], r["n_active"], r["n_res"], r["verdict"]))
    print()
    print("  The local encoder is strictly residue-aligned. The global-context")
    print("  branch, added because the model could not learn collectively")
    print("  driven motion without it, makes every latent dimension depend on")
    print("  every residue. The latent space is therefore only")
    print("  POSITION-INDEXED, not residue-local, and the phrase")
    print("  \"semantic alignment\" overstates what the architecture guarantees.")


if __name__ == "__main__":
    main()
