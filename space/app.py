"""Gradio demo: coupling indicator (TC) computed from a PDB structure.

Anyone can type a PDB ID and see, in about a second, which parts of the chain
are the ones where treating residues as independently moving is a poor
approximation.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from pdb_io import burial_proxy, ca_trace
from tc_profile import contact_map, fetch_pdb, gnm_correlation, tc_profile

EXAMPLE_IDS = ["1UBQ", "1CRN", "1PGB", "3CHY", "1TIT"]
WINDOW = 7


def compute(pdb_id: str):
    """Return (figure, summary text) for one PDB id."""
    pdb_id = (pdb_id or "").strip().upper()
    if len(pdb_id) != 4 or not pdb_id.isalnum():
        return None, "Please enter a 4-character PDB ID, for example 1UBQ."

    try:
        path = fetch_pdb(pdb_id)
        ca, names, _ = ca_trace(path)
    except Exception as exc:  # noqa: BLE001
        return None, "Could not fetch or parse %s (%s)." % (pdb_id, repr(exc)[:80])

    if len(ca) < 30:
        return None, ("%s has only %d resolved C-alpha atoms in its first chain; "
                      "the indicator needs a longer chain." % (pdb_id, len(ca)))

    contacts = contact_map(ca)
    corr = gnm_correlation(contacts)
    prof = tc_profile(corr, WINDOW)
    burial = burial_proxy(ca, cutoff=10.0)

    n = len(ca)
    resids = np.arange(1, n + 1)
    burial_n = (burial - burial.min()) / max(1e-9, burial.max() - burial.min())
    peak = int(np.nanargmax(prof))

    fig, axes = plt.subplots(1, 3, figsize=(14, 3.9))

    ax = axes[0]
    im = ax.imshow(contacts, cmap="Greys", interpolation="nearest")
    ax.set_title("contact map (C-alpha < 8 A)\nmean degree %.1f" % contacts.sum(1).mean())
    ax.set_xlabel("residue")
    ax.set_ylabel("residue")
    fig.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[1]
    im = ax.imshow(corr, cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
    ax.set_title("GNM correlation matrix")
    ax.set_xlabel("residue")
    ax.set_ylabel("residue")
    fig.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[2]
    ax.plot(resids, prof, "-o", ms=3, color="#2f855a", label="TC (window %d)" % WINDOW)
    ax.plot(resids, burial_n * float(np.nanmax(prof)), "-", lw=1.1,
            color="#2b6cb0", alpha=0.7, label="burial proxy (scaled)")
    ax.axvspan(max(1, peak + 1 - 3), min(n, peak + 1 + 3), color="#e53e3e", alpha=0.13)
    ax.set_xlabel("residue number")
    ax.set_ylabel("TC (nats)")
    ax.set_title("coupling indicator")
    ax.legend(fontsize=8, loc="lower left")
    ax.grid(alpha=0.25)

    fig.suptitle("%s  -  %d residues in the first chain" % (pdb_id, n), fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.93])

    summary = (
        "**%s** - %d residues\n\n"
        "- TC range: **%.3f - %.3f** nats (mean %.3f)\n"
        "- highest coupling around residue **%d** (%s)\n"
        "- residues with no contacts in their window can give a near-singular "
        "block and are reported as NaN\n\n"
        "High TC = the residues in that window move together, so treating them "
        "as independent costs the most. This is a *structural indicator*, not a "
        "binding-site predictor - on five holo structures, ligand-binding "
        "residues had **lower** TC than the chain average."
        % (pdb_id, n, float(np.nanmin(prof)), float(np.nanmax(prof)),
           float(np.nanmean(prof)), peak + 1, names[peak])
    )
    return fig, summary


def build_demo():
    import gradio as gr

    with gr.Blocks(title="TC coupling profile") as demo:
        gr.Markdown(
            "# TC coupling profile\n"
            "**Mean field is almost always good enough. This shows where "
            "\"almost\" fails.**\n\n"
            "Type a PDB ID. The indicator is the block total correlation of the "
            "Gaussian Network Model correlation matrix,\n\n"
            "`TC(B) = -1/2 ln det R_BB`\n\n"
            "computed per residue window straight from the structure - no "
            "simulation, no training, no experimental data."
        )
        with gr.Row():
            inp = gr.Textbox(value="1UBQ", label="PDB ID", max_lines=1)
            btn = gr.Button("Compute", variant="primary")
        plot = gr.Plot(label="structure and profile")
        out = gr.Markdown()

        gr.Examples(examples=[[i] for i in EXAMPLE_IDS], inputs=[inp])
        gr.Markdown(
            "Source, caveats and the multi-protein validation: "
            "[github.com/recklesswater/tc-coupling-profile]"
            "(https://github.com/recklesswater/tc-coupling-profile)\n\n"
            "If a structure fails to download, RCSB may be intermittently "
            "unreachable from this Space."
        )

        btn.click(compute, inputs=inp, outputs=[plot, out])
        inp.submit(compute, inputs=inp, outputs=[plot, out])
    return demo


if __name__ == "__main__":
    build_demo().launch()
