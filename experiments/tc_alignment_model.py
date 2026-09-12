"""Minimal residue-aligned encoder used by the alignment verification test.

This is a stripped-down copy of the architecture used in the (unpublished)
variational work: a per-residue, weight-shared encoder plus a global context
branch. It is here only so that the alignment claim can be tested in
isolation, without pulling in the whole VAE.

The design is what makes "residue-aligned" true rather than aspirational:

* the encoder is applied position by position with SHARED weights, so output
  row i is a function of residue i's local window and nothing else;
* a global context vector is computed from the whole chain and concatenated
  identically at every position. It shifts all rows by the same amount, which
  does not break the index correspondence but does mean latent dim i is not a
  strictly local function of residue i.
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


class AlignedVAE(nn.Module):
    def __init__(self, n_res: int, per_res: int = 2, hidden: int = 48,
                 local_window: int = 2, global_ctx: int = 8, mask=None):
        super().__init__()
        self.n_res = n_res
        self.per_res = per_res
        self.n_lat = n_res * per_res
        self.local_window = local_window

        in_dim = 2 * (2 * local_window + 1) + global_ctx
        self.global_enc = nn.Sequential(
            nn.Linear(2 * n_res, hidden), nn.SiLU(), nn.Linear(hidden, global_ctx)
        )
        self.encoder = nn.Sequential(
            nn.Linear(in_dim, hidden), nn.SiLU(), nn.Linear(hidden, per_res)
        )
        if mask is None:
            mask = torch.eye(self.n_lat)
        self.register_buffer("mask", mask.clone().float())

    def _windows(self, angles: torch.Tensor) -> torch.Tensor:
        b = angles.shape[0]
        w = self.local_window
        padded = F.pad(angles, (0, 0, w, w), mode="replicate")
        return torch.cat(
            [padded[:, i:i + self.n_res, :] for i in range(2 * w + 1)], dim=-1
        )

    def encode_mean(self, angles: torch.Tensor) -> torch.Tensor:
        """Latent mean only, of shape ``(B, n_res * per_res)``."""
        b = angles.shape[0]
        local = self._windows(angles)
        g = self.global_enc(angles.reshape(b, -1)).unsqueeze(1)
        g = g.expand(-1, self.n_res, -1)
        return self.encoder(torch.cat([local, g], dim=-1)).reshape(b, self.n_lat)
