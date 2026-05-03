"""
plot_EDAC_C.py – legge bch_results.csv e plotta il grafico BCH correction rate.

Uso:
    # Compila e lancia il test (genera bch_results.csv)
    g++ -std=c++17 -O2 -DBCH_NO_MAIN -o test_bch test_bch.cpp gf2_poly.cpp
    ./test_bch

    # Poi plotta
    python3 plot.py
"""

import csv
import sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

CSV_FILE = "bch_results.csv"
FLOOR    = 5e-5   # valore minimo mostrato sull'asse log per gli zero esatti

# ── Legge i dati ──────────────────────────────────────────────────────────────
try:
    with open(CSV_FILE) as f:
        reader = csv.DictReader(f)
        rows   = list(reader)
except FileNotFoundError:
    print(f"Errore: {CSV_FILE} non trovato. Esegui prima ./test_bch")
    sys.exit(1)

flips = [int(r["flips"])        for r in rows]
rates = [float(r["success_rate"]) for r in rows]

# Per la scala log gli zero diventano FLOOR (ma li marchiamo diversamente)
is_zero    = [r == 0.0 for r in rates]
plot_rates = [FLOOR if z else r for r, z in zip(rates, is_zero)]

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 5))
ax.set_title("BCH Correction Success Rate vs. Injected Bit-flips", fontsize=13)

# Linea e punti normali
ax.plot(flips, plot_rates, color="#1A5F9E", linewidth=1.8, zorder=2)
ax.scatter(
    [f for f, z in zip(flips, is_zero) if not z],
    [r for r, z in zip(plot_rates, is_zero) if not z],
    color="#1A5F9E", s=40, zorder=3, label="T=2"
)

# Punti a zero: triangolo verso il basso per indicare "fuori scala"
ax.scatter(
    [f for f, z in zip(flips, is_zero) if z],
    [FLOOR * 1.5 for z in is_zero if z],
    marker="v", color="#1A5F9E", s=60, zorder=3, label="0 (zero esatto)"
)

ax.set_yscale("log")
ax.set_ylim(FLOOR * 0.5, 5)
ax.yaxis.set_major_formatter(ticker.LogFormatterSciNotation(base=10))
ax.set_xlabel("Number of Injected Bit-flips", fontsize=11)
ax.set_ylabel("Success Rate [-]",             fontsize=11)
ax.set_xticks(flips)
ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.6)
ax.legend(fontsize=10)

plt.tight_layout()
# plt.savefig("bch_plot.png", dpi=150)
# print("Salvato: bch_plot.png")
plt.show()
