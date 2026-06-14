#!/usr/bin/env bash
# Render reports/model_scorecard.md -> reports/model_scorecard.pdf via Quarto+typst
# (no LaTeX install needed). Run after a pipeline run refreshes the scorecard.
set -euo pipefail
cd "$(dirname "$0")/.."
tmp="$(mktemp -d)/scorecard.qmd"
cat > "$tmp" <<'YAML'
---
format:
  typst:
    papersize: us-letter
    margin:
      x: 1.4cm
      y: 1.6cm
    fontsize: 9pt
    toc: true
    toc-title: Contents
    toc-depth: 2
---
YAML
cat reports/model_scorecard.md >> "$tmp"
quarto render "$tmp" --to typst --quiet
cp "${tmp%.qmd}.pdf" reports/model_scorecard.pdf
echo "wrote reports/model_scorecard.pdf"
