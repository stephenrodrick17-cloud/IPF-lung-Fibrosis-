#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT_DIR"

mkdir -p logs data/raw data/processed results/tables results/figures
: > logs/pipeline.log

LOGFILE="$ROOT_DIR/logs/pipeline.log"

log() {
  echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOGFILE"
}

log "Starting IPF matrisome pipeline"

# Step 1
log "Running Step 1: GEO metadata extraction"
Rscript scripts/01_download_geo.R

# Step 2
log "Running Step 2: Limma differential expression"
Rscript scripts/02_deg_limma.R

# Step 3
log "Running Step 3: Matrisome filtering"
python scripts/03_matrisome_filter.py

# Step 4
log "Running Step 4: DEG volcano and heatmap"
Rscript scripts/04_deg_volcano_heatmap.R

# Step 5
log "Running Step 5: Common DEGs intersection"
python scripts/05_common_degs.py

# Step 6
log "Running Step 6: Functional enrichment"
python scripts/06_enrichment.py

# Step 7
log "Running Step 7: PPI network and hub genes"
python scripts/07_ppi_hubs.py

# Step 8
log "Running Step 8: Regulatory networks"
python scripts/08_regulatory_networks.py

# Step 9
log "Running Step 9: CTD chemical-gene interactions"
python scripts/09_ctd_chemicals.py

# Step 10
log "Running Step 10: Validation in GSE135065"
Rscript scripts/10_validation_GSE135065.R

log "Pipeline complete"
