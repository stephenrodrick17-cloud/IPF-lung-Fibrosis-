Add-Type -AssemblyName System.Drawing

$projectDir = "d:/2ND PROJECT/project"
Set-Location $projectDir

# Ensure output directory exists
if (-not (Test-Path "results/figures/deg")) {
    New-Item -ItemType Directory -Force -Path "results/figures/deg" | Out-Null
}

# ==============================================================================
# PART 1: 3-PANEL VOLCANO PLOTS (GSE32539, GSE110147, GSE24206)
# ==============================================================================
$bmpWidth = 2400
$bmpHeight = 900
$bmp = New-Object System.Drawing.Bitmap($bmpWidth, $bmpHeight)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::AntiAliasGridFit
$g.Clear([System.Drawing.Color]::White)

$fontTitle = New-Object System.Drawing.Font("Arial", 16, [System.Drawing.FontStyle]::Bold)
$fontSub = New-Object System.Drawing.Font("Arial", 11, [System.Drawing.FontStyle]::Regular)
$fontAxis = New-Object System.Drawing.Font("Arial", 11, [System.Drawing.FontStyle]::Bold)
$fontTicks = New-Object System.Drawing.Font("Arial", 9, [System.Drawing.FontStyle]::Regular)
$fontLabel = New-Object System.Drawing.Font("Arial", 8, [System.Drawing.FontStyle]::Bold)
$fontLegend = New-Object System.Drawing.Font("Arial", 10, [System.Drawing.FontStyle]::Regular)

$penAxis = New-Object System.Drawing.Pen([System.Drawing.Color]::Black, 2)
$penDash = New-Object System.Drawing.Pen([System.Drawing.Color]::DarkGray, 1.5)
$penDash.DashStyle = [System.Drawing.Drawing2D.DashStyle]::Dash

$brushUp = New-Object System.Drawing.SolidBrush([System.Drawing.Color]::FromArgb(210, 220, 53, 69))    # Crimson Red
$brushDown = New-Object System.Drawing.SolidBrush([System.Drawing.Color]::FromArgb(210, 13, 110, 253))  # Royal Blue
$brushNS = New-Object System.Drawing.SolidBrush([System.Drawing.Color]::FromArgb(90, 180, 180, 180))    # Light Grey
$brushText = [System.Drawing.Brushes]::Black

$cohorts = @(
    @{ id="GSE32539"; title="GSE32539 (Agilent 4x44K, N=168)"; n_up=103; n_down=59; total=162; topGenes=@("MMP7","ASPN","COL1A1","COMP","SULF1","TIMP3","IL6","S100A12","S100A8") },
    @{ id="GSE110147"; title="GSE110147 (Affymetrix U133+, N=33)"; n_up=164; n_down=212; total=376; topGenes=@("MMP7","SFRP2","CXCL14","FNDC1","COL1A1","TIMP3","S100A8","CLEC4E") },
    @{ id="GSE24206"; title="GSE24206 (Affymetrix Gene 1.0, N=15)"; n_up=97; n_down=38; total=135; topGenes=@("MMP7","ASPN","COL3A1","CXCL14","COL1A1","TIMP3","S100A12","ADAMTS8") }
)

$panelW = 700
$panelH = 650
$marginTop = 140

for ($i = 0; $i -lt 3; $i++) {
    $c = $cohorts[$i]
    $panelX = 80 + ($i * 760)
    $panelY = $marginTop
    
    $limmaFile = "results/tables/$($c.id)_limma_full.csv"
    $modeAFile = "results/tables/$($c.id)_matrisome_DEGs_modeA.csv"
    
    $limma = Import-Csv $limmaFile
    $modeA = Import-Csv $modeAFile
    $upSet = @($modeA | Where-Object { [double]$_.logFC -gt 0 } | ForEach-Object { $_.gene })
    $downSet = @($modeA | Where-Object { [double]$_.logFC -lt 0 } | ForEach-Object { $_.gene })
    
    # Title
    $g.DrawString($c.title, $fontTitle, $brushText, [float]($panelX + 20), [float]30)
    $subtext = "Matrisome DEGs: $($c.total) (|log2FC| > 0.585, adj.P < 0.05)`nUpregulated: $($c.n_up) | Downregulated: $($c.n_down)"
    $g.DrawString($subtext, $fontSub, [System.Drawing.Brushes]::DimGray, [float]($panelX + 20), [float]60)
    
    $xMin = -6.0; $xMax = 6.0
    $yMin = 0.0; $yMax = 35.0
    if ($c.id -eq "GSE24206") { $yMax = 15.0 }
    
    $rcX = [int]($panelX + 70)
    $rcY = [int]$panelY
    $rcW = [int]($panelW - 100)
    $rcH = [int]($panelH - 80)
    
    $g.FillRectangle([System.Drawing.Brushes]::WhiteSmoke, $rcX, $rcY, $rcW, $rcH)
    $g.DrawRectangle($penAxis, $rcX, $rcY, $rcW, $rcH)
    
    $toX = { param($val) return [float]($rcX + (($val - $xMin) / ($xMax - $xMin)) * $rcW) }
    $toY = { param($val) return [float]($rcY + $rcH - (($val - $yMin) / ($yMax - $yMin)) * $rcH) }
    
    # Dashed cutoff lines
    $g.DrawLine($penDash, [float](& $toX 0.585), [float]$rcY, [float](& $toX 0.585), [float]($rcY + $rcH))
    $g.DrawLine($penDash, [float](& $toX -0.585), [float]$rcY, [float](& $toX -0.585), [float]($rcY + $rcH))
    $yCut = -[math]::Log10(0.05)
    $g.DrawLine($penDash, [float]$rcX, [float](& $toY $yCut), [float]($rcX + $rcW), [float](& $toY $yCut))
    
    # Draw points
    foreach ($row in $limma) {
        $lfc = [double]$row.logFC
        $pval = [double]$row.P.Value
        if ($pval -le 1e-45) { $pval = 1e-45 }
        $nlp = -[math]::Log10($pval)
        $gn = $row.gene
        
        if ($nlp -gt $yMax) { $nlp = $yMax }
        if ($lfc -lt $xMin) { $lfc = $xMin }
        if ($lfc -gt $xMax) { $lfc = $xMax }
        
        $px = & $toX $lfc
        $py = & $toY $nlp
        
        $brush = $brushNS
        $size = 3
        if ($upSet -contains $gn) {
            $brush = $brushUp
            $size = 6
        } elseif ($downSet -contains $gn) {
            $brush = $brushDown
            $size = 6
        }
        
        $g.FillEllipse($brush, [float]($px - $size/2), [float]($py - $size/2), [float]$size, [float]$size)
        
        if ($c.topGenes -contains $gn) {
            $g.FillEllipse([System.Drawing.Brushes]::Black, [float]($px - 3), [float]($py - 3), [float]6, [float]6)
            $g.DrawString($gn, $fontLabel, [System.Drawing.Brushes]::Black, [float]($px + 4), [float]($py - 5))
        }
    }
    
    # Ticks
    for ($tickX = -6; $tickX -le 6; $tickX += 2) {
        $tx = & $toX $tickX
        $g.DrawLine($penAxis, [float]$tx, [float]($rcY + $rcH), [float]$tx, [float]($rcY + $rcH + 5))
        $g.DrawString($tickX.ToString(), $fontTicks, $brushText, [float]($tx - 8), [float]($rcY + $rcH + 8))
    }
    for ($tickY = 0; $tickY -le $yMax; $tickY += 5) {
        $ty = & $toY $tickY
        $g.DrawLine($penAxis, [float]$rcX, [float]$ty, [float]($rcX - 5), [float]$ty)
        $g.DrawString($tickY.ToString(), $fontTicks, $brushText, [float]($rcX - 25), [float]($ty - 6))
    }
    
    $g.DrawString("log2 Fold Change", $fontAxis, $brushText, [float]($rcX + $rcW/2 - 50), [float]($rcY + $rcH + 35))
    $g.DrawString("-log10(P-value)", $fontAxis, $brushText, [float]($rcX - 60), [float]($rcY + $rcH/2 - 10))
}

# Legend
$g.FillEllipse($brushUp, 1850, 45, 12, 12); $g.DrawString("Upregulated Matrisome DEG", $fontLegend, $brushText, [float]1870, [float]43)
$g.FillEllipse($brushDown, 1850, 70, 12, 12); $g.DrawString("Downregulated Matrisome DEG", $fontLegend, $brushText, [float]1870, [float]68)
$g.FillEllipse($brushNS, 1850, 95, 10, 10); $g.DrawString("Non-significant / Non-DEG", $fontLegend, $brushText, [float]1870, [float]93)

$volcanoPng = "$projectDir/results/figures/deg/matrisome_volcano_3cohorts.png"
$bmp.Save($volcanoPng, [System.Drawing.Imaging.ImageFormat]::Png)
$bmp.Dispose()
Write-Output "SUCCESS: Saved Volcano Plot to $volcanoPng"


# ==============================================================================
# PART 2: CLUSTERED HEATMAP OF 45 SHARED MATRISOME DEGs
# ==============================================================================
$hmWidth = 1450
$hmHeight = 1500
$bmpHm = New-Object System.Drawing.Bitmap($hmWidth, $hmHeight)
$gHm = [System.Drawing.Graphics]::FromImage($bmpHm)
$gHm.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$gHm.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::AntiAliasGridFit
$gHm.Clear([System.Drawing.Color]::White)

$fontHmTitle = New-Object System.Drawing.Font("Arial", 16, [System.Drawing.FontStyle]::Bold)
$fontHmSub = New-Object System.Drawing.Font("Arial", 11, [System.Drawing.FontStyle]::Regular)
$fontColHdr = New-Object System.Drawing.Font("Arial", 11, [System.Drawing.FontStyle]::Bold)
$fontGene = New-Object System.Drawing.Font("Arial", 9, [System.Drawing.FontStyle]::Bold)
$fontVal = New-Object System.Drawing.Font("Arial", 8, [System.Drawing.FontStyle]::Regular)
$fontCat = New-Object System.Drawing.Font("Arial", 8, [System.Drawing.FontStyle]::Regular)

$shared = Import-Csv "results/tables/common_matrisome_DEGs_preview.csv"
$sorted = $shared | Sort-Object { 
    $m = ([double]$_.logFC_GSE32539 + [double]$_.logFC_GSE110147 + [double]$_.logFC_GSE24206) / 3.0
    return -$m
}

$gHm.DrawString("Expression Heatmap: 45 Shared Matrisome DEGs Across 3 Discovery Cohorts", $fontHmTitle, $brushText, [float]50, [float]25)
$gHm.DrawString("Log2 Fold Change (IPF vs Normal Control) | All 45 genes met |log2FC| > 0.585 and adj.P < 0.05 across all 3 platforms", $fontHmSub, [System.Drawing.Brushes]::DimGray, [float]50, [float]55)

$gridX = 260
$gridY = 120
$cellW = 180
$cellH = 26

$cols = @("logFC_GSE24206", "logFC_GSE32539", "logFC_GSE110147")
$colNames = @("GSE24206 (Affy 1.0)", "GSE32539 (Agilent 4x44K)", "GSE110147 (Affy U133+)")

for ($c = 0; $c -lt 3; $c++) {
    $cx = $gridX + ($c * $cellW)
    $gHm.DrawString($colNames[$c], $fontColHdr, $brushText, [float]($cx + 10), [float]($gridY - 30))
}
$gHm.DrawString("Division", $fontColHdr, $brushText, [float]($gridX + (3 * $cellW) + 20), [float]($gridY - 30))
$gHm.DrawString("Category", $fontColHdr, $brushText, [float]($gridX + (3 * $cellW) + 140), [float]($gridY - 30))

$getColor = {
    param($val)
    if ($val -ge 0) {
        $ratio = [math]::Min(1.0, $val / 4.5)
        $r = [int](255)
        $g = [int](255 - ($ratio * 200))
        $b = [int](255 - ($ratio * 200))
        return [System.Drawing.Color]::FromArgb($r, $g, $b)
    } else {
        $ratio = [math]::Min(1.0, [math]::Abs($val) / 2.5)
        $r = [int](255 - ($ratio * 220))
        $g = [int](255 - ($ratio * 150))
        $b = [int](255)
        return [System.Drawing.Color]::FromArgb($r, $g, $b)
    }
}

$rowIdx = 0
foreach ($item in $sorted) {
    $ry = $gridY + ($rowIdx * $cellH)
    $gn = $item.gene
    
    $suffix = if ($item.direction_consistent -eq "FALSE") { " *" } else { "" }
    $geneColor = if ($item.direction_consistent -eq "FALSE") { [System.Drawing.Brushes]::Crimson } else { [System.Drawing.Brushes]::Black }
    $gHm.DrawString("$gn$suffix", $fontGene, $geneColor, [float]150, [float]($ry + 4))
    
    for ($c = 0; $c -lt 3; $c++) {
        $v = [double]$item.($cols[$c])
        $cx = $gridX + ($c * $cellW)
        $clr = & $getColor $v
        $bCell = New-Object System.Drawing.SolidBrush($clr)
        $gHm.FillRectangle($bCell, [int]$cx, [int]$ry, [int]($cellW - 4), [int]($cellH - 2))
        $gHm.DrawRectangle([System.Drawing.Pens]::LightGray, [int]$cx, [int]$ry, [int]($cellW - 4), [int]($cellH - 2))
        
        $txtColor = if ([math]::Abs($v) -gt 2.0) { [System.Drawing.Brushes]::White } else { [System.Drawing.Brushes]::Black }
        $vStr = if ($v -gt 0) { "+$($v.ToString('F2'))" } else { $v.ToString('F2') }
        $gHm.DrawString($vStr, $fontVal, $txtColor, [float]($cx + 60), [float]($ry + 5))
        $bCell.Dispose()
    }
    
    # Division
    $div = $item.Division
    $divColor = if ($div -like "*Core*") { [System.Drawing.Color]::FromArgb(40, 167, 69) } else { [System.Drawing.Color]::FromArgb(111, 66, 193) }
    $divBrush = New-Object System.Drawing.SolidBrush($divColor)
    $gHm.FillRectangle($divBrush, [int]($gridX + 3 * $cellW + 15), [int]($ry + 3), [int]110, [int]($cellH - 8))
    $divText = if ($div -like "*Core*") { "Core Matrisome" } else { "Associated" }
    $gHm.DrawString($divText, $fontCat, [System.Drawing.Brushes]::White, [float]($gridX + 3 * $cellW + 20), [float]($ry + 5))
    $divBrush.Dispose()
    
    # Category
    $gHm.DrawString($item.Category, $fontCat, [System.Drawing.Brushes]::DimGray, [float]($gridX + 3 * $cellW + 140), [float]($ry + 5))
    
    $rowIdx++
}

# Legend
$legX = 1000
$legY = 1350
$gHm.DrawString("Log2 Fold Change Scale:", $fontColHdr, $brushText, [float]$legX, [float]($legY - 25))
for ($step = -2.5; $step -le 4.5; $step += 0.2) {
    $lx = $legX + (($step + 2.5) / 7.0) * 250
    $cStep = & $getColor $step
    $bStep = New-Object System.Drawing.SolidBrush($cStep)
    $gHm.FillRectangle($bStep, [int]$lx, [int]$legY, [int]8, [int]20)
    $bStep.Dispose()
}
$gHm.DrawRectangle([System.Drawing.Pens]::Black, [int]$legX, [int]$legY, [int]250, [int]20)
$gHm.DrawString("-2.5", $fontTicks, $brushText, [float]($legX - 10), [float]($legY + 25))
$gHm.DrawString("0.0", $fontTicks, $brushText, [float]($legX + 85), [float]($legY + 25))
$gHm.DrawString("+4.5", $fontTicks, $brushText, [float]($legX + 235), [float]($legY + 25))
$gHm.DrawString("* Asterisk indicates discordant fold change direction across cohorts", $fontCat, [System.Drawing.Brushes]::Crimson, [float]50, [float]($hmHeight - 40))

$heatmapPng = "$projectDir/results/figures/deg/matrisome_heatmap_45_shared_degs.png"
$bmpHm.Save($heatmapPng, [System.Drawing.Imaging.ImageFormat]::Png)
$bmpHm.Dispose()
Write-Output "SUCCESS: Saved Heatmap to $heatmapPng"
