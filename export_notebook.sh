#!/bin/bash
# Quick export script for sharing notebook results

echo "🚀 Exporting notebook with all outputs..."

# Export to HTML (most compatible)
echo "📄 Creating HTML version..."
$HOME/miniconda3/bin/conda run -n md-dock jupyter nbconvert \
    --to html \
    --execute \
    visualize_md_docking.ipynb \
    --output visualize_md_docking.html

# Export to PDF (for printing/reports)
echo "📄 Creating PDF version..."
$HOME/miniconda3/bin/conda run -n md-dock jupyter nbconvert \
    --to pdf \
    --execute \
    visualize_md_docking.ipynb \
    --output visualize_md_docking.pdf \
    2>/dev/null || echo "⚠️  PDF export failed (LaTeX not installed)"

echo ""
echo "✅ Export complete!"
echo ""
echo "Files created:"
ls -lh visualize_md_docking.html visualize_md_docking.pdf 2>/dev/null
echo ""
echo "📤 Ready to share:"
echo "  • visualize_md_docking.html (opens in any browser)"
echo "  • visualize_md_docking.pdf (universal PDF format)"
echo ""
echo "💡 Tip: Right-click file → 'Share' or attach to email"
