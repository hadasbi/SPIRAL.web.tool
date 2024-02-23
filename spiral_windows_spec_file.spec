# -*- mode: python ; coding: utf-8 -*-

b = [
    ('C:/Users/hadas.biran/Documents/GitHub/SPIRAL.web.tool/spiral_venv/Lib/site-packages/tables/libblosc2.dll', './tables')
    ]

a = Analysis(
    ['spiral.py'],
    pathex=[],
    binaries=b,
    datas=[ ('ensembl/caenorhabditis_elegans_gene_annotation.tsv', './ensembl'),  ('ensembl/danio_rerio_gene_annotation.tsv', './ensembl'), ('ensembl/drosophila_melanogaster_gene_annotation.tsv', './ensembl'), ('ensembl/homo_sapiens_gene_annotation.tsv', './ensembl'), ('ensembl/mus_musculus_gene_annotation.tsv', './ensembl'), ('ensembl/rattus_norvegicus_gene_annotation.tsv', './ensembl'), ('ensembl/saccharomyces_cerevisiae_gene_annotation.tsv', './ensembl')],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='spiral',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='spiral',
)
