import os
import re
import base64
from django.conf import settings
from django.shortcuts import render
from django.http import HttpResponse
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import status

from .converter_core import (
    detect_format,
    parse_and_extract_atoms,
    export_to_target_format,
    _apply_integrity_audit,
    COMPLETE_TARGET_FORMATS
)
from .phase1_engine import sanitize_and_extract_all_fragments, create_zip_archive_of_all_fragments
from .scientific_engine import QuantumLogParser, DeepContentSniffer, BatchSdfProcessor
from openbabel import openbabel as ob
from openbabel import pybel

def studio_view(request):
    return render(request, "calculator/studio.html", {
        "target_registry": COMPLETE_TARGET_FORMATS
    })


def downloads_view(request):
    return render(request, "calculator/downloads.html")


def about_view(request):
    return render(request, "calculator/about.html")



def safe_int(val, default):
    try:
        if val is None or val == '':
            return default
        return int(val)
    except (ValueError, TypeError):
        return default

@api_view(['POST'])
def api_convert_file(request):
    uploaded_file = request.FILES.get('file')
    raw_content = request.data.get('raw_content', '')
    target_format = str(request.data.get('target_format', 'orca') or 'orca').lower()
    supported_formats = {fmt for group in COMPLETE_TARGET_FORMATS.values() for fmt, _, _ in group}
    supported_formats.update({"inp", "com", "gaussian", "in", "vina", "gromacs"})
    if target_format not in supported_formats:
        return Response({"error": "Unsupported target format."}, status=status.HTTP_400_BAD_REQUEST)
    in_format = str(request.data.get('input_format', '') or '').lower()
    fragment_index = safe_int(request.data.get('fragment_index'), 0)
    convert_mode = str(request.data.get('convert_mode', 'single') or 'single') # 'single' or 'all_zip'
    pdb_filter_mode = str(request.data.get('pdb_filter_mode', 'all') or 'all')
    add_hydrogens = str(request.data.get('add_hydrogens', 'true')).lower() in ['true', '1']
    strip_salts = str(request.data.get('strip_salts', 'true')).lower() in ['true', '1']

    raw_nprocs = safe_int(request.data.get('nprocs'), 8)
    nprocs = max(1, min(128, raw_nprocs if raw_nprocs > 0 else 8))
    raw_mem = safe_int(request.data.get('mem_gb'), 16)
    mem_gb = max(1, min(1024, raw_mem if raw_mem > 0 else 16))

    qm_options = {
        "method": str(request.data.get('method', 'B3LYP') or 'B3LYP'),
        "basis": str(request.data.get('basis', 'def2-SVP') or 'def2-SVP'),
        "charge": safe_int(request.data.get('charge'), 0),
        "mult": safe_int(request.data.get('mult'), 1),
        "calc_type": str(request.data.get('calc_type', 'opt freq') or 'opt freq'),
        "nprocs": nprocs,
        "mem_gb": mem_gb
    }

    filename = "input.xyz"
    if uploaded_file:
        raw_name = os.path.basename(uploaded_file.name)
        # Security Guard: Path traversal & special character filter
        filename = re.sub(r'[^\w\.\-]', '_', raw_name)
        if not in_format:
            in_format = detect_format(filename, "")

        # Giant 50 MB - 500 MB Log Streaming Parser:
        # Avoid loading giant logs in full. Seek directly to the final 4MB and extract the latest geometry.
        if in_format in ["log", "out", "output", "pwout"] and uploaded_file.size > 2 * 1024 * 1024:
            scan_res = QuantumLogParser.reverse_stream_extract_from_file(uploaded_file)
            if scan_res:
                rev_coords, parser_label = scan_res
                ob_mol = ob.OBMol()
                for sym, x, y, z in rev_coords:
                    at = ob_mol.NewAtom()
                    at.SetAtomicNum(ob.GetAtomicNum(sym))
                    at.SetVector(x, y, z)
                mol = pybel.Molecule(ob_mol)
                meta = {
                    "is_2d": False,
                    "is_periodic": False,
                    "parsed_by": f"Streaming Tail Parser ({parser_label})",
                    "formula": mol.formula,
                    "num_atoms": len(mol.atoms),
                    "input_format": in_format,
                    "warnings": [
                        f"Giant Log Tail Parser Active: Scanned the final segment of a {uploaded_file.size / (1024*1024):.1f} MB file and extracted the latest stationary point."
                    ]
                }
                mol = _apply_integrity_audit(mol, meta, strip_salts=strip_salts)
                user_charge = request.data.get('charge')
                user_mult = request.data.get('mult')
                if user_charge in [None, '']:
                    qm_options['charge'] = meta.get('suggested_charge', meta.get('charge', 0))
                if user_mult in [None, '']:
                    qm_options['mult'] = meta.get('suggested_mult', 1)
                out_content, out_filename, export_warnings = export_to_target_format(mol, target_format, qm_options)
                if export_warnings:
                    meta["warnings"].extend(export_warnings)
                viewer_molblock = mol.write("mol")
                return Response({
                    "success": True,
                    "is_zip": False,
                    "filename": out_filename,
                    "output_content": out_content,
                    "detected_input_format": in_format,
                    "target_format": target_format,
                    "metadata": meta,
                    "mol_block_3d": viewer_molblock
                }, status=status.HTTP_200_OK)

        # The giant-log tail scanner above moves the file pointer; rewind before the full read.
        uploaded_file.seek(0)
        file_bytes = uploaded_file.read()
        if not file_bytes.strip():
            return Response({"error": "Empty chemical file provided."}, status=status.HTTP_400_BAD_REQUEST)
    elif raw_content and raw_content.strip():
        if len(raw_content) > 25 * 1024 * 1024:
            return Response({"error": "Raw text input exceeds maximum allowed size (25 MB)."}, status=status.HTTP_400_BAD_REQUEST)
        file_bytes = raw_content.encode('utf-8')
    else:
        return Response({"error": "No file or text provided."}, status=status.HTTP_400_BAD_REQUEST)

    # Deep Content Sniffing & Corruption Guard (Never trust file extension blindly)
    sniffed_fmt, diag_error = DeepContentSniffer.sniff_format(file_bytes, filename)
    if diag_error:
        return Response({"error": diag_error}, status=status.HTTP_400_BAD_REQUEST)
    
    if sniffed_fmt:
        # Auto-correct if file extension was mismatched or missing
        in_format = sniffed_fmt
    elif not in_format:
        in_format = "xyz"

    # Special Case: User wants ALL molecules in reaction scheme bundled into a ZIP file
    if convert_mode == 'all_zip' and (in_format == 'cdxml' or b'CDXML' in file_bytes[:1000]):
        try:
            cdxml_str = file_bytes.decode('utf-8', errors='ignore')
            fragments, stats = sanitize_and_extract_all_fragments(cdxml_str)
            if len(fragments) > 0:
                zip_bytes, zip_filename = create_zip_archive_of_all_fragments(
                    fragments, target_format, export_to_target_format, qm_options
                )
                b64_zip = base64.b64encode(zip_bytes).decode('utf-8')
                return Response({
                    "success": True,
                    "is_zip": True,
                    "filename": zip_filename,
                    "zip_base64": b64_zip,
                    "output_content": f"# Reaction Scheme Multi-Molecule Export\n# {len(fragments)} chemical structures exported into ZIP archive.\n# Target Format: {target_format.upper()}",
                    "metadata": {
                        "formula": f"{len(fragments)} Molecules",
                        "num_atoms": sum(f['atom_count'] for f in fragments),
                        "parsed_by": "Reaction Scheme ZIP Generator"
                    }
                }, status=status.HTTP_200_OK)
        except Exception:
            return Response({"error": "Failed to generate the reaction-scheme archive."}, status=status.HTTP_422_UNPROCESSABLE_ENTITY)

    # Special Case: Multi-Structure SDF batch export to ZIP (Loophole 5)
    if (convert_mode == 'all_zip' or b'$$$$' in file_bytes) and file_bytes.count(b'$$$$') >= 2:
        try:
            sdf_str = file_bytes.decode('utf-8', errors='ignore')
            records = BatchSdfProcessor.split_sdf_records(sdf_str)
            if len(records) > 1:
                zip_bytes, zip_filename, count = BatchSdfProcessor.create_batch_zip(
                    records, target_format, export_to_target_format, qm_options
                )
                b64_zip = base64.b64encode(zip_bytes).decode('utf-8')
                return Response({
                    "success": True,
                    "is_zip": True,
                    "filename": zip_filename,
                    "zip_base64": b64_zip,
                    "output_content": f"# Batch SDF Virtual Screening Library Export\n# {count} chemical structures exported into ZIP archive.\n# Target Format: {target_format.upper()}",
                    "metadata": {
                        "formula": f"{count} Molecules",
                        "num_atoms": count,
                        "parsed_by": "Batch SDF Virtual Screening ZIP Generator"
                    }
                }, status=status.HTTP_200_OK)
        except Exception:
            return Response({"error": "Failed to generate the SDF batch archive."}, status=status.HTTP_422_UNPROCESSABLE_ENTITY)

    # Standard Single / Selected Molecule Conversion
    mol, meta = parse_and_extract_atoms(
        file_bytes, in_format, filename,
        fragment_index=fragment_index,
        pdb_filter_mode=pdb_filter_mode,
        add_hydrogens=add_hydrogens,
        strip_salts=strip_salts
    )
    if mol is None:
        return Response({
            "error": f"Invalid Chemical Structure / Corrupted File: Could not extract valid chemical coordinates as '{in_format.upper()}'. Please check file syntax or ensure the file is not empty or corrupted."
        }, status=status.HTTP_422_UNPROCESSABLE_ENTITY)

    # Smart Spin & Charge: Default to detected electronic state if not explicitly specified
    user_charge = request.data.get('charge')
    user_mult = request.data.get('mult')
    if user_charge is None or user_charge == '':
        qm_options['charge'] = meta.get('detected_charge', 0)
    else:
        try: qm_options['charge'] = int(user_charge)
        except ValueError: qm_options['charge'] = meta.get('detected_charge', 0)

    if user_mult is None or user_mult == '':
        qm_options['mult'] = meta.get('suggested_mult', 1)
    else:
        try: qm_options['mult'] = int(user_mult)
        except ValueError: qm_options['mult'] = meta.get('suggested_mult', 1)

    qm_options['ghost_map'] = meta.get('ghost_map', {})

    try:
        out_content, out_filename, export_warnings = export_to_target_format(mol, target_format, qm_options)
        if export_warnings:
            meta["warnings"].extend(export_warnings)

        viewer_molblock = mol.write("mol")

        # In-Memory Real-Time MedChem & Lipinski Ro5 Descriptors
        medchem_props = None
        try:
            from rdkit import Chem
            from rdkit.Chem import Descriptors, Lipinski, Crippen
            rd_mol = Chem.MolFromMolBlock(viewer_molblock, sanitize=False)
            if rd_mol is not None:
                try:
                    Chem.SanitizeMol(rd_mol)
                except Exception:
                    pass
                mw = round(float(Descriptors.MolWt(rd_mol)), 2)
                logp = round(float(Crippen.MolLogP(rd_mol)), 2)
                hbd = int(Lipinski.NumHDonors(rd_mol))
                hba = int(Lipinski.NumHAcceptors(rd_mol))
                tpsa = round(float(Descriptors.TPSA(rd_mol)), 2)
                rotb = int(Lipinski.NumRotatableBonds(rd_mol))

                violations = 0
                if mw > 500: violations += 1
                if logp > 5.0: violations += 1
                if hbd > 5: violations += 1
                if hba > 10: violations += 1

                medchem_props = {
                    "mw": mw,
                    "logp": logp,
                    "hbd": hbd,
                    "hba": hba,
                    "tpsa": tpsa,
                    "rotb": rotb,
                    "lipinski_violations": violations,
                    "is_druglike": (violations <= 1)
                }
        except Exception:
            pass
        
        return Response({
            "success": True,
            "is_zip": False,
            "filename": out_filename,
            "output_content": out_content,
            "detected_input_format": in_format,
            "target_format": target_format,
            "metadata": meta,
            "medchem": medchem_props,
            "mol_block_3d": viewer_molblock
        }, status=status.HTTP_200_OK)

    except Exception:
        return Response({"error": f"Conversion failed for target format '{target_format}'."}, status=status.HTTP_422_UNPROCESSABLE_ENTITY)


def _site_url() -> str:
    """Canonical public origin. Fixed by configuration (not the request Host header)
    so robots.txt / sitemap.xml are always https and cannot be poisoned by a forged Host."""
    return settings.SITE_URL.rstrip("/")


def robots_txt_view(request):
    lines = [
        "User-agent: *",
        "Allow: /",
        "Disallow: /admin/",
        f"Sitemap: {_site_url()}/sitemap.xml"
    ]
    return HttpResponse("\n".join(lines), content_type="text/plain; charset=utf-8")


def sitemap_xml_view(request):
    base = _site_url()
    pages = [
        ("/", "weekly", "1.0"),
        ("/docking/", "weekly", "1.0"),
        ("/properties/", "monthly", "0.8"),
        ("/downloads/", "weekly", "0.8"),
        ("/about/", "monthly", "0.8"),
    ]
    entries = "\n".join(
        f"  <url>\n    <loc>{base}{path}</loc>\n    <changefreq>{freq}</changefreq>\n    <priority>{prio}</priority>\n  </url>"
        for path, freq, prio in pages
    )
    xml = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">\n'
        f"{entries}\n"
        "</urlset>"
    )
    return HttpResponse(xml, content_type="application/xml; charset=utf-8")


def og_image_view(request):
    """Serve the 1200x630 OpenGraph / Twitter sharing banner card."""
    img_path = os.path.join(settings.BASE_DIR, "calculator", "static", "images", "og-image.png")
    if os.path.exists(img_path):
        with open(img_path, "rb") as f:
            resp = HttpResponse(f.read(), content_type="image/png")
            resp["Cache-Control"] = "public, max-age=86400"
            return resp
    return HttpResponse(status=404)


def google_verification_view(request):
    return HttpResponse("google-site-verification: google889b369925c8069a.html", content_type="text/html; charset=utf-8")




