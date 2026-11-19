from django.http import JsonResponse
from rest_framework import status
from rest_framework.decorators import api_view
from rest_framework.response import Response
from calculator.medchem_engine import analyze_medchem_profile, batch_screen_smiles
from calculator.qm_linter import lint_qm_geometry


def healthz_view(request):
    """Ultra-lightweight health check for Render (no DB or template overhead)."""
    return JsonResponse({"status": "ok", "service": "ilovechem", "storage": "ram-only"})


@api_view(["POST"])
def api_medchem_screen(request):
    """Single or batch SMILES screening with PAINS/Brenk/NIH alerts & BOILED-Egg."""
    if not isinstance(request.data, dict):
        return Response({"error": "JSON object required."}, status=status.HTTP_400_BAD_REQUEST)

    if "smiles_list" in request.data:
        raw_list = request.data.get("smiles_list")
        if isinstance(raw_list, str):
            lines = raw_list.splitlines()
        elif isinstance(raw_list, list):
            lines = [str(x) for x in raw_list]
        else:
            return Response({"error": "'smiles_list' must be a newline string or list of SMILES."}, status=status.HTTP_400_BAD_REQUEST)
        return Response(batch_screen_smiles(lines), status=status.HTTP_200_OK)

    smiles = request.data.get("smiles", "")
    if not isinstance(smiles, str) or not smiles.strip():
        return Response({"error": "Provide a valid 'smiles' string or 'smiles_list'."}, status=status.HTTP_400_BAD_REQUEST)
    profile = analyze_medchem_profile(smiles)
    if profile is None:
        return Response({"error": "Could not parse SMILES structure."}, status=status.HTTP_400_BAD_REQUEST)
    return Response(profile, status=status.HTTP_200_OK)


@api_view(["POST"])
def api_qm_lint(request):
    """Run HPC Pre-Flight Linter & generate a SLURM submission script."""
    if not isinstance(request.data, dict):
        return Response({"error": "JSON object required."}, status=status.HTTP_400_BAD_REQUEST)
    atoms = request.data.get("atoms", [])
    if not isinstance(atoms, list):
        return Response({"error": "'atoms' must be a list of [symbol, x, y, z]."}, status=status.HTTP_400_BAD_REQUEST)
    try:
        parsed_atoms = [(str(a[0]), float(a[1]), float(a[2]), float(a[3])) for a in atoms]
        charge = int(request.data.get("charge", 0))
        mult = int(request.data.get("multiplicity", 1))
        nprocs = int(request.data.get("nprocs", 8))
    except (ValueError, TypeError, IndexError):
        return Response({"error": "Malformed atom coordinates or integer parameters."}, status=status.HTTP_400_BAD_REQUEST)

    engine = str(request.data.get("engine", "orca")).strip().lower()
    return Response(lint_qm_geometry(parsed_atoms, charge=charge, multiplicity=mult, engine=engine, nprocs=nprocs))
