from django.shortcuts import render, redirect
from django.contrib.auth import login
from django.contrib.auth.models import User
from django.views.decorators.http import require_POST
from django.contrib import messages
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import status
from .chem_utils import (
    validate_smiles,
    calculate_properties,
    generate_3d_molblock,
    generate_2d_base64_image
)
from .medchem_engine import analyze_medchem_profile
from .lead_optimizer import generate_bioisostere_analogs

MAX_SMILES_LEN = 2000

def home_view(request):
    """
    Main interactive dashboard view.
    Handles form submission, computes properties and 3D structure,
    and returns context for instant rendering.
    """
    context = {}
    default_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"  # Aspirin
    selected_smiles = default_smiles

    if request.method == "POST":
        raw = request.POST.get("smiles", "")
        selected_smiles = str(raw).strip() if raw is not None else ""

    if selected_smiles:
        context["input_smiles"] = selected_smiles
        if len(selected_smiles) > MAX_SMILES_LEN:
            context["error"] = f"SMILES input exceeds maximum allowed length ({MAX_SMILES_LEN} characters)."
        elif not validate_smiles(selected_smiles):
            context["error"] = f"Invalid SMILES string: '{selected_smiles}'. Please verify the chemical notation."
        else:
            props = calculate_properties(selected_smiles)
            if "error" in props:
                context["error"] = props["error"]
            else:
                context["properties"] = props
                context["img_2d"] = generate_2d_base64_image(selected_smiles)
                structure_res = generate_3d_molblock(selected_smiles)
                if "mol_block" in structure_res:
                    context["mol_block"] = structure_res["mol_block"]
                else:
                    context["structure_warning"] = structure_res.get("error", "3D generation unavailable")

                medchem = analyze_medchem_profile(selected_smiles)
                if medchem:
                    context["medchem"] = medchem

                lead_opt = generate_bioisostere_analogs(selected_smiles)
                if lead_opt and lead_opt.get("analogs"):
                    context["bioisosteres"] = lead_opt["analogs"]

                # Pharma Zero-Storage Invariant: Never persist user SMILES to disk or SQLite database.

    context["recent_queries"] = []
    return render(request, "calculator/index.html", context)


@api_view(['POST'])
def api_calculate_properties(request):
    """
    JSON API Endpoint for external tools and programmatic use.
    POST Body: {"smiles": "..."}
    """
    if not isinstance(request.data, dict):
        return Response({"error": "Invalid JSON object payload."}, status=status.HTTP_400_BAD_REQUEST)

    raw_smiles = request.data.get("smiles", "")
    if not isinstance(raw_smiles, str):
        return Response({"error": "Field 'smiles' must be a string."}, status=status.HTTP_400_BAD_REQUEST)

    smiles = raw_smiles.strip()
    if not smiles:
        return Response(
            {"error": "Please provide a non-empty 'smiles' parameter."},
            status=status.HTTP_400_BAD_REQUEST
        )
    if len(smiles) > MAX_SMILES_LEN:
        return Response(
            {"error": f"SMILES input exceeds maximum length ({MAX_SMILES_LEN} chars)."},
            status=status.HTTP_400_BAD_REQUEST
        )

    if not validate_smiles(smiles):
        return Response(
            {"error": f"Invalid SMILES string: '{smiles}'. Could not parse chemical structure."},
            status=status.HTTP_400_BAD_REQUEST
        )

    props = calculate_properties(smiles)
    if "error" in props:
        return Response(props, status=status.HTTP_400_BAD_REQUEST)

    structure_3d = generate_3d_molblock(smiles)
    mol_block = structure_3d.get("mol_block", "")
    warn_3d = structure_3d.get("error", None)

    img_2d = generate_2d_base64_image(smiles)

    resp_data = {
        "success": True,
        "properties": props,
        "mol_block": mol_block,
        "image_2d": img_2d
    }
    medchem = analyze_medchem_profile(smiles)
    if medchem:
        resp_data["medchem"] = medchem
    if warn_3d:
        resp_data["warning_3d"] = warn_3d

    return Response(resp_data, status=status.HTTP_200_OK)


@require_POST
def guest_login_view(request):
    """
    1-Click Guest Researcher login.
    Enables users to immediately test authenticated features (saved docking boxes,
    lab sessions, API access) without requiring Google/GitHub OAuth setup.
    """
    user, _ = User.objects.get_or_create(
        username='guest_researcher',
        defaults={
            'email': 'guest.researcher@ilovechem.org',
            'first_name': 'Guest',
            'last_name': 'Researcher',
        }
    )
    login(request, user, backend='django.contrib.auth.backends.ModelBackend')
    messages.success(request, "Welcome! You are signed in as Guest Researcher.")
    next_url = request.POST.get('next')
    if not next_url or not next_url.startswith('/'):
        next_url = '/'
    return redirect(next_url)

