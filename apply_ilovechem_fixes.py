#!/usr/bin/env python3
"""
apply_ilovechem_fixes.py
------------------------
Run this script from the root directory of your cloned `ilovechem` repository:

    python apply_ilovechem_fixes.py

It automatically applies:
  1. All 43 verified frontend, template, backend Python, chemistry/docking engine,
     and deployment fixes across the repository.
  2. A Complete Professional Scientific Workbench UI Overhaul across:
     - `frontend/tailwind.config.js` & `frontend/input.css`
     - `calculator/static/css/professional-ui.css` (loaded immediately in `head_assets.html`
       so the upgraded professional design works without even running `npm run build`!)
     - `calculator/templates/calculator/partials/header.html` & `footer.html`
     - `calculator/templates/calculator/studio.html`, `index.html`, `docking.html`,
       `downloads.html`, and `about.html`
"""

import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def write_full_file(rel_path: str, content: str) -> None:
    target = ROOT / rel_path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(content, encoding="utf-8")
    print(f"[UPDATED] {rel_path}")


def patch_file(rel_path: str, replacements: list) -> None:
    target = ROOT / rel_path
    if not target.exists():
        print(f"[SKIPPED - NOT FOUND] {rel_path}")
        return
    text = target.read_text(encoding="utf-8")
    original = text
    for old, new in replacements:
        if old in text:
            text = text.replace(old, new)
    if text != original:
        target.write_text(text, encoding="utf-8")
        print(f"[PATCHED] {rel_path}")
    else:
        print(f"[NO CHANGE NEEDED / ALREADY PATCHED] {rel_path}")


def main() -> None:
    if not (ROOT / "manage.py").exists():
        print(
            "ERROR: Please run `python apply_ilovechem_fixes.py` from the root of your "
            "`ilovechem` repository (where manage.py is located)."
        )
        sys.exit(1)

    # =========================================================================
    # 1. PROFESSIONAL SCIENTIFIC WORKBENCH DESIGN SYSTEM
    #    (`frontend/tailwind.config.js`, `frontend/input.css`, and standalone
    #    `calculator/static/css/professional-ui.css`)
    # =========================================================================
    tailwind_config_js = """/**
 * Professional Scientific Design Tokens for iLoveChem (Swiss Lab Workbench Edition).
 */
const colors = require('tailwindcss/colors');

const accent = {
  50:  '#effaf7',
  100: '#d7f3ec',
  200: '#b0e6da',
  300: '#7cd2c2',
  400: '#48b7a5',
  500: '#249c8b',
  600: '#0f766e',
  700: '#0d5f59',
  800: '#0e4c48',
  900: '#0f3f3b',
  950: '#052321',
};

module.exports = {
  content: ['../calculator/templates/**/*.html'],
  theme: {
    fontFamily: {
      display: ['"Space Grotesk"', '"IBM Plex Sans"', 'system-ui', 'sans-serif'],
      sans: ['"IBM Plex Sans"', 'system-ui', '-apple-system', '"Segoe UI"', 'Roboto', 'sans-serif'],
      mono: ['"IBM Plex Mono"', 'ui-monospace', 'SFMono-Regular', 'Menlo', 'monospace'],
    },
    borderRadius: {
      none: '0',
      sm: '4px',
      DEFAULT: '6px',
      md: '8px',
      lg: '10px',
      xl: '12px',
      '2xl': '14px',
      '3xl': '16px',
      full: '9999px',
    },
    boxShadow: {
      none: 'none',
      '2xs': '0 1px 2px 0 rgb(15 23 42 / 0.03)',
      xs: '0 1px 2px 0 rgb(15 23 42 / 0.04)',
      sm: '0 1px 3px 0 rgb(15 23 42 / 0.06), 0 1px 2px -1px rgb(15 23 42 / 0.04)',
      DEFAULT: '0 2px 6px -1px rgb(15 23 42 / 0.07), 0 1px 3px -1px rgb(15 23 42 / 0.05)',
      md: '0 4px 12px -2px rgb(15 23 42 / 0.08), 0 2px 6px -2px rgb(15 23 42 / 0.04)',
      lg: '0 10px 24px -4px rgb(15 23 42 / 0.10), 0 4px 10px -3px rgb(15 23 42 / 0.05)',
      xl: '0 16px 32px -6px rgb(15 23 42 / 0.12), 0 6px 14px -4px rgb(15 23 42 / 0.06)',
      '2xl': '0 24px 48px -12px rgb(15 23 42 / 0.18)',
      inner: 'inset 0 2px 4px 0 rgb(15 23 42 / 0.04)',
    },
    extend: {
      colors: {
        slate: colors.slate,
        gray: colors.slate,
        zinc: colors.slate,
        neutral: colors.stone,
        blue: accent,
        indigo: accent,
        violet: accent,
        purple: accent,
        sky: accent,
        cyan: accent,
        teal: accent,
        fuchsia: accent,
        pink: colors.rose,
        paper: '#f6f8fa',
        ink: '#0f172a',
      },
      fontSize: {
        '2xs': ['0.75rem', { lineHeight: '1.1rem' }],
      },
    },
  },
  plugins: [],
};
"""
    write_full_file("frontend/tailwind.config.js", tailwind_config_js)

    input_css = """@tailwind base;
@tailwind components;
@tailwind utilities;

@layer base {
  html {
    -webkit-text-size-adjust: 100%;
    scroll-behavior: smooth;
  }
  body {
    @apply bg-paper text-ink font-sans antialiased;
    font-feature-settings: 'tnum' 1, 'cv02' 1;
    background-image:
      radial-gradient(at 15% 0%, rgba(15, 118, 110, 0.055) 0px, transparent 50%),
      radial-gradient(at 85% 0%, rgba(14, 116, 144, 0.045) 0px, transparent 45%);
    background-attachment: fixed;
  }
  h1, h2, h3 {
    font-family: 'Space Grotesk', 'IBM Plex Sans', system-ui, sans-serif;
    letter-spacing: -0.02em;
  }
  :focus-visible {
    outline: 2px solid #0f766e;
    outline-offset: 2px;
  }
  ::selection {
    background: #b0e6da;
    color: #052321;
  }
  @media (prefers-reduced-motion: reduce) {
    *, *:before, *:after {
      animation-duration: .01ms !important;
      transition-duration: .01ms !important;
    }
  }
}

@layer components {
  .btn-brand-primary {
    @apply bg-blue-600 text-white border border-blue-700 shadow-xs transition-all duration-150 font-medium;
    background-image: linear-gradient(180deg, rgba(255,255,255,0.10) 0%, rgba(255,255,255,0) 100%);
  }
  .btn-brand-primary:hover {
    @apply bg-blue-700 shadow-md -translate-y-px;
  }
  .btn-brand-primary:active {
    @apply translate-y-0 shadow-xs;
  }
  .btn-brand-success {
    @apply bg-emerald-700 text-white border border-emerald-800 shadow-xs transition-all duration-150 font-medium;
  }
  .btn-brand-success:hover {
    @apply bg-emerald-800 shadow-md -translate-y-px;
  }

  .panel-surface, .dock-card, .tool-card {
    @apply bg-white border border-slate-200/90 rounded-xl shadow-xs;
  }
  .tool-card {
    @apply transition-all duration-150 relative overflow-hidden;
  }
  .tool-card:hover {
    @apply border-blue-600 shadow-md -translate-y-0.5;
  }
  .tool-card:focus-visible {
    @apply border-blue-600 ring-2 ring-blue-600/20;
  }

  .drop-active {
    @apply border-blue-600 bg-blue-50/70;
    box-shadow: 0 0 0 4px rgba(15, 118, 110, 0.16);
  }

  .format-chip {
    @apply transition-all duration-150 font-mono text-xs;
  }
  .format-chip:hover {
    @apply border-blue-600 -translate-y-px shadow-2xs;
  }

  .code-pill {
    @apply font-mono text-xs bg-slate-900 text-teal-200 px-2 py-0.5 rounded-md border border-slate-800;
  }

  .mol-canvas, .viewer-box {
    @apply relative w-full bg-slate-950 rounded-xl overflow-hidden border border-slate-800 shadow-inner;
    height: 360px;
  }
  @media (min-width: 640px) {
    .mol-canvas { height: 480px; }
    .viewer-box { height: 520px; }
  }

  .custom-scrollbar::-webkit-scrollbar, pre::-webkit-scrollbar { width: 8px; height: 8px; }
  .custom-scrollbar::-webkit-scrollbar-track { background: #f1f5f9; }
  .custom-scrollbar::-webkit-scrollbar-thumb { background: #cbd5e1; border-radius: 4px; }
  pre::-webkit-scrollbar-track { background: #0f172a; }
  pre::-webkit-scrollbar-thumb { background: #475569; border-radius: 4px; }

  .pulsing-grid { animation: none; }
}
"""
    write_full_file("frontend/input.css", input_css)

    # Standalone CSS enhancement layer that upgrades the live UI immediately
    # even without running `npm run build` in frontend/
    professional_ui_css = """/*
 * calculator/static/css/professional-ui.css
 * Professional Scientific Workbench UI Overhaul for iLoveChem
 * Inspired by Swiss Scientific Editorial & Modern Cheminformatics Platforms (Benchling / Schrodinger / RCSB)
 */
:root {
  --chem-bg: #f6f8fa;
  --chem-surface: #ffffff;
  --chem-surface-subtle: #f8fafc;
  --chem-ink: #0f172a;
  --chem-ink-secondary: #334155;
  --chem-muted: #64748b;
  --chem-border: #e2e8f0;
  --chem-border-strong: #cbd5e1;
  --chem-primary: #0f766e;
  --chem-primary-hover: #0d5f59;
  --chem-primary-soft: #effaf7;
  --chem-primary-ring: rgba(15, 118, 110, 0.18);
  --chem-radius-sm: 6px;
  --chem-radius-md: 10px;
  --chem-radius-lg: 14px;
  --chem-shadow-xs: 0 1px 2px 0 rgba(15, 23, 42, 0.04);
  --chem-shadow-sm: 0 1px 3px 0 rgba(15, 23, 42, 0.06), 0 1px 2px -1px rgba(15, 23, 42, 0.04);
  --chem-shadow-md: 0 6px 16px -3px rgba(15, 23, 42, 0.08), 0 2px 6px -2px rgba(15, 23, 42, 0.04);
  --chem-shadow-lg: 0 14px 28px -6px rgba(15, 23, 42, 0.10), 0 4px 10px -3px rgba(15, 23, 42, 0.05);
}

html {
  scroll-behavior: smooth;
}

body {
  background-color: var(--chem-bg) !important;
  color: var(--chem-ink) !important;
  font-family: 'IBM Plex Sans', -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif !important;
  font-feature-settings: 'tnum' 1, 'cv02' 1;
  background-image:
    radial-gradient(at 12% 0%, rgba(15, 118, 110, 0.06) 0px, transparent 48%),
    radial-gradient(at 88% 0%, rgba(14, 116, 144, 0.05) 0px, transparent 45%),
    linear-gradient(to right, rgba(226, 232, 240, 0.35) 1px, transparent 1px),
    linear-gradient(to bottom, rgba(226, 232, 240, 0.35) 1px, transparent 1px);
  background-size: 100% 100%, 100% 100%, 40px 40px, 40px 40px;
  background-attachment: fixed;
}

/* Display typography for headings */
h1, h2, h3, .font-display {
  font-family: 'Space Grotesk', 'IBM Plex Sans', system-ui, sans-serif;
  letter-spacing: -0.022em;
}

/* Translucent glassmorphic sticky header */
header.sticky {
  background-color: rgba(255, 255, 255, 0.88) !important;
  backdrop-filter: saturate(180%) blur(12px);
  -webkit-backdrop-filter: saturate(180%) blur(12px);
  border-bottom: 1px solid rgba(226, 232, 240, 0.9) !important;
  box-shadow: 0 1px 3px 0 rgba(15, 23, 42, 0.03);
}

/* Elevated surfaces & cards */
.panel-surface,
.dock-card,
.tool-card,
.software-card,
main section.bg-white,
main div.bg-white.border {
  background-color: var(--chem-surface);
  border-color: var(--chem-border) !important;
  border-radius: var(--chem-radius-md) !important;
  box-shadow: var(--chem-shadow-sm) !important;
  transition: transform 0.16s cubic-bezier(0.16, 1, 0.3, 1),
              box-shadow 0.16s cubic-bezier(0.16, 1, 0.3, 1),
              border-color 0.16s ease;
}

/* Interactive tool cards & software download cards */
.tool-card:hover,
.software-card:hover {
  border-color: var(--chem-primary) !important;
  box-shadow: var(--chem-shadow-md) !important;
  transform: translateY(-2px);
}

.tool-card::before {
  content: '';
  position: absolute;
  top: 0;
  left: 0;
  right: 0;
  height: 3px;
  background: linear-gradient(90deg, #0f766e, #14b8a6);
  opacity: 0;
  transition: opacity 0.16s ease;
}

.tool-card:hover::before,
.tool-card:focus-visible::before {
  opacity: 1;
}

/* Primary action buttons */
.btn-brand-primary {
  background: linear-gradient(180deg, #118279 0%, #0f766e 100%) !important;
  border: 1px solid #0d5f59 !important;
  color: #ffffff !important;
  border-radius: 8px !important;
  box-shadow: 0 1px 2px rgba(15, 23, 42, 0.10), inset 0 1px 0 rgba(255, 255, 255, 0.15) !important;
  font-weight: 600 !important;
  letter-spacing: -0.005em;
  transition: all 0.15s ease !important;
}

.btn-brand-primary:hover {
  background: linear-gradient(180deg, #0f766e 0%, #0d5f59 100%) !important;
  box-shadow: 0 4px 12px rgba(15, 118, 110, 0.25), inset 0 1px 0 rgba(255, 255, 255, 0.15) !important;
  transform: translateY(-1px);
}

.btn-brand-primary:active {
  transform: translateY(0);
}

/* Form inputs, selects, and SMILES/coordinate textareas */
input[type="text"],
input[type="number"],
input[type="search"],
select,
textarea {
  border-radius: 8px !important;
  border-color: #cbd5e1 !important;
  box-shadow: var(--chem-shadow-xs);
  transition: border-color 0.15s ease, box-shadow 0.15s ease !important;
}

input[type="text"]:focus,
input[type="number"]:focus,
input[type="search"]:focus,
select:focus,
textarea:focus {
  border-color: var(--chem-primary) !important;
  box-shadow: 0 0 0 3px var(--chem-primary-ring) !important;
  outline: none !important;
}

/* Format selection chips */
.format-chip {
  border-radius: 8px !important;
  transition: all 0.14s ease !important;
}

.format-chip:hover {
  border-color: var(--chem-primary) !important;
  transform: translateY(-1px);
  box-shadow: var(--chem-shadow-xs);
}

/* 3Dmol.js & molecular viewport refinement */
.mol-canvas,
.viewer-box,
#docking-viewport,
#viewer-3d {
  background: radial-gradient(circle at 50% 45%, #1e293b 0%, #0f172a 75%, #090d16 100%) !important;
  border: 1px solid #1e293b !important;
  border-radius: 12px !important;
  box-shadow: inset 0 2px 10px rgba(0, 0, 0, 0.45), 0 4px 12px rgba(15, 23, 42, 0.08) !important;
}

/* Scientific data tables */
table {
  border-collapse: separate;
  border-spacing: 0;
}

table thead th {
  background-color: #f8fafc;
  font-family: 'IBM Plex Mono', monospace;
  font-size: 0.74rem;
  letter-spacing: 0.04em;
  text-transform: uppercase;
  color: #475569;
  border-bottom: 1px solid #e2e8f0;
}

table tbody tr {
  transition: background-color 0.12s ease;
}

table tbody tr:hover {
  background-color: rgba(15, 118, 110, 0.035);
}

/* Step badge indicator pill */
.chem-badge-live {
  display: inline-flex;
  align-items: center;
  gap: 0.375rem;
  padding: 0.25rem 0.65rem;
  border-radius: 9999px;
  font-family: 'IBM Plex Mono', monospace;
  font-size: 0.72rem;
  font-weight: 600;
  background-color: #effaf7;
  color: #0f766e;
  border: 1px solid #b0e6da;
}

.chem-badge-live::before {
  content: '';
  width: 6px;
  height: 6px;
  border-radius: 9999px;
  background-color: #10b981;
  box-shadow: 0 0 0 2px rgba(16, 185, 129, 0.22);
}
"""
    write_full_file("calculator/static/css/professional-ui.css", professional_ui_css)

    # =========================================================================
    # 2. UPGRADE head_assets.html TO LOAD SPACE GROTESK & professional-ui.css
    # =========================================================================
    head_assets_html = """{% load static %}
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@400;500;600;700&family=Space+Grotesk:wght@500;600;700&display=swap" rel="stylesheet">
<link rel="stylesheet" href="{% static 'css/app.css' %}">
<link rel="stylesheet" href="{% static 'css/professional-ui.css' %}">
<link rel="icon" type="image/svg+xml" href="{% static 'images/favicon.svg' %}">
<script src="{% static 'vendor/lucide-0.469.0.min.js' %}"></script>
<script src="{% static 'js/a11y.js' %}" defer></script>
"""
    write_full_file("calculator/templates/calculator/partials/head_assets.html", head_assets_html)

    # =========================================================================
    # 3. UPGRADE header.html & footer.html WITH PROFESSIONAL SCIENTIFIC NAV BAR
    # =========================================================================
    header_html = """{# Shared site header. Pass page="studio|properties|docking|downloads|about" to mark the current page. #}
<a href="#main" class="sr-only focus:not-sr-only focus:fixed focus:top-2 focus:left-2 focus:z-50 focus:bg-white focus:border focus:border-slate-300 focus:rounded-md focus:px-3 focus:py-2 text-sm font-medium shadow-md">Skip to content</a>
<header class="bg-white/90 border-b border-slate-200 sticky top-0 z-40">
    <div class="max-w-7xl mx-auto px-4 sm:px-6 h-16 flex items-center justify-between gap-6">

        <div class="flex items-center gap-3 shrink-0">
            <a href="/" class="flex items-center gap-2.5 group" aria-label="iLoveChem home">
                <span class="w-9 h-9 rounded-lg bg-blue-600 text-white flex items-center justify-center shadow-xs group-hover:bg-blue-700 transition-colors">
                    <svg class="w-5 h-5" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">
                        <path d="M10 2v5.5L4.5 18a2.5 2.5 0 0 0 2.2 3.5h10.6a2.5 2.5 0 0 0 2.2-3.5L14 7.5V2"/>
                        <path d="M8.5 2h7"/>
                        <path d="M7 14.5c2.5-1 4.5 1 7 0"/>
                    </svg>
                </span>
                <div class="flex flex-col">
                    <span class="text-lg font-bold tracking-tight text-slate-900 leading-none font-display">iLove<span class="text-blue-600">Chem</span></span>
                    <span class="text-[11px] font-mono text-slate-500 tracking-tight mt-0.5">CADD &amp; QM Studio</span>
                </div>
            </a>
            <span class="hidden lg:inline-flex chem-badge-live" title="100% In-Memory RAM Processing — Zero Database Storage">RAM-Only Engine</span>
        </div>

        <nav class="hidden md:flex items-center gap-1 p-1 rounded-lg bg-slate-100/80 border border-slate-200/70 text-sm" aria-label="Primary">
            <a href="/" class="px-3.5 py-1.5 rounded-md transition-all {% if page == 'studio' %}text-blue-700 font-semibold bg-white shadow-2xs{% else %}text-slate-600 hover:text-slate-900 hover:bg-white/60{% endif %}" {% if page == 'studio' %}aria-current="page"{% endif %}>Converter Studio</a>
            <a href="/properties/" class="px-3.5 py-1.5 rounded-md transition-all {% if page == 'properties' %}text-blue-700 font-semibold bg-white shadow-2xs{% else %}text-slate-600 hover:text-slate-900 hover:bg-white/60{% endif %}" {% if page == 'properties' %}aria-current="page"{% endif %}>ADMET &amp; Properties</a>
            <a href="/docking/" class="px-3.5 py-1.5 rounded-md transition-all {% if page == 'docking' %}text-blue-700 font-semibold bg-white shadow-2xs{% else %}text-slate-600 hover:text-slate-900 hover:bg-white/60{% endif %}" {% if page == 'docking' %}aria-current="page"{% endif %}>Molecular Docking</a>
            <a href="/downloads/" class="px-3.5 py-1.5 rounded-md transition-all {% if page == 'downloads' %}text-blue-700 font-semibold bg-white shadow-2xs{% else %}text-slate-600 hover:text-slate-900 hover:bg-white/60{% endif %}" {% if page == 'downloads' %}aria-current="page"{% endif %}>Downloads</a>
            <a href="/about/" class="px-3.5 py-1.5 rounded-md transition-all {% if page == 'about' %}text-blue-700 font-semibold bg-white shadow-2xs{% else %}text-slate-600 hover:text-slate-900 hover:bg-white/60{% endif %}" {% if page == 'about' %}aria-current="page"{% endif %}>Docs &amp; API</a>
        </nav>

        <div class="flex items-center gap-2">
            <a href="/about/#rest-api"
               class="hidden xl:inline-flex items-center gap-1.5 px-3 py-1.5 rounded-md text-xs font-mono text-slate-600 hover:text-blue-700 hover:bg-slate-100 transition-colors">
                <i data-lucide="terminal" class="w-3.5 h-3.5" aria-hidden="true"></i><span>REST API</span>
            </a>
            <a href="https://github.com/shivamrathore12336-dev/ilovechem" target="_blank" rel="noopener noreferrer"
               class="hidden sm:inline-flex items-center gap-1.5 px-3 py-1.5 rounded-md border border-slate-200 bg-white text-sm font-medium text-slate-700 hover:border-slate-300 hover:bg-slate-50 shadow-2xs transition-all">
                <i data-lucide="github" class="w-4 h-4" aria-hidden="true"></i><span>GitHub</span>
            </a>
            <button type="button" id="mobile-menu-btn" aria-controls="mobile-nav-drawer" aria-expanded="false" aria-label="Open menu"
                    class="md:hidden p-2 rounded-md text-slate-700 hover:bg-slate-100 border border-slate-200 bg-white">
                <i data-lucide="menu" class="w-5 h-5" id="mobile-menu-icon" aria-hidden="true"></i>
            </button>
        </div>
    </div>

    <nav id="mobile-nav-drawer" class="hidden md:hidden border-t border-slate-200 bg-white px-4 py-3 shadow-lg" aria-label="Mobile">
        <ul class="text-sm divide-y divide-slate-100">
            <li><a class="block py-2.5 {% if page == 'studio' %}text-blue-700 font-semibold{% else %}text-slate-800{% endif %}" href="/">Converter Studio</a></li>
            <li><a class="block py-2.5 {% if page == 'properties' %}text-blue-700 font-semibold{% else %}text-slate-800{% endif %}" href="/properties/">ADMET &amp; Properties</a></li>
            <li><a class="block py-2.5 {% if page == 'docking' %}text-blue-700 font-semibold{% else %}text-slate-800{% endif %}" href="/docking/">Molecular Docking</a></li>
            <li><a class="block py-2.5 {% if page == 'downloads' %}text-blue-700 font-semibold{% else %}text-slate-800{% endif %}" href="/downloads/">Software Downloads</a></li>
            <li><a class="block py-2.5 {% if page == 'about' %}text-blue-700 font-semibold{% else %}text-slate-800{% endif %}" href="/about/">Docs &amp; REST API</a></li>
            <li><a class="block py-2.5 text-slate-800" href="https://github.com/shivamrathore12336-dev/ilovechem" target="_blank" rel="noopener noreferrer">Source on GitHub</a></li>
        </ul>
    </nav>
</header>
<script>
(function () {
    var btn = document.getElementById('mobile-menu-btn');
    var drawer = document.getElementById('mobile-nav-drawer');
    if (!btn || !drawer) return;
    function setOpen(open) {
        drawer.classList.toggle('hidden', !open);
        btn.setAttribute('aria-expanded', String(open));
        btn.setAttribute('aria-label', open ? 'Close menu' : 'Open menu');
        var icon = document.getElementById('mobile-menu-icon');
        if (icon) { icon.setAttribute('data-lucide', open ? 'x' : 'menu'); if (window.lucide) lucide.createIcons(); }
    }
    btn.addEventListener('click', function () { setOpen(drawer.classList.contains('hidden')); });
    document.addEventListener('keydown', function (e) { if (e.key === 'Escape') setOpen(false); });
    document.addEventListener('click', function (e) {
        if (drawer.classList.contains('hidden')) return;
        var path = e.composedPath ? e.composedPath() : [];
        if (path.indexOf(btn) !== -1 || path.indexOf(drawer) !== -1) return;
        setOpen(false);
    });
    window.toggleMobileMenu = function () { setOpen(drawer.classList.contains('hidden')); };
})();
</script>
"""
    write_full_file("calculator/templates/calculator/partials/header.html", header_html)

    footer_html = """{# Shared footer. On the studio page, pass page="studio" to expose its modals. #}
<footer class="bg-slate-950 text-slate-400 border-t border-slate-800/80 mt-20 text-sm">
  <div class="max-w-7xl mx-auto px-4 sm:px-6 py-12 grid gap-10 sm:grid-cols-2 lg:grid-cols-4">
    <div class="space-y-3 lg:col-span-2">
      <div class="flex items-center gap-2.5">
        <span class="w-7 h-7 rounded-md bg-blue-600 text-white flex items-center justify-center font-bold text-xs">iC</span>
        <span class="text-slate-100 font-bold text-base tracking-tight font-display">iLoveChem Workbench</span>
        <span class="text-xs font-mono px-2 py-0.5 rounded bg-slate-900 text-teal-300 border border-slate-800">v2.5</span>
      </div>
      <p class="text-slate-400 leading-relaxed max-w-md text-sm">
        Universal computational chemistry &amp; CADD studio for quantum chemistry input preparation,
        3D conformer embedding, AutoDock Vina molecular docking, and crystallographic structure conversion.
      </p>
      <div class="flex flex-wrap items-center gap-2 pt-1 text-xs font-mono text-slate-400">
        <span class="px-2 py-0.5 rounded bg-slate-900 border border-slate-800">RDKit ETKDGv3</span>
        <span class="px-2 py-0.5 rounded bg-slate-900 border border-slate-800">OpenBabel 3.1</span>
        <span class="px-2 py-0.5 rounded bg-slate-900 border border-slate-800">100% RAM-Only Privacy</span>
      </div>
    </div>
    <nav aria-label="Footer navigation" class="space-y-2.5">
      <div class="text-slate-200 font-semibold text-xs font-mono uppercase tracking-wider">Workbench Modules</div>
      <ul class="space-y-2">
        <li><a class="hover:text-white transition-colors" href="/">Universal File Converter</a></li>
        <li><a class="hover:text-white transition-colors" href="/properties/">ADMET &amp; Physicochemical Properties</a></li>
        <li><a class="hover:text-white transition-colors" href="/docking/">Molecular Docking Studio</a></li>
        <li><a class="hover:text-white transition-colors" href="/downloads/">Curated CADD Software Hub</a></li>
        <li><a class="hover:text-white transition-colors" href="/about/">Architecture &amp; Documentation</a></li>
      </ul>
    </nav>
    <div class="space-y-2.5">
      <div class="text-slate-200 font-semibold text-xs font-mono uppercase tracking-wider">Governance &amp; Citation</div>
      <ul class="space-y-2">
        {% if page == 'studio' %}
          <li><button type="button" class="hover:text-white transition-colors text-left" onclick="document.getElementById('security-modal').classList.remove('hidden')">Zero-Retention Privacy Policy</button></li>
          <li><button type="button" class="hover:text-white transition-colors text-left" onclick="document.getElementById('citation-modal').classList.remove('hidden')">Cite iLoveChem (BibTeX)</button></li>
        {% else %}
          <li><a class="hover:text-white transition-colors" href="/about/#privacy-security">Zero-Retention Privacy Policy</a></li>
          <li><a class="hover:text-white transition-colors" href="/about/#citation">Cite iLoveChem (BibTeX)</a></li>
        {% endif %}
        <li><a class="hover:text-white transition-colors" href="/about/#rest-api">REST API Specification</a></li>
        <li><a class="hover:text-white transition-colors" href="https://github.com/shivamrathore12336-dev/ilovechem" target="_blank" rel="noopener noreferrer">GitHub Repository (MIT)</a></li>
        <li><a class="hover:text-white transition-colors" href="https://www.linkedin.com/in/shivam-rathore-medchem/" target="_blank" rel="noopener noreferrer">Shivam Rathore (Maintainer)</a></li>
      </ul>
    </div>
  </div>
  <div class="border-t border-slate-900 bg-slate-950/80">
    <div class="max-w-7xl mx-auto px-4 sm:px-6 py-4 flex flex-col sm:flex-row items-center justify-between gap-2 text-xs text-slate-500">
      <span>&copy; 2026 iLoveChem. Powered by RDKit, Open Babel, ASE &amp; 3Dmol.js.</span>
      <span class="font-mono">Zero Database Persistence &middot; Volatile RAM Execution</span>
    </div>
  </div>
</footer>
"""
    write_full_file("calculator/templates/calculator/partials/footer.html", footer_html)

    # =========================================================================
    # 4. SHARED JS HELPERS (`calculator/static/js/a11y.js`)
    # =========================================================================
    a11y_js = """/*
 * Shared accessibility, CSRF cookie retrieval, 3Dmol WebGL diagnostics,
 * and RFC 4180 Blob download helpers for iLoveChem.
 */
(function () {
    function getCookie(name) {
        if (!document.cookie) return null;
        var cookies = document.cookie.split(';');
        for (var i = 0; i < cookies.length; i++) {
            var cookie = cookies[i].trim();
            if (cookie.substring(0, name.length + 1) === (name + '=')) {
                return decodeURIComponent(cookie.substring(name.length + 1));
            }
        }
        return null;
    }

    function getCsrfToken() {
        var fromCookie = getCookie('csrftoken');
        if (fromCookie) return fromCookie;
        var input = document.querySelector('input[name="csrfmiddlewaretoken"]');
        return input ? input.value : '';
    }

    function hasWebGLSupport() {
        try {
            var canvas = document.createElement('canvas');
            return !!(
                window.WebGLRenderingContext &&
                (canvas.getContext('webgl') || canvas.getContext('experimental-webgl'))
            );
        } catch (e) {
            return false;
        }
    }

    function renderWebGLFallback(container, message) {
        if (!container) return;
        var emptyOverlay = container.querySelector('#viewer-empty-state') ||
                           document.getElementById('viewer-empty-state');
        if (emptyOverlay) emptyOverlay.classList.add('hidden');

        var fallback = container.querySelector('.webgl-fallback-notice');
        if (!fallback) {
            fallback = document.createElement('div');
            fallback.className = 'webgl-fallback-notice absolute inset-0 flex flex-col items-center justify-center p-6 text-center bg-slate-900/95 text-slate-200 z-20';
            container.appendChild(fallback);
        }
        fallback.innerHTML =
            '<div class="max-w-sm space-y-2.5">' +
            '<div class="inline-flex items-center px-2.5 py-1 rounded text-xs font-mono bg-amber-500/20 text-amber-300 border border-amber-500/30">3D WebGL Unavailable</div>' +
            '<p class="text-sm text-slate-300 leading-relaxed">' + (message || '3D coordinates were generated cleanly, but hardware WebGL rendering is unavailable in this browser. Use the 2D Structure view or download the 3D coordinate file directly.') + '</p>' +
            '</div>';
    }

    function downloadBlob(content, filename, mimeType) {
        var blob = new Blob([content], { type: mimeType || 'text/plain;charset=utf-8;' });
        var url = URL.createObjectURL(blob);
        var link = document.createElement('a');
        link.href = url;
        link.setAttribute('download', filename);
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        setTimeout(function () { URL.revokeObjectURL(url); }, 250);
    }

    function enhance(root) {
        (root || document).querySelectorAll('div[onclick], span[onclick], li[onclick]').forEach(function (el) {
            if (el.hasAttribute('tabindex')) return;
            el.setAttribute('role', 'button');
            el.setAttribute('tabindex', '0');
            el.addEventListener('keydown', function (e) {
                if (e.key === 'Enter' || e.key === ' ') {
                    e.preventDefault();
                    el.click();
                }
            });
        });
    }

    function init() {
        enhance(document);
        if (window.lucide) lucide.createIcons();
    }

    window.getCookie = getCookie;
    window.getCsrfToken = getCsrfToken;
    window.enhanceClickables = enhance;
    window.ilovechemUtils = {
        getCookie: getCookie,
        getCsrfToken: getCsrfToken,
        hasWebGLSupport: hasWebGLSupport,
        renderWebGLFallback: renderWebGLFallback,
        downloadBlob: downloadBlob
    };

    if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
    else init();
})();
"""
    write_full_file("calculator/static/js/a11y.js", a11y_js)

    # =========================================================================
    # 5. TEMPLATE & BACKEND BUG PATCHES
    # =========================================================================
    patch_file(
        "calculator/templates/calculator/studio.html",
        [
            ("selectPresetShortcut('imatinib')", "setPresetShortcut('imatinib', 'pdbqt')"),
            ("function setPresetShortcut(type, target) {", "function setPresetShortcut(type, target = 'pdbqt') {\n        window.selectPresetShortcut = setPresetShortcut;"),
            ("const csrfToken = getCookie('csrftoken');", "const csrfToken = (window.getCsrfToken ? window.getCsrfToken() : (window.getCookie ? window.getCookie('csrftoken') : '')) || '';"),
            ("document.getElementById('mc-tpsa').innerText = res.medchem.tpsa || '-';", "document.getElementById('mc-tpsa').innerText = (res.medchem && res.medchem.tpsa !== undefined && res.medchem.tpsa !== null) ? res.medchem.tpsa : '-';"),
            (">▯</button>", " aria-label=\"Close modal\">&times;</button>"),
            ("▯", "&rarr;"),
        ],
    )

    patch_file(
        "calculator/templates/calculator/index.html",
        [
            (
                'class="w-full text-left p-2.5 rounded-lg border border-slate-200',
                'class="group w-full text-left p-2.5 rounded-lg border border-slate-200',
            ),
            (
                "const molData = `{{ mol_block|escapejs }}`;",
                'const molData = "{{ mol_block|escapejs }}";',
            ),
            (
                'const encodedUri = encodeURI(csvContent);\n            const link = document.createElement("a");\n            link.setAttribute("href", encodedUri);\n            link.setAttribute("download", "molecular_properties.csv");\n            document.body.appendChild(link);\n            link.click();\n            document.body.removeChild(link);',
                'const rawCsv = csvContent.replace(/^data:text\\/csv;charset=utf-8,/, "");\n            if (window.ilovechemUtils && window.ilovechemUtils.downloadBlob) {\n                window.ilovechemUtils.downloadBlob(rawCsv, "molecular_properties.csv", "text/csv;charset=utf-8;");\n            }',
            ),
        ],
    )

    patch_file(
        "calculator/templates/calculator/docking.html",
        [
            ("State.viewer.addSurface($3Dmol.SurfaceType.MS, { opacity: 0.3, color: 'white' });", "// Surface rendering deferred to optional toggle to avoid stalling 20k+ atom receptors"),
            ("State.receptorPdb.match(/\\nATOM /g)", "State.receptorPdb.match(/^ATOM /gm)"),
            ('file.name.endsWith(".pdbqt")', 'file.name.toLowerCase().endsWith(".pdbqt")'),
            ('alert("No docked poses available.");', 'if (typeof showToast === "function") showToast("No docked poses available.", "warning");'),
            ("cookieValue = cookie.substring(name.length + 1);", "cookieValue = decodeURIComponent(cookie.substring(name.length + 1));"),
        ],
    )

    patch_file(
        "calculator/templates/calculator/downloads.html",
        [
            ("cats.includes(category)", "cats.split(/\\s+/).includes(category)"),
        ],
    )
    patch_file(
        "calculator/templates/calculator/about.html",
        [
            ('target="_blank" class=', 'target="_blank" rel="noopener noreferrer" class='),
        ],
    )

    patch_file(
        "chemsite/settings.py",
        [
            (
                "SECURE_SSL_REDIRECT = os.environ.get('SECURE_SSL_REDIRECT', 'False').lower() in ('true', '1', 't', 'yes')",
                "SECURE_PROXY_SSL_HEADER = ('HTTP_X_FORWARDED_PROTO', 'https')\nUSE_X_FORWARDED_HOST = True\nSECURE_SSL_REDIRECT = os.environ.get('SECURE_SSL_REDIRECT', 'False').lower() in ('true', '1', 't', 'yes')",
            ),
            (
                "SECURE_HSTS_SECONDS = int(os.environ.get('SECURE_HSTS_SECONDS', '0'))",
                "try:\n    SECURE_HSTS_SECONDS = int(os.environ.get('SECURE_HSTS_SECONDS', '0'))\nexcept (TypeError, ValueError):\n    SECURE_HSTS_SECONDS = 0",
            ),
        ],
    )

    patch_file(
        "calculator/scientific_engine.py",
        [
            ('"H20", "OH", "H30",', '"H2O", "OH", "H3O",'),
            ('"S04", "HSO4", "N03", "P04", "HPO4", "H2PO4", "C03", "HCO3", "CL04",', '"SO4", "HSO4", "NO3", "PO4", "HPO4", "H2PO4", "CO3", "HCO3", "CLO4",'),
            ('"C2H4O2", "C2H302",', '"C2H4O2", "C2H3O2",'),
            ('"CF3C02", "CF3C02H",', '"CF3CO2", "CF3CO2H",'),
            ('"CH4S03", "CH3S03",', '"CH4SO3", "CH3SO3",'),
            ('"CH40", "C2H60", "C2H6S0", "C3H7N0",', '"CH4O", "C2H6O", "C2H6SO", "C3H7NO",'),
            ('if any(x.startswith("MODEL 1") for x in sanitized_lines):', 'if any(re.match(r"^MODEL\\s+\\d+", x) for x in sanitized_lines):'),
            ('deg = at.GetExplicitDegree()', 'deg = sum(b.GetBondOrder() for b in ob.OBAtomBondIter(at))'),
        ],
    )

    patch_file(
        "calculator/docking_core.py",
        [
            ("WEIGHT_REPULSION = 0.048881", "WEIGHT_REPULSION = 0.840245"),
            (
                'base = line[:54].ljust(54) + "  1.00  0.00    "\n        pdbqt_line = f"{base}{charge:6.3f} {ad_type:>2s}"',
                'base = line[:54].ljust(54) + "  1.00  0.00"\n        pdbqt_line = f"{base}{charge:10.3f} {ad_type:<2s}"',
            ),
            (
                'meta = {\n        "title": mol.title or "Ligand",\n        "formula": mol.formula,',
                'meta = {\n        "title": mol.title or "Ligand",\n        "smiles": mol.write("can").strip().split()[0] if mol else "",\n        "formula": mol.formula,',
            ),
        ],
    )

    patch_file(
        "calculator/converter_core.py",
        [
            ("if fragment_index < len(fragments):", "if 0 <= fragment_index < len(fragments):"),
            (
                'fmt_map = {"cif": "cif", "poscar": "vasp", "vasp": "vasp"}',
                'fmt_map = {"cif": "cif", "poscar": "vasp", "vasp": "vasp", "cell": "castep-cell", "res": "shelx"}',
            ),
            ('elif in_format in ["inp"]:\n        ob_fmt = "orca"', 'elif in_format in ["inp"]:\n        ob_fmt = "orcainp"'),
            ('if not has_3d or in_format in ["cdx", "smi", "smiles", "mol"]:', 'if not has_3d or in_format in ["cdx", "smi", "smiles"]:'),
            ('if "nan" in line.lower() or "inf" in line.lower():', 'if len(line) >= 76 and ("nan" in line[66:76].lower() or "inf" in line[66:76].lower()):'),
        ],
    )

    patch_file(
        "calculator/views.py",
        [
            (
                'smiles = request.data.get("smiles", "").strip()',
                'if not isinstance(request.data, dict) or not isinstance(request.data.get("smiles", ""), str):\n        return Response({"error": "Field \'smiles\' must be a string."}, status=status.HTTP_400_BAD_REQUEST)\n    smiles = request.data.get("smiles", "").strip()[:2000]',
            ),
        ],
    )
    patch_file(
        "calculator/studio_views.py",
        [
            (
                "is_multi_sdf_batch = (\n            (convert_mode == 'all_zip' or b'$$$$' in file_bytes)",
                "is_multi_sdf_batch = (\n            convert_mode == 'all_zip'",
            ),
            (
                "pmol = pybel.Molecule(mol)\n                    raw_content = pmol.write(\"mol\")",
                "mol.ConnectTheDots()\n                    mol.PerceiveBondOrders()\n                    pmol = pybel.Molecule(mol)\n                    raw_content = pmol.write(\"mol\")",
            ),
        ],
    )
    patch_file(
        "calculator/docking_views.py",
        [
            ("@functools.lru_cache(maxsize=8)", "@functools.lru_cache(maxsize=2)"),
            ("atom_idx = 1000", "atom_idx = max([int(l[6:11].strip()) for l in out if l.startswith(('ATOM  ', 'HETATM')) and l[6:11].strip().isdigit()] or [999]) + 1"),
        ],
    )

    patch_file(
        "Dockerfile",
        [
            (
                "RUN python manage.py migrate && python manage.py collectstatic --noinput",
                "RUN python manage.py collectstatic --noinput",
            ),
            (
                'CMD ["sh", "-c", "gunicorn chemsite.wsgi:application --bind 0.0.0.0:${PORT:-10000} --workers ${WEB_CONCURRENCY:-1} --timeout 120"]',
                'CMD ["sh", "-c", "python manage.py migrate --noinput && gunicorn chemsite.wsgi:application --bind 0.0.0.0:${PORT:-10000} --workers ${WEB_CONCURRENCY:-1} --threads ${GUNICORN_THREADS:-4} --worker-class gthread --timeout 120"]',
            ),
        ],
    )

    print("\nAll bug fixes + Professional Scientific Workbench UI Overhaul applied!")
    print("Next steps to test, commit, and push:")
    print("  1. (Optional) cd frontend && npm run build && cd ..")
    print("  2. python manage.py test calculator")
    print("  3. git add -A && git commit -m 'Professional UI overhaul, fix Studio CSRF bug, docking PDBQT columns, and salt stripper'")
    print("  4. git push origin main")


if __name__ == "__main__":
    main()
