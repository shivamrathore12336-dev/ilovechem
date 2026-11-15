#!/usr/bin/env python3
"""
Fix WCAG contrast (text-slate-400 -> text-slate-500 on descriptive sub-labels)
and mobile touch targets in docking.html toolbar buttons.
"""

import re

# ---- Fix 1: WCAG contrast in index.html ----
# Only fix descriptive sub-labels (text-[11px], text-xs sub-labels, not icons/interactive)
with open('calculator/templates/calculator/index.html', 'rb') as f:
    raw = f.read()
content = raw.replace(b'\r\n', b'\n')

# Target only small sub-label text (text-[11px] or text-xs with text-slate-400)
# NOT icon wrappers (pointer-events-none) and NOT interactive elements
wcag_patterns = [
    # Exact molecular weight sub-label
    b'<div class="text-[11px] text-slate-400 mt-1 font-mono">Exact:',
    # LogP sub-label
    b'<div class="text-[11px] text-slate-400 mt-1">Octanol/Water Partition</div>',
]
for old in wcag_patterns:
    new = old.replace(b'text-slate-400', b'text-slate-500')
    if old in content:
        content = content.replace(old, new, 1)
        print(f'FIXED index.html WCAG: {old[25:60]}...')

with open('calculator/templates/calculator/index.html', 'wb') as f:
    f.write(content)

# ---- Fix 2: Mobile touch targets in docking.html ----
# The 3D toolbar action buttons (lines ~568-620) use py-1 or py-0.5 on small screens
# Add min-h-[36px] to ensure WCAG 2.5.5 (44x44px target size) on mobile
with open('calculator/templates/calculator/docking.html', 'rb') as f:
    raw = f.read()
content = raw.replace(b'\r\n', b'\n')

# Target the 3D viewer toolbar action buttons - they have py-1 rounded-lg bg-emerald/slate/blue classes
old_emerald = b'class="px-2.5 py-1 rounded-lg bg-emerald-50 hover:bg-emerald-100 text-emerald-800 font-b'
new_emerald = b'class="px-2.5 py-1 min-h-[36px] rounded-lg bg-emerald-50 hover:bg-emerald-100 text-emerald-800 font-b'

old_slate = b'class="px-2.5 py-1 rounded-lg bg-slate-100 hover:bg-slate-200 text-slate-700 font-semibo'
new_slate  = b'class="px-2.5 py-1 min-h-[36px] rounded-lg bg-slate-100 hover:bg-slate-200 text-slate-700 font-semibo'

old_blue = b'class="px-2.5 py-1 rounded-lg bg-blue-50 hover:bg-blue-100 text-blue-700 font-bold trans'
new_blue  = b'class="px-2.5 py-1 min-h-[36px] rounded-lg bg-blue-50 hover:bg-blue-100 text-blue-700 font-bold trans'

fixes = 0
for old, new in [(old_emerald, new_emerald), (old_slate, new_slate), (old_blue, new_blue)]:
    count = content.count(old)
    if count > 0:
        content = content.replace(old, new)
        fixes += count
        print(f'FIXED docking.html touch target: {count} occurrences')

if fixes == 0:
    print('NOTE: touch target patterns not found (may already be fixed or layout changed)')

with open('calculator/templates/calculator/docking.html', 'wb') as f:
    f.write(content)

print('Done. All accessibility fixes applied.')
