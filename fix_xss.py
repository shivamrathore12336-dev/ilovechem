#!/usr/bin/env python3
"""Fix remaining XSS and other bugs in templates."""

# Fix docking.html - filter XSS
with open('calculator/templates/calculator/docking.html', 'rb') as f:
    raw = f.read()
content = raw.replace(b'\r\n', b'\n')

old = b"No interactions detected for filter '${filter}'."
new = b"No interactions detected for filter '${escapeHtml(filter)}'."
if old in content:
    content = content.replace(old, new, 1)
    print('FIXED: docking.html filter XSS')
else:
    print('ERROR: docking.html filter XSS pattern not found')

with open('calculator/templates/calculator/docking.html', 'wb') as f:
    f.write(content)

# Fix studio.html - meta.electronic_state XSS (line 2178)
with open('calculator/templates/calculator/studio.html', 'rb') as f:
    raw = f.read()
content = raw.replace(b'\r\n', b'\n')

old = b'>&#x26A1; ${meta.electronic_state}</span>`'
new = b'>&#x26A1; ${escapeHtml(meta.electronic_state)}</span>`'
if old in content:
    content = content.replace(old, new, 1)
    print('FIXED: studio.html electronic_state XSS')
else:
    # Try alternate encoding
    old2 = b'>\u26a1 ${meta.electronic_state}</span>`'
    if old2 in content:
        content = content.replace(old2, new.replace(b'&#x26A1;', b'\u26a1'.encode('utf-8')), 1)
        print('FIXED: studio.html electronic_state XSS (emoji encoding)')
    else:
        idx = content.find(b'electronic_state')
        print(f'studio.html electronic_state found at {idx}')
        print(repr(content[idx-60:idx+80]))

with open('calculator/templates/calculator/studio.html', 'wb') as f:
    f.write(content)

print('Done.')
