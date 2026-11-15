#!/usr/bin/env python3
"""Fix studio.html electronic_state XSS."""

with open('calculator/templates/calculator/studio.html', 'rb') as f:
    raw = f.read()
content = raw.replace(b'\r\n', b'\n')

# The ⚡ emoji is UTF-8: \xe2\x9a\xa1
old_line = b'\xe2\x9a\xa1 ${meta.electronic_state}</span>`'
new_line = b'\xe2\x9a\xa1 ${escapeHtml(meta.electronic_state)}</span>`'

if old_line in content:
    content = content.replace(old_line, new_line, 1)
    print('FIXED: studio.html electronic_state XSS')
else:
    print('NOT FOUND - check raw bytes')
    idx = content.find(b'meta.electronic_state}</span>')
    print(repr(content[idx-5:idx+50]))

with open('calculator/templates/calculator/studio.html', 'wb') as f:
    f.write(content)
print('Done.')
