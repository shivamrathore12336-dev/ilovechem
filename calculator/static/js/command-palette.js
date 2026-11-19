/* calculator/static/js/command-palette.js */
(function () {
    const ACTIONS = [
        { title: 'Universal File Converter Studio', sub: '24 formats \u00b7 ORCA, Gaussian, PDBQT, CIF, SDF', url: '/' },
        { title: 'ADMET, PAINS & BOILED-Egg Properties', sub: 'Lipinski Ro5, QED, TPSA, structural alerts', url: '/properties/' },
        { title: 'Molecular Docking Studio (AutoDock Vina)', sub: 'RCSB PDB fetch, grid box, 2D interaction map', url: '/docking/' },
        { title: 'Curated CADD Software Downloads Hub', sub: 'PyMOL, ORCA 6, Avogadro 2, ChimeraX, Vina', url: '/downloads/' },
        { title: 'REST API & Zero-Retention Architecture Docs', sub: 'Endpoints, BibTeX citation, privacy specs', url: '/about/#rest-api' },
        { title: 'Preset: Imatinib \u2192 AutoDock PDBQT', sub: 'Load kinase inhibitor sample in Converter', action: function() { window.location.href = '/?preset=imatinib'; } },
        { title: 'Preset: COX-2 (1CX2) + Aspirin Docking', sub: 'Jump straight to Docking Studio', url: '/docking/' }
    ];

    var modal = null;
    var input = null;
    var list = null;

    function buildModal() {
        modal = document.createElement('div');
        modal.id = 'ilovechem-cmd-palette';
        modal.className = 'hidden fixed inset-0 z-50 bg-slate-950/60 items-start justify-center pt-20 px-4';
        modal.style.display = 'none';
        modal.setAttribute('role', 'dialog');
        modal.setAttribute('aria-modal', 'true');
        modal.setAttribute('aria-label', 'Quick command palette');
        modal.innerHTML =
            '<div class="w-full max-w-lg bg-white dark:bg-slate-800 rounded-xl border border-slate-200 dark:border-slate-700 shadow-2xl overflow-hidden">' +
            '  <div class="flex items-center gap-2 px-4 py-3 border-b border-slate-200 dark:border-slate-700">' +
            '    <span class="font-mono text-xs px-1.5 py-0.5 rounded bg-slate-100 dark:bg-slate-700 text-slate-600 dark:text-slate-300">\u2318K</span>' +
            '    <input id="cmd-palette-input" type="text" placeholder="Jump to tool, format, or preset..." ' +
            '           class="w-full text-sm text-slate-900 dark:text-slate-100 bg-transparent border-0 focus:outline-none focus:ring-0" />' +
            '    <button type="button" id="cmd-palette-close" class="text-xs font-mono text-slate-500 hover:text-slate-800 dark:hover:text-slate-200">ESC</button>' +
            '  </div>' +
            '  <ul id="cmd-palette-list" class="max-h-72 overflow-y-auto divide-y divide-slate-100 dark:divide-slate-700"></ul>' +
            '</div>';
        document.body.appendChild(modal);

        input = document.getElementById('cmd-palette-input');
        list = document.getElementById('cmd-palette-list');

        document.getElementById('cmd-palette-close').addEventListener('click', closePalette);
        modal.addEventListener('click', function (e) { if (e.target === modal) closePalette(); });
        input.addEventListener('input', function () { renderItems(input.value); });
        renderItems('');
    }

    function renderItems(query) {
        var q = (query || '').trim().toLowerCase();
        var matches = ACTIONS.filter(function (item) {
            return !q || item.title.toLowerCase().includes(q) || item.sub.toLowerCase().includes(q);
        });
        list.innerHTML = '';
        if (!matches.length) {
            list.innerHTML = '<li class="px-4 py-6 text-center text-sm text-slate-500">No matching chemistry tools found.</li>';
            return;
        }
        matches.forEach(function (item) {
            var li = document.createElement('li');
            var btn = document.createElement('button');
            btn.type = 'button';
            btn.className = 'w-full text-left px-4 py-3 hover:bg-teal-50/70 dark:hover:bg-teal-900/30 transition-colors flex flex-col gap-0.5';
            btn.innerHTML =
                '<span class="text-sm font-semibold text-slate-900 dark:text-slate-100">' + item.title + '</span>' +
                '<span class="text-xs text-slate-500 dark:text-slate-400">' + item.sub + '</span>';
            btn.addEventListener('click', function () {
                closePalette();
                if (item.action) item.action();
                else if (item.url) window.location.href = item.url;
            });
            li.appendChild(btn);
            list.appendChild(li);
        });
    }

    function openPalette() {
        if (!modal) buildModal();
        modal.style.display = 'flex';
        modal.classList.remove('hidden');
        input.value = '';
        renderItems('');
        setTimeout(function () { input.focus(); }, 20);
    }

    function closePalette() {
        if (modal) {
            modal.style.display = 'none';
            modal.classList.add('hidden');
        }
    }

    document.addEventListener('keydown', function (e) {
        if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === 'k') {
            e.preventDefault();
            if (modal && !modal.classList.contains('hidden')) closePalette();
            else openPalette();
        } else if (e.key === 'Escape' && modal && !modal.classList.contains('hidden')) {
            closePalette();
        }
    });

    window.openCommandPalette = openPalette;
})();
