/* calculator/static/js/lab-history.js */
(function () {
    const STORAGE_KEY = 'ilovechem_private_notebook_v1';
    const MAX_ENTRIES = 20;

    function getEntries() {
        try {
            return JSON.parse(localStorage.getItem(STORAGE_KEY) || '[]');
        } catch (e) {
            return [];
        }
    }

    function saveEntry(entry) {
        if (!entry || typeof entry !== 'object') return;
        const list = getEntries().filter(item => item.signature !== entry.signature);
        list.unshift(Object.assign({ ts: new Date().toISOString() }, entry));
        try {
            localStorage.setItem(STORAGE_KEY, JSON.stringify(list.slice(0, MAX_ENTRIES)));
        } catch (e) {
            // Ignore quota errors in private browsing mode
        }
    }

    function clearEntries() {
        try {
            localStorage.removeItem(STORAGE_KEY);
        } catch (e) {}
    }

    window.iLoveChemNotebook = {
        getEntries: getEntries,
        saveEntry: saveEntry,
        clearEntries: clearEntries
    };
})();
