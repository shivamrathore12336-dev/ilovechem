/*
 * Several cards are <div onclick="...">. Make them reachable and operable by keyboard
 * without rewriting every card: add role/tabindex and map Enter/Space to click.
 */
(function () {
    function enhance(root) {
        root.querySelectorAll('div[onclick], span[onclick], li[onclick]').forEach(function (el) {
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
    function init() { enhance(document); if (window.lucide) lucide.createIcons(); }
    if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', init);
    else init();
    window.enhanceClickables = enhance;
})();
